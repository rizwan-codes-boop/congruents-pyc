/**
 * @file solver.c
 * @brief Synchronous galaxy calculations with OpenMP dispatch.
 *
 * Calculates galaxy properties, transport, steady-state electron populations,
 * emission, internal absorption and diagnostics. Buffers crossing the API
 * belong to Python. Galaxies run in parallel; energy loops remain sequential
 * within each worker. No Python callbacks, file I/O, lookup-table generation
 * or nested OpenMP teams occur here.
 */
#include "solver.h"
#include "solver_runtime.h"
#include <stdatomic.h>
#ifndef CG_SERIAL_AUDIT
#include <omp.h>
#else
/* Select worker zero in the serial-only sanitizer build. */
static int omp_get_thread_num(void) { return 0; }
#endif
static double q_p_inject=2.2;
static double T_CR_lims__GeV[2]={1e-3,1e8};
static double E_CRe_lims__GeV[2]={1e-3+5.10998950e-4,1e8+5.10998950e-4};
/* Build the native logarithmic grid with explicit fused rounding. */
static int logspace_array(size_t n,double lo,double hi,double *out) {
    double a=log(lo),step=(log(hi)-a)/(n-1);
    /* Explicit single rounding pins the reference solver-bin convention on
     * compilers/platforms with different implicit contraction defaults. */
    for(size_t i=0;i<n;i++) out[i]=exp(fma(step,(double)i,a));
    return 0;
}
#define hcubature_v cg_cubature
#define malloc cg_allocate
#define free cg_release
#define free2D cg_free2D
#include "physics/worker_physics.h"
#undef free2D
#undef free
#undef malloc
#undef hcubature_v

static atomic_flag busy=ATOMIC_FLAG_INIT;
/* Return the interface version used to validate Python declarations. */
CG_API unsigned cg_solver_abi(void) { return 5; }
/* Report whether this library was compiled with active OpenMP workers. */
CG_API int cg_solver_openmp_enabled(void) {
#ifdef CG_SERIAL_AUDIT
    return 0;
#else
    return 1;
#endif
}
/* Validate a buffer of finite nonnegative or strictly positive values. */
static int values(size_t n,const double *v,int positive) {
    if(!v) return 0;
    for(size_t i=0;i<n;i++) if(!isfinite(v[i]) || (positive?v[i]<=0:v[i]<0)) return 0;
    return 1;
}
/* Require a bounded-size, positive and strictly increasing numerical axis. */
static int axis(size_t n,const double *v) {
    if(n<2 || n>4096 || !values(n,v,1)) return 0;
    for(size_t i=1;i<n;i++) if(v[i]<=v[i-1]) return 0;
    return 1;
}
/* Validate table dimensions, axes and nonnegative kernel values before dispatch. */
static int table_valid(const cg_table_input *t,int one) {
    return t && axis(t->nx,t->x) && (one?t->ny==1:axis(t->ny,t->y)) &&
        t->ny<=4096 && t->nx*t->ny<=4194304 && values(t->nx*t->ny,t->values,0);
}
/* Serialize this private GSL error-handler scope, not the workers within it. */
typedef void (*galaxy_job)(size_t,void *);
/* Run one galaxy job with a private error boundary and guaranteed allocation cleanup. */
static void execute_job(cg_worker *w,size_t i,galaxy_job job,void *data) {
    w->count=0;w->status=0;w->phase="galaxy";cg_active=w;
    if(!setjmp(w->escape)) job(i,data);
    cg_cleanup(w);cg_active=NULL;
}
/* Dispatch disjoint galaxy rows over OpenMP workers under a private GSL error scope. */
static int batch(int threads,size_t n,int *statuses,galaxy_job job,void *data) {
    if(threads<1 || !n || n>100000 || !statuses) return CG_INVALID;
    if(atomic_flag_test_and_set(&busy)) return CG_INVALID;
    gsl_error_handler_t *previous=gsl_set_error_handler_off();
    int team=(size_t)threads<n?threads:(int)n;
    cg_worker *workers=calloc(team,sizeof(*workers));
    if(!workers) { gsl_set_error_handler(previous);atomic_flag_clear(&busy);return CG_ALLOC; }
    #pragma omp parallel for num_threads(team) schedule(guided)
    for(size_t i=0;i<n;i++) {
        cg_worker *w=workers+omp_get_thread_num();
        execute_job(w,i,job,data);
        statuses[i]=w->status;
    }
    free(workers);gsl_set_error_handler(previous);atomic_flag_clear(&busy);
    return CG_OK;
}

/* Reference first region plus the property-only part of its second region.
 * Catalogue columns: z, Mstar Msun, Re kpc, SFR Msun/yr.
 * Property columns match Python's PROPERTY_NAMES exactly. */
typedef struct { const double *rows; double *out; } property_job;
/* Derive gas density, magnetic fields, geometry and dust temperature for one galaxy. */
static void properties_one(size_t i,void *data) {
    property_job *job=data;const double *r=job->rows+4*i;double *p=job->out+10*i;
    const double area=M_PI*pow(r[2]*1e3,2),stars=r[1]/(2*area),sfr=r[3]/(2*area);
    const double gas=Sigma_gas_Shi_iKS__Msolpcm2(sfr,stars),sigma=sigma_gas_Yu__kmsm1(r[3]);
    const double h=pow(sigma,2)/(M_PI*4.302e-3*(gas+sigma/sigma_star_Bezanson__kmsm1(r[1],r[2])*stars));
    const double nh=gas/(1.4*m_H__kg*2*h)*Msol__kg/pow(pc__cm,3);
    const double vai=1000*((sigma/sqrt(2.))/10.)/2.;
    const double b=sqrt(4*M_PI*1e-4*(nh*1.4/1.17)*1.17*m_H__kg*1e3)*vai*1e5;
    double v[10]={h,nh,b,sigma,area,gas,sfr,stars,Tdust__K(r[0],r[3],r[1]),
                  log10(r[3]/r[1])>-10?b/3:b/1.5};
    if(!values(10,v,1)) cg_fail(CG_NUMERIC);
    for(size_t k=0;k<10;k++) p[k]=v[k];
}
/* Validate catalogue buffers and dispatch the parallel galaxy-property calculation. */
CG_API int cg_properties(int threads,size_t n,const double *rows,double *out,int *status) {
    if(!n || n>100000 || !rows || !out) return CG_INVALID;
    for(size_t i=0;i<n;i++) if(!isfinite(rows[4*i]) || rows[4*i]<0 || rows[4*i]>20 || !values(3,rows+4*i+1,1)) return CG_INVALID;
    property_job job={rows,out};return batch(threads,n,status,properties_one,&job);
}

/* Reference transport region. The primary and proton arrays are computed here
 * as part of each galaxy job; secondary injection stays in the spectrum job. */
typedef struct { size_t ne; const double *rows,*props,*kinetic;double *out,*cp; } transport_job;
/* Calculate the energy-dependent streaming speed, capped at the speed of light. */
static double streaming(double t,double mass,double nh,double ion,double norm,double vai,double ula) {
    return fmin(vai*(1+2.3e-3*pow(sqrt(t*t+2*mass*t),1.2)*pow(nh/1e3,1.5)*ion*2/(ula/10*norm/2e-7)),c__cmsm1/1e5);
}
/* Compute one galaxy's injection normalization, calorimetry and diffusion arrays. */
static void transport_one(size_t i,void *data) {
    transport_job *job=data;size_t ne=job->ne;
    const double *r=job->rows+4*i,*p=job->props+10*i,*t=job->kinetic;
    double *o=job->out+7*ne*i;
    double norm=C_norm_E(2.2,m_p__GeV,1e8),norme=C_norm_E(2.2,m_e__GeV,1e5);
    double power=r[3]*1.321680e-2*.1*1e51*erg__GeV/yr__s;
    double ula=p[3]/sqrt(2.),vai=1000*(ula/10)/2,la=p[0]/8.;
    double collision=1/(p[1]*1.4/1.17*40e-27*.5*c__cmsm1),d0=vai*la*1e5*pc__cm;
    double loss=1/(1/collision+1/(pow(p[0]*pc__cm,2)/d0));
    double cn=power*loss/(norm*2*p[4]*2*p[0]*pow(pc__cm,3));
    double cp=power*collision/norm;job->cp[i]=cp;
    for(size_t j=0;j<ne;j++) {
        double vs=streaming(t[j],m_p__GeV,p[1],1,cn,vai,ula);
        o[ne+j]=vs*la*1e5*pc__cm;
        double tau=9.9*p[5]/1e3*p[0]/1e2*1e27/o[ne+j];
        double gam=41.2*p[0]/1e2*vs/1e3*1e27/o[ne+j];
        o[j]=1-1/(gsl_sf_hyperg_0F1(.25/1.25,tau/pow(1.25,2))+tau/gam*gsl_sf_hyperg_0F1(2.25/1.25,tau/pow(1.25,2)));
        if(i==10) o[j]*=.1;
    }
    for(size_t j=0;j<ne;j++) {
        o[2*ne+j]=streaming(t[j],m_e__GeV,p[1],1,cn,vai,ula)*la*1e5*pc__cm;
        o[3*ne+j]=streaming(t[j],m_e__GeV,p[1]/1e3,1e4,(1-o[0])*cn,vai,ula)*la*1e5*pc__cm;
        o[4*ne+j]=J(t[j],.2*power/norme,2.2,m_e__GeV,1e5);
        o[5*ne+j]=0.;
        o[6*ne+j]=J(t[j],cp,2.2,m_p__GeV,1e8)*o[j];
    }
    if(!values(7*ne,o,0)||!isfinite(cp)||cp<=0) cg_fail(CG_NUMERIC);
}
/* Validate transport inputs and dispatch galaxy jobs into Python-owned output buffers. */
CG_API int cg_transport(int threads,size_t n,size_t ne,const double *rows,const double *props,const double *t,double *out,double *cp,int *status) {
    if(!n||n>100000||!axis(ne,t)||!values(n*4,rows,0)||!values(n*10,props,1)||!out||!cp) return CG_INVALID;
    transport_job job={ne,rows,props,t,out,cp};return batch(threads,n,status,transport_one,&job);
}

/* Combine reusable field planes inside the galaxy worker. The family order is
 * 3000/4000/7500/UV, all CMB temperatures, then all FIR temperatures.
 * Preserve reversed reference temperature weights and emission rounding. */
/* Combine field kernels using galaxy dilution and the retained temperature weights. */
static gsl_spline_object_2D combine(const cg_run *a,const cg_table_input *tables,const double *r,const double *p,int rounded) {
    double temps[2]={T_0_CMB__K*(1+r[0]),p[8]};
    size_t starts[2]={4,4+a->ncmb},counts[2]={a->ncmb,a->nfir},idx[2];double w[2];
    for(size_t k=0;k<2;k++) {
        size_t j=0;while(j<counts[k] && tables[starts[k]+j].temperature<=temps[k]) j++;
        if(j==0||j==counts[k]) cg_fail(CG_INVALID);
        idx[k]=starts[k]+j-1;
        double lo=tables[idx[k]].temperature,hi=tables[idx[k]+1].temperature;
        double index=(double)(j-1);
        w[k]=(index+(temps[k]-lo)/(hi-lo))-index;
    }
    (void)rounded;
    double d[5]={C_dil(u_rad_BB__GeVcmm3(3000),L3000K__Lsol(r[1]),r[2],p[0]),
      C_dil(u_rad_BB__GeVcmm3(4000),L4000K__Lsol(r[1]),r[2],p[0]),
      C_dil(u_rad_BB__GeVcmm3(7500),L7500K__Lsol(r[3]),r[2],p[0]),
      C_dil(u_rad_UVMattis__GeVcmm3(),LUV__Lsol(r[3]),r[2],p[0]),
      C_dil(u_rad_modBB__GeVcmm3(p[8]),LFIR__Lsol(r[3]),r[2],p[0])};
    size_t len=tables[0].nx*tables[0].ny;double *z=cg_allocate(len*sizeof(double));
    for(size_t j=0;j<len;j++) {
        z[j]=d[0]*tables[0].values[j]+d[1]*tables[1].values[j]+d[2]*tables[2].values[j]+d[3]*tables[3].values[j];
        z[j]+=w[0]*tables[idx[0]].values[j]+(1-w[0])*tables[idx[0]+1].values[j];
        z[j]+=d[4]*(w[1]*tables[idx[1]].values[j]+(1-w[1])*tables[idx[1]+1].values[j]);
    }
    gsl_spline_object_2D result=gsl_so2D(tables[0].nx,tables[0].ny,tables[0].x,tables[0].y,z);
    cg_release(z);return result;
}

/* Energy-weighted log integration of a linearly interpolated exported array.
 * Separate from the serial observer luminosity calculation in Python. */
typedef struct { gsl_spline_object_1D spectrum; } energy_data;
/* Evaluate E-squared times the interpolated spectrum for log-energy integration. */
static int energy_integrand(unsigned dim,size_t n,const double *x,void *data,unsigned fdim,double *v) {
    energy_data *d=data;for(size_t j=0;j<n;j++) {double e=exp(x[j]);v[j]=e*e*gsl_so1D_eval(d->spectrum,e);}return 0;
}
/* Integrate energy-weighted spectral power over the supplied energy grid. */
static double energy_integral(size_t n,const double *x,const double *y) {
    energy_data data={gsl_so1D(n,x,y)};double lo=log(x[0]),hi=log(x[n-1]),v,err;
    cg_cubature(1,energy_integrand,&data,1,&lo,&hi,100000,0,1e-6,ERROR_INDIVIDUAL,&v,&err);
    gsl_so1D_free(data.spectrum);return v;
}

/* Solve disc then halo populations, compute emission and absorption, and fill diagnostics. */
static void spectra_one(size_t i,void *data) {
    const cg_run *a=data;size_t ne=a->ne,np=a->np;
    const double *r=a->rows+4*i,*p=a->props+10*i,*e=a->electron,*ph=a->photon,*t=a->kinetic;
    double *tr=a->transport+7*ne*i,*el=a->electrons+4*ne*i,*em=a->emission+15*np*i;
    gsl_spline_object_2D ic=combine(a,a->ic,r,p,1),gamma=combine(a,a->gamma,r,p,0);
    gsl_spline_object_2D bs=gsl_so2D(a->bs->nx,a->bs->ny,a->bs->x,a->bs->y,a->bs->values);
    gsl_spline_object_1D sy=gsl_so1D(a->sy->nx,a->sy->x,a->sy->values),fc=gsl_so1D(ne,t,tr);
    for(size_t j=0;j<ne;j++) tr[5*ne+j]=q_e(t[j],p[1],a->cp[i],1e8,fc);
    gsl_spline_object_1D dd=gsl_so1D(ne,e,tr+2*ne),dh=gsl_so1D(ne,e,tr+3*ne);
    gsl_spline_object_1D q1=gsl_so1D(ne,e,tr+4*ne),q2=gsl_so1D(ne,e,tr+5*ne),ss[4];
    double bounds[2]={E_CRe_lims__GeV[0],E_CRe_lims__GeV[1]};
    cg_active->phase="disc matrix and solve";
    CRe_steadystate_solve(1,bounds,(int)a->ns,p[1],p[2],p[0],1,&gamma,bs,dd,q1,q2,ss,ss+1);
    double *escape=cg_allocate(2*ne*sizeof(double));
    for(size_t j=0;j<ne;j++) for(size_t k=0;k<2;k++) escape[k*ne+j]=gsl_so1D_eval(ss[k],e[j])/tau_diff__s(e[j],p[0],dd);
    gsl_spline_object_1D qh1=gsl_so1D(ne,e,escape),qh2=gsl_so1D(ne,e,escape+ne);
    cg_active->phase="halo matrix and solve";
    CRe_steadystate_solve(2,bounds,(int)a->ns,p[1]/1000,p[9],50*p[0],1,&gamma,bs,dh,qh1,qh2,ss+2,ss+3);
    for(size_t k=0;k<4;k++) for(size_t j=0;j<ne;j++) el[k*ne+j]=gsl_so1D_eval(ss[k],e[j]);
    double radiation[6]={T_0_CMB__K*(1+r[0]),p[8],r[1],r[3],r[2],p[0]};
    double limits[2]={a->target_low,a->target_high};
    cg_active->phase="emission and internal absorption";
    for(size_t j=0;j<np;j++) {
        double energy=ph[j],tau=tau_FF_MK(energy,p[6],1e4),atten=exp(-tau);
        em[j]=eps_IC_3(energy,ic,ss[0])*atten;em[np+j]=eps_IC_3(energy,ic,ss[1])*atten;
        em[2*np+j]=eps_BS_3(energy,p[1],bs,ss[0])*atten;em[3*np+j]=eps_BS_3(energy,p[1],bs,ss[1])*atten;
        em[4*np+j]=eps_SY_4(energy,p[2],sy,ss[0])*atten;em[5*np+j]=eps_SY_4(energy,p[2],sy,ss[1])*atten;
        em[6*np+j]=eps_IC_3(energy,ic,ss[2]);em[7*np+j]=eps_IC_3(energy,ic,ss[3]);
        em[8*np+j]=eps_SY_4(energy,p[9],sy,ss[2]);em[9*np+j]=eps_SY_4(energy,p[9],sy,ss[3]);
        em[10*np+j]=eps_FF(energy,r[2],1e4,tau);em[11*np+j]=tau;
        em[12*np+j]=eps_pi(energy,p[1],a->cp[i],1e8,fc);
        em[13*np+j]=eps_pi_fcal1(energy,p[1],a->cp[i],1e8,fc);
        em[14*np+j]=q_nu(energy,p[1],a->cp[i],1e8,fc);
        a->internal[i*np+j]=tau_gg_gal_BW(energy,dndEphot_total__cmm3GeVm1,radiation,limits,p[0]);
    }
    if(!values(7*ne,tr,0)||!values(4*ne,el,0)||!values(15*np,em,0)||!values(np,a->internal+i*np,0)) cg_fail(CG_NUMERIC);
    if(!a->diagnostics) return;
    cg_active->phase="diagnostics";
    double *loss=a->loss+12*ne*i,*critical=a->critical+10*i,*budget=a->budget+16*i,*radio=a->radio+5*i;
    for(size_t z=0;z<2;z++) {
        double nh=z?p[1]/1000:p[1],b=z?p[9]:p[2],h=z?50*p[0]:p[0];
        gsl_spline_object_1D diffusion=z?dh:dd;
        for(size_t j=0;j<ne;j++) {
            loss[(5*z)*ne+j]=tau_sync__s(e[j],b);
            /* At the exact lower bound the diagnostic integral has zero width. */
            loss[(5*z+1)*ne+j]=e[j]<=bounds[0]?INFINITY:tau_BS_fulltest__s(e[j],bounds,nh,bs);
            loss[(5*z+2)*ne+j]=e[j]<=bounds[0]?INFINITY:tau_IC_fulltest__s(e[j],bounds,gamma);
            loss[(5*z+3)*ne+j]=tau_diff__s(e[j],h,diffusion);
            loss[(5*z+4)*ne+j]=tau_ion__s(e[j],nh);
        }
        double ec=sqrt(2*1.49e9*m_e__g*c__cmsm1/(3*b*e__esu))*M_PI*m_e__GeV;
        critical[5*z]=ec/tau_BS_fulltest__s(ec,bounds,nh,bs);
        critical[5*z+1]=ec/tau_sync__s(ec,b);
        critical[5*z+2]=ec/tau_IC_fulltest__s(ec,bounds,gamma);
        critical[5*z+3]=ec/(z?tau_plasma__s(ec,nh):tau_ion__s(ec,nh));
        critical[5*z+4]=ec/tau_diff__s(ec,h,diffusion);
    }
    for(size_t j=0;j<ne;j++) {
        loss[10*ne+j]=1/(p[1]*1.4/1.17*40e-27*.5*c__cmsm1);
        loss[11*ne+j]=pow(p[0]*pc__cm,2)/tr[ne+j];
    }
    for(size_t j=0;j<2*ne;j++) a->escape[2*ne*i+j]=escape[j];
    for(size_t k=0;k<2;k++) {
        budget[k]=energy_integral(ne,t,tr+(4+k)*ne);
        size_t start=k?9:2,indices[5]={4+k,k,2+k,8+k,6+k};
        for(size_t j=0;j<5;j++) budget[start+j]=energy_integral(np,ph,em+indices[j]*np)/budget[k];
        budget[start+5]=0;budget[start+6]=energy_integral(ne,t,escape+k*ne)/budget[k];
    }
    double re=1.49e9*h__GeVs,rtau=tau_FF_MK(re,p[6],1e4),factor=1.49e9*h__Js*h__GeVs;
    for(size_t k=0;k<4;k++) radio[k]=eps_SY_4(re,k<2?p[2]:p[9],sy,ss[k])*(k<2?exp(-rtau):1)*factor;
    radio[4]=eps_FF(re,r[2],1e4,rtau)*factor;
    for(size_t j=0;j<12*ne;j++) if(isnan(loss[j])||loss[j]<0) cg_fail(CG_NUMERIC);
    if(!values(10,critical,0)||!values(16,budget,0)||!values(5,radio,0)||!values(2*ne,escape,0)) cg_fail(CG_NUMERIC);
}

/* Validate shared tables and run complete independent galaxy-spectrum workers. */
CG_API int cg_spectra(int threads,const cg_run *a,int *status) {
    if(!a||!a->n||a->n>100000||a->ns<4||a->ns>500||!axis(a->ne,a->electron)||!axis(a->ne,a->kinetic)||!axis(a->np,a->photon)||
       !values(4*a->n,a->rows,0)||!values(10*a->n,a->props,1)||!values(a->n,a->cp,1)||!values(7*a->n*a->ne,a->transport,0)||
       !a->electrons||!a->emission||!a->internal||!a->ic||!a->gamma||a->ncmb<2||a->nfir<2||a->ncmb>4096||a->nfir>4096||
       !isfinite(a->target_low)||!isfinite(a->target_high)||a->target_low<=0||a->target_high<=a->target_low||
       !table_valid(a->bs,0)||!table_valid(a->sy,1)) return CG_INVALID;
    if(a->diagnostics && (!a->loss||!a->critical||!a->budget||!a->radio||!a->escape)) return CG_INVALID;
    for(size_t family=0;family<2;family++) {
        const cg_table_input *ts=family?a->gamma:a->ic;
        for(size_t j=0;j<4+a->ncmb+a->nfir;j++) {
            if(!table_valid(ts+j,0)||ts[j].nx!=ts[0].nx||ts[j].ny!=ts[0].ny) return CG_INVALID;
            for(size_t k=0;k<ts[0].nx;k++) if(ts[j].x[k]!=ts[0].x[k]) return CG_INVALID;
            for(size_t k=0;k<ts[0].ny;k++) if(ts[j].y[k]!=ts[0].y[k]) return CG_INVALID;
            if(j>=4 && (!isfinite(ts[j].temperature)||ts[j].temperature<=0)) return CG_INVALID;
            if((j>4 && j<4+a->ncmb) || j>4+a->ncmb) if(ts[j].temperature<=ts[j-1].temperature) return CG_INVALID;
        }
    }
    for(size_t k=0;k<2;k++) if(fabs(a->electron[k?a->ne-1:0]/E_CRe_lims__GeV[k]-1)>1e-12) return CG_INVALID;
    return batch(threads,a->n,status,spectra_one,(void *)a);
}
