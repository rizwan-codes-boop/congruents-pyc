/**
 * @file steady_state_native.h
 * @brief Finite-volume assembly and steady-state electron solves.
 *
 * Defines transition, injection, escape and continuous-loss integrals used
 * to assemble the disc/halo kinetic equations for two electron populations.
 * Uses checked worker-owned allocations and the runtime linear solver;
 * reference energy-grid rounding is retained.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */

/* Exact integrals of the bilinear rate table at fixed initial energy.
 * Split at every transfer-energy knot. Within a segment P is linear, so its
 * outgoing, incoming and lost-energy moments are analytic positive sums.
 * This avoids both an unresolved Thomson peak and cancellation of Ei-Ef. */
static void table_transition_moments(gsl_spline_object_2D table,double initial,
                                    double low,double high,double *outgoing,double *incoming,double *loss) {
    gsl_spline2d *s=table.spline;
    size_t nx=s->interp_object.xsize,ny=s->interp_object.ysize;
    if(initial<s->yarr[0] || initial>s->yarr[ny-1]) return;
    low=fmax(low,s->xarr[0]);high=fmin(high,s->xarr[nx-1]);
    if(high<=low) return;
    size_t iy=gsl_interp_bsearch(s->yarr,initial,0,ny-1);
    double wy=(initial-s->yarr[iy])/(s->yarr[iy+1]-s->yarr[iy]);
    size_t ix=gsl_interp_bsearch(s->xarr,low,0,nx-1);
    long double o=0.,in=0.,lost=0.;
    while(low<high && ix+1<nx) {
        double right=fmin(high,s->xarr[ix+1]),width=right-low;
        double zl=(1-wy)*s->zarr[iy*nx+ix]+wy*s->zarr[(iy+1)*nx+ix];
        double zr=(1-wy)*s->zarr[iy*nx+ix+1]+wy*s->zarr[(iy+1)*nx+ix+1];
        double slope=(zr-zl)/(s->xarr[ix+1]-s->xarr[ix]);
        double pl=fmax(0.,zl+slope*(low-s->xarr[ix]));
        double pr=fmax(0.,zl+slope*(right-s->xarr[ix]));
        o+=(long double)initial*width*.5*(pl+pr);
        lost+=(long double)width/6.*((2*low+right)*pl+(low+2*right)*pr);
        double fl=initial-low,fr=initial-right;
        in+=(long double)width/6.*((2*fl+fr)*pl+(fl+2*fr)*pr);
        low=right;ix++;
    }
    *outgoing+=(double)o;*incoming+=(double)in;*loss+=(double)lost;
}

/* Only the initial log energy remains numerical. The transfer moment is exact
 * for the supplied table, including its zero-outside-domain convention. */
static int transition_moments(unsigned ndim,size_t npts,const double *x,void *data,unsigned fdim,double *v) {
    struct F_int_data *f=data;
    for(size_t j=0;j<npts;j++) {
        double initial=exp(x[j*ndim]),low=fmax(0.,initial-f->final_high);
        double high=initial-fmax(m_e__GeV,f->final_low);
        double outgoing=0.,incoming=0.,lost=0.;
        table_transition_moments(f->gso2D_BS,initial,low,high,&outgoing,&incoming,&lost);
        double factor=c__cmsm1*f->n_H__cmm3*mb__cm2;
        outgoing*=factor;incoming*=factor;lost*=factor;
        for(int k=0;k<f->n_gso2D;k++)
            table_transition_moments(f->gso_2D_radfield[k],initial,low,high,&outgoing,&incoming,&lost);
        v[j*fdim]=f->same_cell?lost:outgoing;
        v[j*fdim+1]=x[j*ndim]*v[j*fdim];
        if(fdim==4) {
            v[j*fdim+2]=f->same_cell?0.:incoming;
            v[j*fdim+3]=x[j*ndim]*v[j*fdim+2];
        }
    }
    return 0;
}

/* Two outgoing moments for transitions below the tracked energy range. */
static int F_Gamma_i0_2D_2_log_E(unsigned ndim,size_t n,const double *x,void *p,unsigned fdim,double *v) {
    return transition_moments(ndim,n,x,p,fdim,v);
}

/* Outgoing/incoming energy moments, or their stable same-cell difference. */
static int F_Gamma_2D_4_log_E(unsigned ndim,size_t n,const double *x,void *p,unsigned fdim,double *v) {
    return transition_moments(ndim,n,x,p,fdim,v);
}

/* Evaluate escape and continuous-loss moments for the electron matrix. */
static int F_EdotDE_2_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j * fdim + 0] = exp(x[j*ndim+0]) * gsl_so1D_eval( fdata_in.gso_1D_D__cm2sm1, exp(x[j*ndim+0]) )/pow(fdata_in.h__pc*pc__cm,2) -
                             fdata_in.E_func( exp(x[j*ndim+0]), fdata_in.B__G, fdata_in.n_H__cmm3, fdata_in.h__pc );
        fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
    }
    return 0;
}

/* Evaluate the energy-squared injection term in logarithmic electron-energy coordinates. */
static int F_QE2_i_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j] = pow(exp(x[j*ndim+0]),2) * gsl_so1D_eval( fdata_in.gso_1D_Q, exp(x[j*ndim+0]) );
    }
    return 0;
}

/* Add a coefficient multiplying slope_j*Delta(log E), using B5. */
static void add_slope(double *row,int j,int n,double coefficient) {
    if(j==0) { row[0]-=coefficient;row[1]+=coefficient; }
    else if(j==n-1) { row[n-2]-=coefficient;row[n-1]+=coefficient; }
    else { row[j-1]-=.5*coefficient;row[j+1]+=.5*coefficient; }
}

/* Solve only the contiguous supported range, removing TRAILING zero sources.
 * A wholly zero source must return zero without a singular or zero-size solve. */
static void solve_population(int n,double **matrix,double *source,const double *energy,double *out) {
    int active=n;
    while(active>0 && source[active-1]==0.) active--;
    for(int i=0;i<n+2;i++) out[i]=0.;
    if(!active) return;
    double *a=malloc((size_t)active*active*sizeof(double));
    double *x=malloc(active*sizeof(double));
    for(int i=0;i<active;i++) for(int j=0;j<active;j++) a[i*active+j]=matrix[i][j];
    cg_active->population_energy=energy+1;
    cg_linear_solve(active,a,source,x);
    cg_active->population_energy=NULL;
    for(int i=0;i<active;i++) out[i+1]=x[i]/energy[i+1];
    /* Avoid log(0) and extrapolated infinities after a flagged negative cell.
     * Use constant continuation only in that unresolved boundary case. */
    out[0]=out[1];
    if(active>1 && out[1]>0. && out[2]>0.)
        out[0]=exp(log(out[1])+(log(out[2])-log(out[1]))*
                    log(energy[0]/energy[1])/log(energy[2]/energy[1]));
    /* No incoming population at the upper boundary. */
    out[n+1]=0.;
    free(x);free(a);
}

/* Assemble and solve the coupled primary/secondary steady-state electron equations. */
static int CRe_steadystate_solve( int structure, double E_e_lim__GeV[2], int n_E, double n_H__cmm3, double B__G, double h__pc,
    unsigned int n_gso2D, gsl_spline_object_2D * gso_2D_radfields, gsl_spline_object_2D gso2D_BS, gsl_spline_object_1D gso_1D_D__cm2sm1,
    gsl_spline_object_1D gso_1D_Q_inject_1, gsl_spline_object_1D gso_1D_Q_inject_2,
    gsl_spline_object_1D * qe_1_so_1D, gsl_spline_object_1D * qe_2_so_1D )
{

    if(n_E<4 || n_E>500) cg_fail(CG_INVALID);
    double E__GeV[n_E+1];
    logspace_array( n_E+1, E_e_lim__GeV[0], E_e_lim__GeV[1], E__GeV );

    double DeltalogE = log(E_e_lim__GeV[1]/E_e_lim__GeV[0])/n_E;

    int i,j;

    double lnE_i__GeV[n_E];
    for (i = 0; i < n_E; ++i)
    {
        lnE_i__GeV[i] = (log(E__GeV[i])+log(E__GeV[i+1]))/2.;
    }

    double E_out__GeV[n_E+2];
    for (i = 0; i < n_E; ++i)
    {
        E_out__GeV[i+1] = exp((log(E__GeV[i])+log(E__GeV[i+1]))/2.);
    }
    E_out__GeV[0] = E__GeV[0];
    E_out__GeV[n_E+1] = E__GeV[n_E];

    double xmin[1], xmax[1];
    double xmin2D[2], xmax2D[2];
    double res, abserr;
    double res2[2], abserr2[2];
    double res4[4], abserr4[4];

    struct F_int_data fdata;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.B__G = B__G;
    fdata.h__pc = h__pc;
    fdata.n_gso2D = n_gso2D;
    fdata.gso_2D_radfield = gso_2D_radfields;
    fdata.gso2D_BS = gso2D_BS;
    fdata.gso_1D_D__cm2sm1 = gso_1D_D__cm2sm1;

    double Gamma_i0[n_E];
    double Gamma_i0_prime[n_E];
    xmin2D[1] = 0.;xmax2D[1] = 1.;
    fdata.final_low=m_e__GeV;fdata.final_high=E__GeV[0];fdata.same_cell=0;
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_Gamma_i0_2D_2_log_E, &fdata, 1, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        Gamma_i0[i] = res2[0]/DeltalogE;
        Gamma_i0_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_i0[i];
    }

    double Gamma_ii[n_E];
    double Gamma_ii_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        fdata.final_low=E__GeV[i];fdata.final_high=E__GeV[i+1];fdata.same_cell=1;
        hcubature_v( 4, F_Gamma_2D_4_log_E, &fdata, 1, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res4, abserr4 );
        Gamma_ii[i] = (res4[0]/DeltalogE) - (res4[2]/DeltalogE) ;
        Gamma_ii_prime[i] = (res4[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * (res4[0]/DeltalogE)) -
                            (res4[3]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * (res4[2]/DeltalogE));
    }

    double **Gamma_ij = malloc(sizeof *Gamma_ij * n_E);
    if (Gamma_ij){for (i = 0; i < n_E; i++){Gamma_ij[i] = malloc(sizeof *Gamma_ij[i] * n_E);}}
    double **Gamma_ij_prime = malloc(sizeof *Gamma_ij_prime * n_E);
    if (Gamma_ij_prime){for (i = 0; i < n_E; i++){Gamma_ij_prime[i] = malloc(sizeof *Gamma_ij_prime[i] * n_E);}}
    double **Gamma_ji_prime = malloc(sizeof *Gamma_ji_prime * n_E);
    if (Gamma_ji_prime){for (i = 0; i < n_E; i++){Gamma_ji_prime[i] = malloc(sizeof *Gamma_ji_prime[i] * n_E);}}
    double **Gamma_ji = malloc(sizeof *Gamma_ji * n_E);
    if (Gamma_ji){for (i = 0; i < n_E; i++){Gamma_ji[i] = malloc(sizeof *Gamma_ji[i] * n_E);}}

    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        for (j = 0; j < i; ++j)
        {
            fdata.final_low=E__GeV[j];fdata.final_high=E__GeV[j+1];fdata.same_cell=0;

            hcubature_v( 4, F_Gamma_2D_4_log_E, &fdata, 1, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res4, abserr4 );

            Gamma_ij[i][j] = res4[0]/DeltalogE;
            Gamma_ij_prime[i][j] = res4[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_ij[i][j];
            Gamma_ji[i][j] = res4[2]/DeltalogE;
            Gamma_ji_prime[i][j] = res4[3]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_ji[i][j];
        }
        for (j = i; j < n_E; ++j)
        {
            Gamma_ij[i][j] = 0.;
            Gamma_ij_prime[i][j] = 0.;
            Gamma_ji[i][j] = 0.;
            Gamma_ji_prime[i][j] = 0.;
        }
    }

    double Gamma_i[n_E];
    double Gamma_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        Gamma_i[i] = Gamma_i0[i] + Gamma_ii[i];
        Gamma_i_prime[i] = Gamma_i0_prime[i] + Gamma_ii_prime[i];
        for (j = 0; j < i; ++j)
        {
            Gamma_i[i] += Gamma_ij[i][j];
            Gamma_i_prime[i] += Gamma_ij_prime[i][j];
        }
    }

    double Edot_i[n_E+1];
    if (structure == 1)
    {
        fdata.E_func = dEdtm1_total_disc__GeVsm1;
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_disc__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc );
        }
    }
    else if (structure == 2)
    {
        fdata.E_func = dEdtm1_total_halo__GeVsm1;
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_halo__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc );
        }
    }
    else
    {
        printf("No valid structure specified in CRe steady state solver!");
    }

    double D_i[n_E];
    double D_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_EdotDE_2_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        D_i[i] = res2[0]/DeltalogE;
        D_i_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * D_i[i];
    }

    double **M = malloc(sizeof *M * n_E);
    for(i=0;i<n_E;i++) {
        M[i]=malloc(n_E*sizeof **M);
        for(j=0;j<n_E;j++) M[i][j]=0.;
    }

    /* Assemble directly from B1-B5, with the LOWER face reconstructed as
     * alpha_i - slope_i*Delta/2. The plus sign printed in B6/B17 and used by
     * the old stencil evaluates the wrong face of the upwind cell.
     * At the top edge there is no incoming population from E > Emax.
     * The bottom edge is free outflow, not a reflecting boundary. */
    for(i=0;i<n_E;i++) {
        M[i][i] -= D_i[i]+Gamma_i[i];
        M[i][i] += Edot_i[i]/DeltalogE;
        add_slope(M[i],i,n_E,-(D_i_prime[i]+Gamma_i_prime[i])
                  -Edot_i[i]/(2.*DeltalogE));
        if(i+1<n_E) {
            M[i][i+1] -= Edot_i[i+1]/DeltalogE;
            add_slope(M[i],i+1,n_E,Edot_i[i+1]/(2.*DeltalogE));
        }
        for(j=i+1;j<n_E;j++) {
            M[i][j] += Gamma_ji[j][i];
            add_slope(M[i],j,n_E,Gamma_ji_prime[j][i]);
        }
    }

    double Q_i_1[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_1;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_1[i] = -1.*res/DeltalogE;
    }

    double Q_i_2[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_2;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_2[i] = -1.*res/DeltalogE;
    }

    double q_e_1[n_E+2],q_e_2[n_E+2];
    solve_population(n_E,M,Q_i_1,E_out__GeV,q_e_1);
    solve_population(n_E,M,Q_i_2,E_out__GeV,q_e_2);

    *qe_1_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_1 );
    *qe_2_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_2 );

    free2D( n_E, Gamma_ij );
    free2D( n_E, Gamma_ij_prime );
    free2D( n_E, Gamma_ji );
    free2D( n_E, Gamma_ji_prime );

    free2D( n_E, M );

    return 0;

}
