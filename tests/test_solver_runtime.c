/* Standalone checked-runtime audit: no Python interpreter or OpenMP needed. */
#define CG_TESTING
#include "../csrc/solver.c"
#include <string.h>
#include <assert.h>

/* Deliberately return invalid values to exercise the integration failure path. */
static int invalid_integrand(unsigned dim,size_t n,const double *x,void *p,unsigned fd,double *v) {
    (void)dim;(void)x;(void)p;
    for(size_t i=0;i<n*fd;i++) v[i]=NAN;
    return 0;
}
/* Exercise allocation, integration and linear-solve recovery under injected failures. */
static int exercise(size_t fail_at, int mode) {
    /* Heap storage remains defined across setjmp/longjmp. */
    cg_worker *w=calloc(1,sizeof(*w));
    assert(w);
    w->phase="runtime audit";w->fail_allocation=fail_at;cg_active=w;
    if(!setjmp(w->escape)) {
        double axis[2]={1.,2.}, values[2]={2.,4.}, plane[4]={2.,3.,3.,4.};
        gsl_spline_object_1D a=gsl_so1D(2,axis,values);
        gsl_spline_object_2D b=gsl_so2D(2,2,axis,axis,plane);
        assert(gsl_so1D_eval(a,1.5)==3.);
        assert(gsl_so2D_eval(b,1.5,1.5)==3.);
        assert(isnan(gsl_so1D_eval(a,.5)));
        assert(gsl_so2D_eval(b,.5,1.5)==0.);
        double matrix[4]={2.,0.,0.,4.}, rhs[2]={4.,8.}, solution[2];
        if(mode==1) matrix[3]=0.; /* Singular LU must not abort the process. */
        cg_linear_solve(2,matrix,rhs,solution);
        assert(solution[0]==2. && solution[1]==2.);
        void *extra=cg_allocate(128);
        cg_release(extra);
        gsl_so1D_free(a);gsl_so2D_free(b);
        if(mode==2) {
            double low=0.,high=1.,v,error;
            assert(cg_cubature(1,invalid_integrand,NULL,1,&low,&high,100,0.,1e-8,
                               ERROR_INDIVIDUAL,&v,&error)==1);
            assert(isnan(v));
        }
    }
    int status=w->status;
    cg_cleanup(w);
    assert(w->count==0);
    cg_cleanup(w); /* Cleanup is idempotent. */
    free(w);cg_active=NULL;
    return status;
}
/* Check the real worker dispatcher and table-combination convention using
 * constant synthetic planes, independent of scientific output fixtures. */
/* Check combined field tables and interpolation weights inside a worker boundary. */
static void combination_job(size_t i,void *data) {
    double axes[2]={1.,2.},zero[4]={0.,0.,0.,0.},low[4]={2.,2.,2.,2.},high[4]={10.,10.,10.,10.};
    cg_table_input tables[8];
    for(size_t j=0;j<8;j++) tables[j]=(cg_table_input){2,2,axes,axes,zero,0.};
    tables[4].values=low;tables[4].temperature=T_0_CMB__K;
    tables[5].values=high;tables[5].temperature=T_0_CMB__K+.5;
    tables[6].temperature=29.;tables[7].temperature=34.;
    cg_run args={0};args.ncmb=2;args.nfir=2;
    double row[4]={*(int *)data?1.:0.,1e10,2.,1.};
    double props[10]={100.,1.,1e-5,10.,1.,1.,1.,1.,30.,1e-6};
    gsl_spline_object_2D result=combine(&args,tables,row,props,0);
    /* At an exact lower temperature node the reference uses the upper plane. */
    assert(gsl_so2D_eval(result,1.5,1.5)==10.);
}
/* Verify native property thread equivalence and table-combination failure recovery. */
static void worker_audit(void) {
    double rows[8]={0.,1e10,2.,1.,.01,1e9,1.,.1},one[20],many[20];
    int status[2];
    assert(cg_properties(1,2,rows,one,status)==0 && !status[0] && !status[1]);
    assert(cg_properties(4,2,rows,many,status)==0 && !status[0] && !status[1]);
    assert(memcmp(one,many,sizeof(one))==0);
    assert(cg_properties(1,0,NULL,NULL,NULL)==CG_INVALID);
    assert(cg_spectra(1,NULL,NULL)==CG_INVALID);
    for(int repeat=0;repeat<10;repeat++) {
        int bad=0;
        assert(batch(2,2,status,combination_job,&bad)==0 && !status[0] && !status[1]);
        bad=1;
        assert(batch(2,2,status,combination_job,&bad)==0 && status[0]==CG_INVALID && status[1]==CG_INVALID);
    }
}
/* Run the native allocation, numerical-error and worker regression audit. */
int main(void) {
    worker_audit();
    gsl_error_handler_t *previous=gsl_set_error_handler_off();
    for(int repeat=0;repeat<100;repeat++) {
        assert(exercise(0,0)==0);
        /* 2 spline1D allocations, 3 spline2D, 2 LU, 1 raw buffer. */
        for(size_t fail=1;fail<=8;fail++) assert(exercise(fail,0)==2);
        assert(exercise(0,1)==3);
        assert(exercise(0,2)==3);
        assert(exercise(0,0)==0); /* Failure never poisons the next call. */
    }
    gsl_set_error_handler(previous);
    puts("Runtime audit passed: 1200 success/failure/recovery exercises.");
    return 0;
}
