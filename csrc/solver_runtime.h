/**
 * @file solver_runtime.h
 * @brief Private numerical support, memory ownership and worker recovery.
 *
 * Defines thread-local worker state, tracked allocation/cleanup, checked
 * interpolation and integration, and GSL LU solving. Numerical failures are
 * recorded for the current galaxy; failed outputs must not be used. Included
 * only by solver.c; OpenMP dispatch is implemented there, not in this header.
 */
#ifndef CG_SOLVER_RUNTIME_H
#define CG_SOLVER_RUNTIME_H
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <setjmp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_linalg.h>
#include <cubature.h>

typedef void (*cg_destructor)(void *);
typedef struct {
    jmp_buf escape;
    int status;
    const char *phase;
    size_t count;
#ifdef CG_TESTING
    size_t fail_allocation;
#endif
    void *owned[8192];
    cg_destructor destroy[8192];
} cg_worker;
static _Thread_local cg_worker *cg_active;

/* Record the current galaxy's failure and unwind to its cleanup boundary. */
static void cg_fail(int status) {
    if(getenv("CONGRUENTS_DEBUG")) fprintf(stderr,"solver failure %d in %s\n",status,cg_active->phase);
    cg_active->status = status;
    longjmp(cg_active->escape, 1);
}
/* Register a native allocation and destructor with the current worker. */
static void *cg_own(void *p, cg_destructor destroy) {
    if (!p) cg_fail(2);
#ifdef CG_TESTING
    if (cg_active->fail_allocation && --cg_active->fail_allocation == 0) {
        destroy(p);
        cg_fail(2);
    }
#endif
    if (cg_active->count == 8192) {
        destroy(p);
        cg_fail(2);
    }
    size_t i = cg_active->count++;
    cg_active->owned[i] = p;
    cg_active->destroy[i] = destroy;
    return p;
}
/* Release one tracked allocation without leaving a dangling cleanup entry. */
static void cg_release(void *p) {
    if (!p) return;
    for (size_t i=cg_active->count; i>0; --i)
        if (cg_active->owned[i-1] == p) {
            cg_active->destroy[i-1](p);
            cg_active->owned[i-1] = NULL;
            return;
        }
}
/* Destroy remaining worker allocations in reverse order after success or failure. */
static void cg_cleanup(cg_worker *w) {
    for (size_t i=w->count; i>0; --i)
        if (w->owned[i-1]) {
            w->destroy[i-1](w->owned[i-1]);
            w->owned[i-1]=NULL;
        }
    w->count=0;
}
/* Allocate bytes and register ownership for automatic worker cleanup. */
static void *cg_allocate(size_t n) { return cg_own(malloc(n), free); }
/* Release the tracked rows and outer pointer of a two-dimensional array. */
static void cg_free2D(int n, double **p) {
    for (int i=0;i<n;i++) cg_release(p[i]);
    cg_release(p);
}

/* Override only the adapter's spline wrappers, not the scientific reference. */
#define CONGRUENTS_GSL_DECS_H
typedef struct { gsl_spline *spline; gsl_interp_accel *acc; double x_lim[2]; } gsl_spline_object_1D;
typedef struct {
    gsl_spline2d *spline;
    gsl_interp_accel *xacc, *yacc;
    double x_lim[2], y_lim[2];
} gsl_spline_object_2D;
/* Destroy a one-dimensional GSL spline through the worker ownership registry. */
static void cg_free_spline(void *p) { gsl_spline_free(p); }
/* Destroy a two-dimensional GSL spline through the worker ownership registry. */
static void cg_free_spline2d(void *p) { gsl_spline2d_free(p); }
/* Destroy a GSL interpolation accelerator through the worker ownership registry. */
static void cg_free_acc(void *p) { gsl_interp_accel_free(p); }
/* Destroy a GSL vector through the worker ownership registry. */
static void cg_free_vector(void *p) { gsl_vector_free(p); }
/* Destroy a GSL permutation through the worker ownership registry. */
static void cg_free_permutation(void *p) { gsl_permutation_free(p); }
/* Create a worker-owned linear spline and its private interpolation accelerator. */
static gsl_spline_object_1D gsl_so1D(size_t n,const double *x,const double *z) {
    gsl_spline_object_1D o;
    o.acc=cg_own(gsl_interp_accel_alloc(),cg_free_acc);
    o.spline=cg_own(gsl_spline_alloc(gsl_interp_linear,n),cg_free_spline);
    if(gsl_spline_init(o.spline,x,z,n)) cg_fail(3);
    o.x_lim[0]=x[0]; o.x_lim[1]=x[n-1];
    return o;
}
/* Evaluate a one-dimensional spline using the runtime's domain and error policy. */
static double gsl_so1D_eval(gsl_spline_object_1D o,double x) {
    double y=NAN;
    /* Do not jump out of a cubature callback: it must free its own workspace. */
    if((x<o.x_lim[0] || x>o.x_lim[1]) && getenv("CONGRUENTS_DEBUG"))
        fprintf(stderr,"spline range %.17g [%.17g, %.17g]\n",x,o.x_lim[0],o.x_lim[1]);
    if(x<o.x_lim[0] || x>o.x_lim[1] ||
       gsl_spline_eval_e(o.spline,x,o.acc,&y)) return NAN;
    return y;
}
/* Release a one-dimensional spline and its tracked accelerator. */
static void gsl_so1D_free(gsl_spline_object_1D o) {
    cg_release(o.spline); cg_release(o.acc);
}
/* Create a worker-owned bilinear spline with independent axis accelerators. */
static gsl_spline_object_2D gsl_so2D(size_t nx,size_t ny,const double *x,const double *y,const double *z) {
    gsl_spline_object_2D o;
    o.xacc=cg_own(gsl_interp_accel_alloc(),cg_free_acc);
    o.yacc=cg_own(gsl_interp_accel_alloc(),cg_free_acc);
    o.spline=cg_own(gsl_spline2d_alloc(gsl_interp2d_bilinear,nx,ny),cg_free_spline2d);
    if(gsl_spline2d_init(o.spline,x,y,z,nx,ny)) cg_fail(3);
    o.x_lim[0]=x[0];o.x_lim[1]=x[nx-1];o.y_lim[0]=y[0];o.y_lim[1]=y[ny-1];
    return o;
}
/* Create a two-dimensional spline and attach its temperature coordinate. */
static gsl_spline_object_2D gsl_so2D_temp(size_t nx,size_t ny,const double *x,const double *y,
                                        const double xl[2],const double yl[2],const double *z) {
    gsl_spline_object_2D o=gsl_so2D(nx,ny,x,y,z);
    o.x_lim[0]=xl[0];o.x_lim[1]=xl[1];o.y_lim[0]=yl[0];o.y_lim[1]=yl[1];
    return o;
}
/* Evaluate a bilinear spline; return zero outside its tabulated domain. */
static double gsl_so2D_eval(gsl_spline_object_2D o,double x,double y) {
    double v=NAN;
    if(x<o.x_lim[0] || x>o.x_lim[1] || y<o.y_lim[0] || y>o.y_lim[1]) return 0.;
    if(gsl_spline2d_eval_e(o.spline,x,y,o.xacc,o.yacc,&v)) return NAN;
    return v;
}
/* Release a two-dimensional spline and both tracked accelerators. */
static void gsl_so2D_free(gsl_spline_object_2D o) {
    cg_release(o.spline); cg_release(o.xacc); cg_release(o.yacc);
}
typedef struct { integrand_v function; void *data; } cg_integrand;
/* Validate batched integrand values while allowing cubature to clean up on failure. */
static int cg_integrand_call(unsigned ndim,size_t n,const double *x,void *p,unsigned dim,double *v) {
    cg_integrand *f=p;
    int status=cg_active->status;
    if(!status) status=f->function(ndim,n,x,f->data,dim,v);
    if(!status) for(size_t i=0;i<n*dim;i++) if(!isfinite(v[i])) {
        status=1;
        break;
    }
    if(status) {
        /* Legacy cubature leaks regions on callback failure. Mark the entire
         * galaxy invalid, but let cubature leave through its normal cleanup.
         * These zeros are NEVER exposed as valid scientific results. */
        cg_active->status=3;
        for(size_t i=0;i<n*dim;i++) v[i]=0.;
    }
    return 0;
}
/* Run checked adaptive integration and propagate non-finite or failed results. */
static int cg_cubature(unsigned dim,integrand_v f,void *data,unsigned ndim,const double *a,const double *b,
                       size_t maxeval,double abs,double rel,error_norm norm,double *v,double *err) {
    cg_integrand checked={f,data};
    int status=cg_active->status;
    if(!status) status=hcubature_v(dim,cg_integrand_call,&checked,ndim,a,b,maxeval,abs,rel,norm,v,err);
    if(cg_active->status) status=1;
    if(!status) for(unsigned i=0;i<dim;i++) if(!isfinite(v[i])) status=1;
    if(status) {
        /* Return normally through nested callbacks' enclosing functions so
         * GNU trampoline cleanup is not skipped by a non-local jump. */
        cg_active->status=3;
        for(unsigned i=0;i<dim;i++) { v[i]=NAN;err[i]=NAN; }
        return 1;
    }
    return 0;
}
/* Solve a finite linear system by GSL LU decomposition and clip negative populations. */
static int cg_linear_solve(int n,double *a,double *b,double *out) {
    for(int i=0;i<n*n;i++) if(!isfinite(a[i])) cg_fail(3);
    for(int i=0;i<n;i++) if(!isfinite(b[i])) cg_fail(3);
    gsl_matrix_view m=gsl_matrix_view_array(a,n,n);
    gsl_vector_view rhs=gsl_vector_view_array(b,n);
    gsl_vector *x=cg_own(gsl_vector_alloc(n),cg_free_vector);
    gsl_permutation *p=cg_own(gsl_permutation_alloc(n),cg_free_permutation);
    int sign;
    int lu=gsl_linalg_LU_decomp(&m.matrix,p,&sign);
    if(!lu) lu=gsl_linalg_LU_solve(&m.matrix,p,&rhs.vector,x);
    if(lu) {
        if(getenv("CONGRUENTS_DEBUG")) fprintf(stderr,"LU status %d\n",lu);
        cg_fail(3);
    }
    for(int i=0;i<n;i++) {
        double v=gsl_vector_get(x,i);
        if(!isfinite(v)) {
            if(getenv("CONGRUENTS_DEBUG")) fprintf(stderr,"nonfinite LU result\n");
            cg_fail(3);
        }
        out[i]=fmax(0.,v);
    }
    cg_release(p);cg_release(x);
    return 0;
}
#endif
