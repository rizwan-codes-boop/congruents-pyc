/**
 * @file solver.h
 * @brief Public C interface and buffer layouts for the Python wrapper.
 *
 * Declares synchronous galaxy-property, transport and spectrum entry points,
 * interface version checks and error codes. Contiguous buffers are borrowed
 * from the caller and must remain valid until the call returns. Shapes,
 * units and component ordering are documented with the structures below.
 */
#ifndef CG_SOLVER_H
#define CG_SOLVER_H
#include <stddef.h>
#define CG_API __attribute__((visibility("default")))
enum { CG_OK=0, CG_INVALID=1, CG_ALLOC=2, CG_NUMERIC=3 };
/* Borrowed contiguous float64 buffers; owners must outlive synchronous calls.
 * Tables use values[ny][nx]. temperature is kelvin for CMB/FIR, zero otherwise.
 */
typedef struct { size_t nx,ny; const double *x,*y,*values; double temperature; } cg_table_input;
/* Catalogue rows[n][4]: z, mass Msun, radius kpc, SFR Msun/year.
 * props[n][10]: h pc,nH cm^-3,B G,sigma km/s,area pc²,SigmaGas,
 * SigmaSFR,SigmaStar,Tdust K,Bhalo G. Surface units follow README.
 * transport[n][7][ne]: fcal,Dp,Ddisc,Dhalo,Qprimary,Qsecondary,Nproton.
 * electrons[n][4][ne]: primary/secondary disc, primary/secondary halo.
 * emission[n][15][np]: IC1d,IC2d,BS1d,BS2d,SY1d,SY2d,IC1h,IC2h,
 * SY1h,SY2h,FF,tauFF,pi,pi_fcal1,nu. Disc components include FF absorption.
 * loss[n][12][ne]: disc SY/BS/IC/DI/IO, halo same, proton PP/DI (seconds).
 * critical[n][10] GeV/s; budget[n][16] powers then fractions;
 * radio[n][5] W/Hz; escape[n][2][ne] GeV^-1 s^-1; internal[n][np] dimensionless.
 * Diagnostic pointers may be NULL only when diagnostics=0.
 * Failed rows must be discarded; statuses[n] is authoritative.
 */
typedef struct {
    size_t n,ne,np,ns,ncmb,nfir;
    const double *rows,*props,*kinetic,*electron,*photon,*cp;
    double *transport,*electrons,*emission,*internal;
    const cg_table_input *ic,*gamma,*bs,*sy;
    double target_low,target_high;
    int diagnostics;
    double *loss,*critical,*budget,*radio,*escape;
} cg_run;
/* Return the interface version used to validate Python declarations. */
CG_API unsigned cg_solver_abi(void);
/* Report whether this library was compiled with active OpenMP workers. */
CG_API int cg_solver_openmp_enabled(void);
/* Validate catalogue buffers and dispatch the parallel galaxy-property calculation. */
CG_API int cg_properties(int,size_t,const double *,double *,int *);
/* Validate transport inputs and dispatch galaxy jobs into Python-owned output buffers. */
CG_API int cg_transport(int,size_t,size_t,const double *,const double *,const double *,double *,double *,int *);
/* Validate shared tables and run complete independent galaxy-spectrum workers. */
CG_API int cg_spectra(int,const cg_run *,int *);
#endif
