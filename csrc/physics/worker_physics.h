/**
 * @file worker_physics.h
 * @brief Internal declarations and includes for galaxy-worker physics.
 *
 * Collects shared integration state, reference energy bounds and forward
 * declarations, then includes the per-process implementations. This is a
 * private implementation header for solver.c, not the public Python-C API.
 * The included physics functions are defined directly in their headers.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
#include "constants.h"
#include <gsl/gsl_sf_hyperg.h>

/* The threshold proton energy for pion production, in GeV. */
#define T_p_th (2. * m_pi0__GeV + pow(m_pi0__GeV, 2)/(2. * m_p__GeV))

extern double q_p_inject;
extern double T_CR_lims__GeV[2];
extern double E_CRe_lims__GeV[2];

/* The lower and upper bounds for the proton energy in the simulation, in GeV. */
double T_p_low = 1.e-3;
double T_p_high = 1.e8;
double T_p_norm__GeV[2] = { 0., 1.e8 };

/* The lower and upper bounds for the electron energy in the simulation, in GeV. */
double E_e_low_norm = 1.e-3;
double E_e_high_norm = 1.e8;
double K_pi = 0.17;
struct F_int_data
{
    double E_e__GeV;
    double E_f__GeV;
    double n_H__cmm3;
    double B__G;
    double h__pc;
    int n_gso2D;
    gsl_spline_object_2D * gso_2D_radfield;
    gsl_spline_object_2D gso2D_BS;
    gsl_spline_object_1D gso_1D_Q;
    gsl_spline_object_1D gso_1D_D__cm2sm1;

    double (*n_phot)( double *, double );
    double (*E_func)( double, double, double, double );

    double *n_phot_params;
};

/* Normalize the pion-decay differential gamma-ray cross section, including enhancement. */
static double A_max( double T_p );

/* Assemble and solve the coupled primary/secondary steady-state electron equations. */
static int CRe_steadystate_solve( int structure, double E_e_lim__GeV[2], int n_E, double n_H__cmm3, double B__G, double h__pc,
    unsigned int n_gso2D, gsl_spline_object_2D * gso_2D_radfields, gsl_spline_object_2D gso2D_BS, gsl_spline_object_1D gso_1D_D__cm2sm1,
    gsl_spline_object_1D gso_1D_Q_inject_1, gsl_spline_object_1D gso_1D_Q_inject_2,
    gsl_spline_object_1D * qe_1_so_1D, gsl_spline_object_1D * qe_2_so_1D );

/* Compute dimensionless radiation dilution using the adopted 2*pi*Re-squared geometry. */
static double C_dil( double u_rad__GeVcmm3, double L_obs__Lsol, double re__kpc, double h__pc );

/* Integrate the energy-weighted injection shape to obtain its normalization. */
static double C_norm_E( double q, double m, double T_cutoff );

/* Return the kinematic neutral-pion energy in the centre-of-mass frame, GeV. */
static double E_pi_CM( double T_p );

/* Transform the maximum neutral-pion energy to the laboratory frame, GeV. */
static double E_pi_LAB_max( double T_p );

/* Evaluate escape and continuous-loss moments for the electron matrix. */
static int F_EdotDE_2_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval );

/* Evaluate four log-energy moments of nonlocal electron transitions for matrix assembly. */
static int F_Gamma_2D_4_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval );

/* Evaluate transition-rate moments for the lower-boundary electron matrix contribution. */
static int F_Gamma_i0_2D_2_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval );

/* Evaluate the energy-squared injection term in logarithmic electron-energy coordinates. */
static int F_QE2_i_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval );

/* Return the unit step, taking the value one at zero. */
static double Heaviside( double arg );

/* Evaluate the momentum-power-law particle spectrum with an exponential kinetic-energy cutoff. */
static double J( double T, double C, double q, double m, double T_cutoff );

/* Assign the 3000-K fraction of the old-stellar luminosity in solar units. */
static double L3000K__Lsol( double M_star__Msol );

/* Assign the 4000-K fraction of the old-stellar luminosity in solar units. */
static double L4000K__Lsol( double M_star__Msol );

/* Assign the 7500-K fraction of the young-stellar luminosity in solar units. */
static double L7500K__Lsol( double SFR__Msolyrm1 );

/* Assign the absorbed stellar luminosity to the far-infrared component. */
static double LFIR__Lsol( double SFR__Msolyrm1 );

/* Assign the ultraviolet fraction of the young-stellar luminosity in solar units. */
static double LUV__Lsol( double SFR__Msolyrm1 );

/* Estimate absorbed stellar luminosity in solar units from star-formation rate. */
static double Labs__Lsol( double SFR__Msolyrm1 );

/* Estimate the observed old-stellar luminosity in solar units from stellar mass. */
static double Lobs_old__Lsol( double M_star__Msol );

/* Estimate the observed young-stellar luminosity from the piecewise SFR fit. */
static double Lobs_young__Lsol( double SFR__Msolyrm1 );

/* Convert the bremsstrahlung cross-section table into an electron transition rate. */
static double P_BS__GeVm1sm1( double E_e__GeV, double E_f__GeV, double n_H__cmm3, gsl_spline_object_2D gso2D_BS  );

/* Read the IC transition rate per final electron energy from its transfer table. */
static double P_IC__GeVm1sm1( double E_e__GeV, double E_f__GeV, gsl_spline_object_2D gso2D_IC );

/* Return the neutral-pion centre-of-mass momentum in natural GeV units. */
static double P_pi_CM( double T_p );

/* Infer gas surface density from SFR and stellar surface densities. */
static double Sigma_gas_Shi_iKS__Msolpcm2( double Sigma_SFR__Msolyrm1pcm2, double Sigma_star__Msolpcm2 );

/* Estimate dust temperature in kelvin from redshift and specific star-formation rate. */
static double Tdust__K( double z, double SFR__Msolyrm1, double M_star__Msol );

/* Return the centre-of-mass speed as a fraction of the speed of light. */
static double beta_CM( double T_p );

/* Represent diffusive escape as negative E D/h-squared in GeV/s. */
static double dEdtm1_diff__GeVsm1( double E_e__GeV, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1 );

/* Return neutral-medium electron ionisation dE/dt in GeV/s (negative for losses). */
static double dEdtm1_ion__GeVsm1( double E_e__GeV, double n_H__cmm3 );

/* Return plasma electron dE/dt in GeV/s using the local plasma frequency. */
static double dEdtm1_plasma__GeVsm1( double E_e__GeV, double n_H__cmm3 );

/* Return synchrotron electron dE/dt in GeV/s for a magnetic field in gauss. */
static double dEdtm1_sync__GeVsm1( double E_e__GeV, double B__G );

/* Sum continuous synchrotron and neutral ionisation losses for the disc solver. */
static double dEdtm1_total_disc__GeVsm1( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc );

/* Sum continuous synchrotron and plasma losses for the halo solver. */
static double dEdtm1_total_halo__GeVsm1( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc );

/* Sum IC and bremsstrahlung transitions per logarithmic final-electron-energy interval. */
static double dGammadlogEf_total__logGeVm1sm1( double E_e__GeV, double E_f__GeV, double n_gso2D, gsl_spline_object_2D * gso_2D_radfield, double n_H__cmm3, gsl_spline_object_2D gso2D_BS );

/* Evaluate blackbody photon number density for kelvin temperature and GeV energy. */
static double dndE_BB__cmm3GeVm1( double T__K, double E_phot__GeV );

/* Evaluate the modified-blackbody FIR photon number density in cm^-3 GeV^-1. */
static double dndE_modBB__cmm3GeVm1( double T__K, double E_phot__GeV );

/* Evaluate the piecewise UV photon number density in cm^-3 GeV^-1. */
static double dndEphot_UVMattis__cmm3GeVm1(  double *params, double E_phot__GeV );

/* Sum diluted stellar, UV and FIR photon fields with the redshifted CMB. */
static double dndEphot_total__cmm3GeVm1( double *params, double E_phot__GeV );

/* Return the differential pion-decay photon cross section for GeV energies. */
static double dsig_dEg( double T_p, double E_gam );

/* Integrate bremsstrahlung emission over electrons at the supplied gas density. */
static double eps_BS_3( double E_gam__GeV, double n_H__cmm3, gsl_spline_object_2D gso2D_BS, gsl_spline_object_1D qess_so );

/* Calculate thermal free-free emission using the emitting radius and optical depth. */
static double eps_FF( double E_gam__GeV, double Re__kpc, double T_e__K, double tau_ff );

/* Integrate the IC emission kernel over the solved electron population. */
static double eps_IC_3( double E_gam__GeV, gsl_spline_object_2D gso2D_IC, gsl_spline_object_1D qess_so );

/* Integrate the synchrotron kernel over a solved electron population at GeV photon energy. */
static double eps_SY_4( double E_gam__GeV, double B__G, gsl_spline_object_1D sync_x_so, gsl_spline_object_1D qess_so );

/* Integrate pion-decay photon emission with energy-dependent proton calorimetry. */
static double eps_pi( double E_gam__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal );

/* Integrate the diagnostic pion-decay spectrum assuming complete proton calorimetry. */
static double eps_pi_fcal1( double E_gam__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal );

/* Assemble the electron-neutrino decay spectrum from its piecewise branches. */
static double f_nu_e( double x );

/* Assemble the secondary muon-neutrino decay spectrum from its piecewise branches. */
static double f_nu_mu2( double x );

/* Return the proton collision centre-of-mass Lorentz factor. */
static double gam_CM( double T_p );

/* Integrate secondary electron injection at the requested kinetic energy in GeV. */
static double q_e( double T_e__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal );

/* Integrate the pion and muon decay contributions to the neutrino production spectrum. */
static double q_nu( double E_nu__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal );

/* Calculate pion production at GeV energy, returning zero outside calorimetry coverage. */
static double q_pi( double E_pi__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal );

/* Return the squared proton collision centre-of-mass energy in GeV squared. */
static double s( double T_p );

/* Estimate gas velocity dispersion in km/s from star-formation rate. */
static double sigma_gas_Yu__kmsm1( double SFR__Msolyrm1 );

/* Evaluate the Breit-Wheeler photon-pair cross section in millibarns above threshold. */
static double sigma_gg_BW__mb( double E1__GeV, double E2__GeV );

/* Evaluate the fitted inelastic proton-proton cross section in millibarns. */
static double sigma_inel__mb( double T_p );

/* Estimate stellar velocity dispersion in km/s from mass and effective radius. */
static double sigma_star_Bezanson__kmsm1( double M_star__Msol, double Re__kpc );

/* Integrate bremsstrahlung energy transfer and return its cooling time in seconds. */
static double tau_BS_fulltest__s( double E_e__GeV, double E_e__GeV_lims[2], double n_H__cmm3, gsl_spline_object_2D gso2D_BS );

/* Calculate dimensionless free-free optical depth from photon energy, SFR density and gas temperature. */
static double tau_FF_MK( double E_gam__GeV, double Sigma_SFR__Msolyrm1pcm2, double T_e__K );

/* Integrate the tabulated IC energy transfer and return its cooling time in seconds. */
static double tau_IC_fulltest__s( double E_e__GeV, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so );

/* Return the diffusion escape time h-squared/D in seconds. */
static double tau_diff__s( double E_e__GeV, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1 );

/* Integrate internal pair-production optical depth over the target radiation spectrum. */
static double tau_gg_gal_BW( double E_gam__GeV, double (*n_phot)(double *, double), double *n_phot_params, double E_phot__GeV_lims[2], double h_pc );

/* Convert the ionisation energy-loss rate to a positive cooling time in seconds. */
static double tau_ion__s( double E_e__GeV, double n_H__cmm3 );

/* Convert plasma energy losses to a cooling time in seconds. */
static double tau_plasma__s( double E_e__GeV, double n_H__cmm3 );

/* Convert synchrotron losses to an electron cooling time in seconds. */
static double tau_sync__s( double E_e__GeV, double B__G );

/* Return integrated blackbody energy density in GeV/cm^3. */
static double u_rad_BB__GeVcmm3( double T__K );

/* Return the normalization energy density of the adopted UV spectrum in GeV/cm^3. */
static double u_rad_UVMattis__GeVcmm3( void );

/* Return integrated modified-blackbody energy density in GeV/cm^3. */
static double u_rad_modBB__GeVcmm3( double T__K );

#include "gal_rad.h"
#include "ionisation.h"
#include "diffusion.h"
#include "synchrotron.h"
#include "inverse_Compton.h"
#include "bremsstrahlung.h"
#include "CR_funcs.h"
#include "hadronic.h"
#include "free_free.h"
#include "absorption.h"
#include "galaxy.h"
#include "steady_state_native.h"
