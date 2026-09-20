/**
 * @file constants.h
 * @brief Physical constants and unit conversions for native calculations.
 *
 * Centralises particle masses, cross sections and conversion factors used by
 * the C physics kernels. Symbol suffixes identify units; consult the comments
 * below for the adopted values and conventions.
 */
#ifndef CONGRUENTS_PHYSICAL_CONSTANTS_H
#define CONGRUENTS_PHYSICAL_CONSTANTS_H

/*
 * Reconstructed constants header.
 *
 * Units follow the suffix convention used throughout the original source.
 * Values are CODATA 2018/IAU nominal values, rounded only beyond the precision
 * that affects the published spectra. Keep conversions centralized here.
 */

static const double c__cmsm1 = 2.99792458e10;
static const double h__GeVs = 4.135667696e-24;
static const double h__ergs = 6.62607015e-27;
static const double h__Js = 6.62607015e-34;
static const double k_B__GeVKm1 = 8.617333262e-14;
static const double k_B__ergKm1 = 1.380649e-16;

static const double m_e__GeV = 5.10998950e-4;
static const double m_e__g = 9.1093837015e-28;
static const double m_p__GeV = 9.3827208816e-1;
static const double m_H__g = 1.6735575e-24;
static const double m_H__kg = 1.6735575e-27;
static const double m_mu__GeV = 1.056583755e-1;
static const double m_pi0__GeV = 1.349768e-1;
static const double m_piC__GeV = 1.3957039e-1;
static const double m_d__GeV = 1.87561294257;

static const double sigma_T__mb = 6.6524587321e2;
static const double mb__cm2 = 1.0e-27;
static const double alpha = 7.2973525693e-3;
static const double e__esu = 4.80320471257e-10;
static const double gamma_Euler = 5.7721566490153286e-1;

/* Newton's G in the native galactic-dynamics units used by the model. */
static const double G__pcMsolm1km2sm2 = 4.30091727003628e-3;

static const double erg__GeV = 6.241509074e2;
static const double GeV__erg = 1.602176634e-3;
static const double pc__cm = 3.0856775814913673e18;
static const double Mpc__cm = 3.0856775814913673e24;
static const double yr__s = 3.15576e7;
static const double Msol__g = 1.988409870698051e33;
static const double Msol__kg = 1.988409870698051e30;
static const double Lsol__GeVsm1 = 3.828e33 * 6.241509074e2;

static const double T_0_CMB__K = 2.7255;
static const double arad__GeVcmm3Km4 = 7.5657e-15 * 6.241509074e2;

/* Historical spellings retained for source compatibility. */
#define sigma_T_mb sigma_T__mb
#define mb_cm2 mb__cm2
#define c_cmsm1 c__cmsm1
#define GeV_erg GeV__erg

#endif
