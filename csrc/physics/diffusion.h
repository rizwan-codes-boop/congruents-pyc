/**
 * @file diffusion.h
 * @brief Electron diffusive escape rates and times.
 *
 * Interpolates the supplied diffusion coefficient D in cm^2/s and uses the
 * scale height h in pc to compute h^2/D in seconds. The negative E D/h^2
 * quantity is an escape sink expressed in GeV/s, not radiative cooling.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Represent diffusive escape as negative E D/h-squared in GeV/s. */
static double dEdtm1_diff__GeVsm1( double E_e__GeV, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1 )
{

    double D = gsl_spline_eval( gso_1D_D__cm2sm1.spline, E_e__GeV, gso_1D_D__cm2sm1.acc );
    return -1. * E_e__GeV * D/pow(h__pc * pc__cm,2);
}

/* Return the diffusion escape time h-squared/D in seconds. */
static double tau_diff__s( double E_e__GeV, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1 )
{
    return -1.* E_e__GeV/dEdtm1_diff__GeVsm1( E_e__GeV, h__pc, gso_1D_D__cm2sm1 );
}
