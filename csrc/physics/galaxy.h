/**
 * @file galaxy.h
 * @brief Empirical galaxy-property relations for the property worker.
 *
 * Defines gas and stellar velocity-dispersion estimates and the inverse
 * gas-surface-density relation. Function-name suffixes specify units:
 * dispersions use km/s, masses Msun, radii kpc and surface densities Msun/pc^2.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Estimate gas velocity dispersion in km/s from star-formation rate. */
static double sigma_gas_Yu__kmsm1( double SFR__Msolyrm1 )
{
    return pow(10., 0.2*log10( SFR__Msolyrm1 ) + 1.6 );
}

/* Evaluate the Sersic-dependent virial coefficient for stellar dispersion. */
static double cg_sigma_star_Bezanson__kmsm1_K_nu( double n ){return 73.32/(10.465+pow(n-0.94,2)) + 0.954;}

/* Estimate stellar velocity dispersion in km/s from mass and effective radius. */
static double sigma_star_Bezanson__kmsm1( double M_star__Msol, double Re__kpc )
{

    double n = 1.0;
    return sqrt( G__pcMsolm1km2sm2 * M_star__Msol/( 0.557 * cg_sigma_star_Bezanson__kmsm1_K_nu( n ) * Re__kpc * 1e3 ) );
}

/* Infer gas surface density from SFR and stellar surface densities. */
static double Sigma_gas_Shi_iKS__Msolpcm2( double Sigma_SFR__Msolyrm1pcm2, double Sigma_star__Msolpcm2 )
{
    return pow( 10., 10.28 ) * Sigma_SFR__Msolyrm1pcm2 * pow( Sigma_star__Msolpcm2, -0.48 );
}
