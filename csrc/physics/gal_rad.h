/**
 * @file gal_rad.h
 * @brief Galaxy radiation-field spectra, luminosities and dilution.
 *
 * Defines dust temperature and stellar/FIR luminosity prescriptions, photon
 * spectra and radiation dilution used in field combination and internal
 * absorption. Photon energies are in GeV and number spectra in cm^-3 GeV^-1.
 * Serial lookup-table generation is implemented separately in Python.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Estimate dust temperature in kelvin from redshift and specific star-formation rate. */
static double Tdust__K( double z, double SFR__Msolyrm1, double M_star__Msol )
{
    return 98.*pow(1.+z,-0.065) + 6.9*log10(SFR__Msolyrm1/M_star__Msol);
}

/* Estimate absorbed stellar luminosity in solar units from star-formation rate. */
static double Labs__Lsol( double SFR__Msolyrm1 )
  {
  return pow( 10., 1.096548 * log10(SFR__Msolyrm1/0.56) + 9.710084 );
  }

/* Assign the absorbed stellar luminosity to the far-infrared component. */
static double LFIR__Lsol( double SFR__Msolyrm1 )
  {
  return Labs__Lsol( SFR__Msolyrm1 );
  }

/* Estimate the observed old-stellar luminosity in solar units from stellar mass. */
static double Lobs_old__Lsol( double M_star__Msol )
  {
  return pow( 10., 0.8480565 * log10(M_star__Msol/0.56) + 1.521623 );
  }

/* Assign the 3000-K fraction of the old-stellar luminosity in solar units. */
static double L3000K__Lsol( double M_star__Msol )
  {
  return 0.574 * Lobs_old__Lsol( M_star__Msol );
  }

/* Assign the 4000-K fraction of the old-stellar luminosity in solar units. */
static double L4000K__Lsol( double M_star__Msol )
  {
  return 0.426 * Lobs_old__Lsol( M_star__Msol );
  }

/* Estimate the observed young-stellar luminosity from the piecewise SFR fit. */
static double Lobs_young__Lsol( double SFR__Msolyrm1 )
  {
  if (log10( SFR__Msolyrm1 ) > -2.6 ){return pow( 10., 0.7969616 * log10(SFR__Msolyrm1/0.56) + 9.007323 );}
  else {return pow( 10., 0.9867204 * log10(SFR__Msolyrm1/0.56) + 9.476592 );}
  }

/* Assign the 7500-K fraction of the young-stellar luminosity in solar units. */
static double L7500K__Lsol( double SFR__Msolyrm1 )
  {
  return 0.763 * Lobs_young__Lsol( SFR__Msolyrm1 );
  }

/* Assign the ultraviolet fraction of the young-stellar luminosity in solar units. */
static double LUV__Lsol( double SFR__Msolyrm1 )
  {
  return 0.237 * Lobs_young__Lsol( SFR__Msolyrm1 );
  }

/* Evaluate the piecewise UV photon number density in cm^-3 GeV^-1. */
static double dndEphot_UVMattis__cmm3GeVm1(  double *params, double E_phot__GeV )
{
    double lambda_micron = c__cmsm1 * h__GeVs/E_phot__GeV * 1.e4;
    if ( lambda_micron > 0.1340 && lambda_micron <= 0.2460 )
    {
        return 2.373/pow(E_phot__GeV,2) * erg__GeV * pow( lambda_micron, -0.6678 );
    }
    else if ( lambda_micron > 0.1100 && lambda_micron <= 0.1340 )
    {
        return 6.825e1/pow(E_phot__GeV,2) * erg__GeV * lambda_micron;
    }
    else if ( lambda_micron > 0.0912 && lambda_micron <= 0.1100 )
    {
        return 1.287e5/pow(E_phot__GeV,2) * erg__GeV * pow( lambda_micron, 4.4172 );
    }
    else
    {
        return 0.;
    }
}

/* Evaluate blackbody photon number density for kelvin temperature and GeV energy. */
static double dndE_BB__cmm3GeVm1( double T__K, double E_phot__GeV )
{
    return (8.*M_PI*pow(E_phot__GeV,2))/pow(h__GeVs * c__cmsm1, 3) * 1./( exp( E_phot__GeV/(k_B__GeVKm1 * T__K) ) - 1. );
}

/* Evaluate the modified-blackbody FIR photon number density in cm^-3 GeV^-1. */
static double dndE_modBB__cmm3GeVm1( double T__K, double E_phot__GeV )
{

    double E_0__GeV = 2e12 * h__GeVs;

    return dndE_BB__cmm3GeVm1( T__K, E_phot__GeV ) * E_phot__GeV/E_0__GeV;

}

/* Return integrated blackbody energy density in GeV/cm^3. */
static double u_rad_BB__GeVcmm3( double T__K )
{
    return arad__GeVcmm3Km4 * pow( T__K, 4 );
}

/* Return integrated modified-blackbody energy density in GeV/cm^3. */
static double u_rad_modBB__GeVcmm3( double T__K )
{

    double E_0__GeV = 2e12 * h__GeVs;
    return 24.8863 * 8. * M_PI/( pow( h__GeVs * c__cmsm1, 3 ) * E_0__GeV ) * pow( k_B__GeVKm1 * T__K, 5 );
}

/* Return the normalization energy density of the adopted UV spectrum in GeV/cm^3. */
static double u_rad_UVMattis__GeVcmm3( void )
{
    return 4450.1668;
}

/* Compute dimensionless radiation dilution using the adopted 2*pi*Re-squared geometry. */
static double C_dil( double u_rad__GeVcmm3, double L_obs__Lsol, double re__kpc, double h__pc )
  {

  return L_obs__Lsol * Lsol__GeVsm1/( u_rad__GeVcmm3 * 2. * M_PI * pow(re__kpc * 1e3 * pc__cm,2) * c__cmsm1 );
  }

/* Sum diluted stellar, UV and FIR photon fields with the redshifted CMB. */
static double dndEphot_total__cmm3GeVm1( double *params, double E_phot__GeV )
{
    double Cdil_3000, Cdil_4000, Cdil_7500, Cdil_UV, Cdil_FIR;
    double nphot_params[1];
    double dndEphot_total__cmm3GeVm1 = 0.;

    Cdil_FIR = C_dil( u_rad_modBB__GeVcmm3( params[1] ),  LFIR__Lsol( params[3] ), params[4], params[5] );

    Cdil_3000 = C_dil( u_rad_BB__GeVcmm3( 3000. ), L3000K__Lsol( params[2] ), params[4], params[5] );

    Cdil_4000 = C_dil( u_rad_BB__GeVcmm3( 4000. ),  L4000K__Lsol( params[2] ), params[4], params[5] );

    Cdil_7500 = C_dil( u_rad_BB__GeVcmm3( 7500. ), L7500K__Lsol( params[3] ), params[4], params[5] );

    Cdil_UV = C_dil( u_rad_UVMattis__GeVcmm3(), LUV__Lsol( params[3] ), params[4], params[5] );

    dndEphot_total__cmm3GeVm1 = dndE_BB__cmm3GeVm1( params[0], E_phot__GeV ) + Cdil_FIR * dndE_modBB__cmm3GeVm1( params[1], E_phot__GeV ) +
    Cdil_3000 * dndE_BB__cmm3GeVm1( 3000., E_phot__GeV ) + Cdil_4000 * dndE_BB__cmm3GeVm1( 4000., E_phot__GeV ) +
    Cdil_7500 * dndE_BB__cmm3GeVm1( 7500., E_phot__GeV ) + Cdil_UV * dndEphot_UVMattis__cmm3GeVm1( nphot_params, E_phot__GeV );

    return dndEphot_total__cmm3GeVm1;
}
