/**
 * @file free_free.h
 * @brief Thermal free-free optical depth and photon emission.
 *
 * Defines the worker-side free-free absorption and emission prescriptions.
 * Inputs include photon energy in GeV, electron temperature in kelvin,
 * SFR surface density in Msun/yr/pc^2 and emitting radius in kpc.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Calculate dimensionless free-free optical depth from photon energy, SFR density and gas temperature. */
static double tau_FF_MK( double E_gam__GeV, double Sigma_SFR__Msolyrm1pcm2, double T_e__K )
{
    double Z = 1.;
    double Phi = 4.2e60;
    double nu_gam__Hz = E_gam__GeV/h__GeVs;
    double g_ff;

    if ( nu_gam__Hz > 1e9 )
    {

        g_ff = log( exp(5.960 - sqrt(3)/M_PI * log( Z*( nu_gam__Hz/1e9 )/pow(T_e__K/1e4, 3./2.) ))+ exp(1.) );
    }
    else
    {

        g_ff = sqrt(3.)/M_PI*( log( pow( 2.* k_B__ergKm1 * T_e__K , 1.5 )/( M_PI * Z * pow(e__esu,2) * pow(m_e__g,0.5) * nu_gam__Hz ) ) - 5.*gamma_Euler/2. );
    }

    double alpha_nuff = (4. * pow(e__esu,6))/(3. * c__cmsm1 * h__ergs) * pow( (2.*M_PI)/(3.*pow(m_e__g,3)*k_B__ergKm1*T_e__K), 0.5 ) *
                        pow(Z, 2) * (1.-exp( -(h__ergs* nu_gam__Hz)/(k_B__ergKm1*T_e__K) ) )/pow(nu_gam__Hz,3) * g_ff;

    double alpha_B = 2.54e-13 * pow(Z,2) * pow( T_e__K/1e4 / pow(Z,2), -0.8163 - 0.0208*log( T_e__K/1e4 / pow(Z,2) ) );

    double param_a = (3. * Sigma_SFR__Msolyrm1pcm2/(yr__s * pow(pc__cm,2)) * Phi * alpha_nuff)/(8. * alpha_B);

    return param_a;
}

/* Calculate thermal free-free emission using the emitting radius and optical depth. */
static double eps_FF( double E_gam__GeV, double Re__kpc, double T_e__K, double tau_ff )
{
    if ( tau_ff > 1.e-6 )
    {
        return dndE_BB__cmm3GeVm1( T_e__K, E_gam__GeV ) * c__cmsm1 * 4.*M_PI*pow(Re__kpc*1.e3*pc__cm,2) * (1.-exp(-tau_ff));
    }
    else
    {
        return dndE_BB__cmm3GeVm1( T_e__K, E_gam__GeV ) * c__cmsm1 * 4.*M_PI*pow(Re__kpc*1.e3*pc__cm,2) * (tau_ff - pow(tau_ff,2)/2. + pow(tau_ff,3)/6.);
    }
}
