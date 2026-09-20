/**
 * @file ionisation.h
 * @brief Electron energy losses in neutral gas and ionised plasma.
 *
 * Defines neutral-medium ionisation and plasma loss rates in GeV/s, negative
 * for cooling, and corresponding cooling times in seconds. Electron energies
 * use GeV and gas number densities cm^-3. No photon emission is computed here.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Return neutral-medium electron ionisation dE/dt in GeV/s (negative for losses). */
static double dEdtm1_ion__GeVsm1( double E_e__GeV, double n_H__cmm3 )
{
    double gamma = E_e__GeV/m_e__GeV;
    double mu_ISM = 1.1;
    return -9./4. * c__cmsm1 * sigma_T__mb * mb__cm2 * m_e__GeV * n_H__cmm3 * mu_ISM * (log(gamma) + 0.91 * 2./3. *
           log( m_e__GeV * 1.e9/15.0 ) + 2. * 0.09 * 2./3. * log( m_e__GeV * 1.e9/41.5 ));
}

/* Convert the ionisation energy-loss rate to a positive cooling time in seconds. */
static double tau_ion__s( double E_e__GeV, double n_H__cmm3 )
{
    return -1.* E_e__GeV/dEdtm1_ion__GeVsm1( E_e__GeV, n_H__cmm3 );
}

/* Return plasma electron dE/dt in GeV/s using the local plasma frequency. */
static double dEdtm1_plasma__GeVsm1( double E_e__GeV, double n_H__cmm3 )
{
    double gamma = E_e__GeV/m_e__GeV;
    double mu_ISM = 1.1;
    double nu_p__Hz = sqrt( n_H__cmm3 * mu_ISM/( M_PI * m_e__g ) ) * e__esu;
    return -3./4. * c__cmsm1 * sigma_T__mb * mb__cm2 * m_e__GeV * n_H__cmm3 * mu_ISM * (log(gamma) + 2. * log(m_e__GeV/(h__GeVs * nu_p__Hz)));
}

/* Convert plasma energy losses to a cooling time in seconds. */
static double tau_plasma__s( double E_e__GeV, double n_H__cmm3 )
{
    return -1.* E_e__GeV/dEdtm1_plasma__GeVsm1( E_e__GeV, n_H__cmm3 );
}
