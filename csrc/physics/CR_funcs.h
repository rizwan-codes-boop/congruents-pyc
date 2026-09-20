/**
 * @file CR_funcs.h
 * @brief Combined electron transition rates and continuous energy losses.
 *
 * Combines tabulated inverse Compton and bremsstrahlung transitions for
 * the kinetic solver. Continuous losses combine synchrotron with neutral
 * ionisation in the disc and plasma losses in the halo; dE/dt is in GeV/s.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Sum IC and bremsstrahlung transitions per logarithmic final-electron-energy interval. */
static double dGammadlogEf_total__logGeVm1sm1( double E_e__GeV, double E_f__GeV, double n_gso2D, gsl_spline_object_2D * gso_2D_radfield, double n_H__cmm3, gsl_spline_object_2D gso2D_BS )
{
    double P_total = 0.;
    unsigned short int i;

    if ( E_e__GeV > E_f__GeV )
    {
        for (i = 0; i < n_gso2D; ++i)
        {
            P_total += P_IC__GeVm1sm1( E_e__GeV, E_f__GeV, gso_2D_radfield[i] ) * E_f__GeV;
        }

        P_total += P_BS__GeVm1sm1( E_e__GeV, E_f__GeV, n_H__cmm3, gso2D_BS ) * E_f__GeV;
        return P_total;
    }
    else
    {
        return 0.;
    }
}

/* Sum continuous synchrotron and neutral ionisation losses for the disc solver. */
static double dEdtm1_total_disc__GeVsm1( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc )
{
    return dEdtm1_sync__GeVsm1( E_e__GeV, B__G ) + dEdtm1_ion__GeVsm1( E_e__GeV, n_H__cmm3 );
}

/* Sum continuous synchrotron and plasma losses for the halo solver. */
static double dEdtm1_total_halo__GeVsm1( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc )
{
    return dEdtm1_sync__GeVsm1( E_e__GeV, B__G ) + dEdtm1_plasma__GeVsm1( E_e__GeV, n_H__cmm3 );
}
