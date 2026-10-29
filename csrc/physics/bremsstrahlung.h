/**
 * @file bremsstrahlung.h
 * @brief Bremsstrahlung transition rates, cooling times and photon emission.
 *
 * Uses the supplied cross-section table to calculate electron transition
 * rates, integrate energy losses, and evaluate emission from the solved
 * electron population. Energies are in GeV, gas number densities in cm^-3,
 * and cooling times in seconds. Lookup tables are generated in Python;
 * interpolation and integration here run inside native galaxy workers.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Convert the bremsstrahlung cross-section table into an electron transition rate. */
static double P_BS__GeVm1sm1( double E_e__GeV, double E_f__GeV, double n_H__cmm3, gsl_spline_object_2D gso2D_BS  )
{
    if (E_f__GeV >= m_e__GeV && E_f__GeV < E_e__GeV)
    {
        double E_gam__GeV = E_e__GeV - E_f__GeV;
        return c__cmsm1 * n_H__cmm3 * gsl_so2D_eval( gso2D_BS, E_gam__GeV, E_e__GeV ) * mb__cm2;
    }
    else
    {
        return 0.;
    }
}

/* Return E/|dE/dt| from the complete tabulated loss moment, in seconds.
 * Cooling below the tracked grid is still a physical loss. Integrate down
 * to electron rest energy, not the solver's lower boundary. The legacy
 * bounds argument is retained for internal call compatibility, not used. */
static double tau_BS_fulltest__s(double energy, double bounds[2], double n_H__cmm3, gsl_spline_object_2D table) {
    double outgoing=0.,incoming=0.,loss=0.;
    table_transition_moments(table,energy,0.,energy-m_e__GeV,&outgoing,&incoming,&loss);
    loss*=c__cmsm1*n_H__cmm3*mb__cm2;
    return loss>0.?energy/loss:INFINITY;
}

/* Emission integrated over the solved electron population. */
struct cg_eps_BS_3_fdata_BS
    {
        double E_gam__GeV;
        double n_H__cmm3;
        gsl_spline_object_1D qess_so;
        gsl_spline_object_2D gso2D_BS;
    };

/* Evaluate batched integration samples for bremsstrahlung photon emission; state is passed explicitly. */
static int cg_eps_BS_3_F_BS( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_eps_BS_3_fdata_BS fdata_in = *((struct cg_eps_BS_3_fdata_BS *)fdata);
        for (j = 0; j < npts; ++j)
        {
            fval[j] = exp(x[j*ndim+0]) * gsl_so2D_eval( fdata_in.gso2D_BS, fdata_in.E_gam__GeV, exp(x[j*ndim+0]) ) * mb__cm2 *
                      c__cmsm1 * fdata_in.n_H__cmm3 * gsl_so1D_eval( fdata_in.qess_so, exp(x[j*ndim+0]) );
        }
        return 0;
    }

/* Integrate bremsstrahlung emission over electrons at the supplied gas density. */
static double eps_BS_3( double E_gam__GeV, double n_H__cmm3, gsl_spline_object_2D gso2D_BS, gsl_spline_object_1D qess_so )
{

    double res;
    double abserr;

    double xmin[1] = { log(fmax( E_CRe_lims__GeV[0], E_gam__GeV+m_e__GeV )) };
    double xmax[1] = { log(E_CRe_lims__GeV[1]) };

    struct cg_eps_BS_3_fdata_BS fdata;

    fdata.E_gam__GeV = E_gam__GeV;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.gso2D_BS = gso2D_BS;
    fdata.qess_so = qess_so;

    if (xmin[0] < xmax[0])
    {
        hcubature_v( 1, cg_eps_BS_3_F_BS, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );
    }
    else
    {
        res = 0.;
    }
    return res;
}
