/**
 * @file inverse_Compton.h
 * @brief Tabulated inverse Compton transitions, cooling and emission.
 *
 * Interpolates supplied IC transfer and emission tables, integrates electron
 * energy transfer to obtain cooling times, and computes photon emission
 * from the solved electron population. Energies use GeV and times seconds;
 * the lookup tables themselves are generated separately in Python.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Read the IC transition rate per final electron energy from its transfer table. */
static double P_IC__GeVm1sm1( double E_e__GeV, double E_f__GeV, gsl_spline_object_2D gso2D_IC )
{
    if (E_f__GeV < E_e__GeV)
    {
        double Delta_E__GeV = E_e__GeV - E_f__GeV;
        return gsl_so2D_eval( gso2D_IC, Delta_E__GeV, E_e__GeV );
    }
    else
    {
        return 0.;
    }
}

struct cg_tau_IC_fulltest__s_fdata_IC
    {
        double E_e__GeV;
        gsl_spline_object_2D gso_2D_so;
    };

/* Evaluate batched integration samples for IC energy-loss integration; state is passed explicitly. */
static int cg_tau_IC_fulltest__s_F_IC_out( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_tau_IC_fulltest__s_fdata_IC fdata_in = *((struct cg_tau_IC_fulltest__s_fdata_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] = exp(x[j*ndim+0]) * (fdata_in.E_e__GeV - exp(x[j*ndim+0])) *
                      P_IC__GeVm1sm1( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.gso_2D_so );
        }
        return 0;
    }

/* Integrate the tabulated IC energy transfer and return its cooling time in seconds. */
static double tau_IC_fulltest__s( double E_e__GeV, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so )
{

    struct cg_tau_IC_fulltest__s_fdata_IC fdata;
    double xmin[1], xmax[1];
    double res;
    double abserr;

    fdata.E_e__GeV = E_e__GeV;
    fdata.gso_2D_so = gso_2D_so;

    xmin[0] = log(E_e__GeV_lims[0]);
    xmax[0] = log(E_e__GeV);
    hcubature_v( 1, cg_tau_IC_fulltest__s_F_IC_out, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

    return E_e__GeV/res;

}

/* Emission integrated over the solved electron population. */
struct cg_eps_IC_3_fdata_IC
    {
        double E_gam__GeV;
        gsl_spline_object_1D qess_so;
        gsl_spline_object_2D gso2D_IC;
    };

/* Evaluate batched integration samples for IC photon emission; state is passed explicitly. */
static int cg_eps_IC_3_F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_eps_IC_3_fdata_IC fdata_in = *((struct cg_eps_IC_3_fdata_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] = exp(x[j*ndim+0]) * gsl_so1D_eval( fdata_in.qess_so, exp(x[j*ndim+0]) ) *
                      gsl_so2D_eval( fdata_in.gso2D_IC, fdata_in.E_gam__GeV, exp(x[j*ndim+0]) );
        }
        return 0;
    }

/* Integrate the IC emission kernel over the solved electron population. */
static double eps_IC_3( double E_gam__GeV, gsl_spline_object_2D gso2D_IC, gsl_spline_object_1D qess_so )
{

    double res = 0.;
    double abserr;
    struct cg_eps_IC_3_fdata_IC fdata;

    fdata.gso2D_IC = gso2D_IC;
    fdata.E_gam__GeV = E_gam__GeV;
    fdata.qess_so = qess_so;

    double xmin[1], xmax[1];
    xmin[0] = log(E_CRe_lims__GeV[0]);
    xmax[0] = log(E_CRe_lims__GeV[1]);

    if (xmin[0] < xmax[0])
    {
        hcubature_v( 1, cg_eps_IC_3_F_IC, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );
    }
    return res;
}
