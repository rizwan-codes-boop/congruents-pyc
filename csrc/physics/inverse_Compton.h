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
    if (E_f__GeV >= m_e__GeV && E_f__GeV < E_e__GeV)
    {
        double Delta_E__GeV = E_e__GeV - E_f__GeV;
        return gsl_so2D_eval( gso2D_IC, Delta_E__GeV, E_e__GeV );
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
static double tau_IC_fulltest__s(double energy, double bounds[2], gsl_spline_object_2D table) {
    double outgoing=0.,incoming=0.,loss=0.;
    table_transition_moments(table,energy,0.,energy-m_e__GeV,&outgoing,&incoming,&loss);
    return loss>0.?energy/loss:INFINITY;
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
