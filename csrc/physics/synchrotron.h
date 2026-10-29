/**
 * @file synchrotron.h
 * @brief Synchrotron electron cooling and photon emission.
 *
 * Defines energy-loss rates in GeV/s and cooling times in seconds, and
 * integrates the supplied synchrotron kernel over the solved electron
 * population. Electron/photon energies use GeV and magnetic fields gauss.
 * The shared kernel table is generated separately in Python.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */

/* Return synchrotron electron dE/dt in GeV/s for a magnetic field in gauss. */
static double dEdtm1_sync__GeVsm1( double E_e__GeV, double B__G )
{
    double gamma = E_e__GeV/m_e__GeV;
    double beta = sqrt(1. - 1./pow(gamma,2));
//the older expression was missing the factor of 4/3.
    return -4./3. * sigma_T__mb * mb__cm2 * c__cmsm1 * pow(gamma,2) * pow(beta,2) * pow(B__G, 2)/(8.*M_PI) * erg__GeV;
}

/* Convert synchrotron losses to an electron cooling time in seconds. */
static double tau_sync__s( double E_e__GeV, double B__G )
{
    return -1.* E_e__GeV/dEdtm1_sync__GeVsm1( E_e__GeV, B__G );
}

/* Emission integrated over the solved electron population. */
struct cg_eps_SY_4_fdata_sync
    {
        double xE2;
        gsl_spline_object_1D qess_so;
        gsl_spline_object_1D sync_x_so;
    };

/* Evaluate batched integration samples for synchrotron emission; state is passed explicitly. */
static int cg_eps_SY_4_F_SY( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_eps_SY_4_fdata_sync fdata_in = *((struct cg_eps_SY_4_fdata_sync *)fdata);
        double t;

        for (j = 0; j < npts; ++j)
        {
            t = fdata_in.xE2/pow(exp(x[j*ndim+0]),2);
            /* R(x) is isotropically pitch-angle averaged in Python. Retain
             * beta^2 so its bolometric power matches the cooling expression;
             * the spectral shape itself remains a relativistic approximation. */
            double beta2 = 1. - pow(m_e__GeV/exp(x[j*ndim+0]),2);
            fval[j] = exp(x[j*ndim+0]) * gsl_so1D_eval( fdata_in.qess_so, exp(x[j*ndim+0]) ) * gsl_so1D_eval( fdata_in.sync_x_so, t ) * beta2;
        }
        return 0;
    }

/* Integrate the synchrotron kernel over a solved electron population at GeV photon energy. */
static double eps_SY_4( double E_gam__GeV, double B__G, gsl_spline_object_1D sync_x_so, gsl_spline_object_1D qess_so )
{

    double res = 0.;
    double abserr;
    struct cg_eps_SY_4_fdata_sync fdata;

    /* x uses the critical energy for pitch angle pi/2; R averages the angles. */
    fdata.xE2 = (4.*M_PI*pow(m_e__g,2)*pow(c__cmsm1,3))/(3.*e__esu*B__G*h__ergs) * E_gam__GeV*m_e__GeV;
    fdata.qess_so = qess_so;
    fdata.sync_x_so = sync_x_so;

    double xmin[1], xmax[1];

    xmin[0] = log(E_CRe_lims__GeV[0]);
    xmax[0] = log(E_CRe_lims__GeV[1]);

    hcubature_v( 1, cg_eps_SY_4_F_SY, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

    return (sqrt(3.) * pow(e__esu,3) * B__G)/(h__ergs * m_e__g * pow(c__cmsm1,2)) * res/E_gam__GeV;
}
