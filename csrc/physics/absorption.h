/**
 * @file absorption.h
 * @brief Internal gamma-ray absorption by photon-photon pair production.
 *
 * Defines the Breit-Wheeler cross section and integrates it over the galaxy's
 * target photon spectrum to obtain a dimensionless optical depth. Photon
 * energies are in GeV, cross sections in millibarns, and path lengths in pc.
 * These functions run inside native galaxy workers; EBL attenuation and
 * observer-frame conversion are handled separately in Python.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Evaluate the Breit-Wheeler photon-pair cross section in millibarns above threshold. */
static double sigma_gg_BW__mb( double E1__GeV, double E2__GeV )
{
    if ( E1__GeV * E2__GeV >= pow(m_e__GeV, 2) )
    {
        double beta_hat = sqrt( 1. - pow(m_e__GeV, 2)/( E1__GeV * E2__GeV ) );

        return 3./16. * sigma_T__mb * (1.-pow(beta_hat,2)) *
               (2.*beta_hat * (pow(beta_hat,2)-2.) + (3.-pow(beta_hat,4))*log((1.+beta_hat)/(1.-beta_hat)));
    }
    else
    {
        return 0.;
    }
}

struct cg_tau_gg_gal_BW_fdata_taugg
    {
        double E_gam__GeV;
        double (*n_phot)(double *, double);
        double *n_phot_params;
    };

/* Evaluate batched integration samples for internal photon-pair absorption; state is passed explicitly. */
static int cg_tau_gg_gal_BW_F_taugg( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_tau_gg_gal_BW_fdata_taugg fdata_in = *((struct cg_tau_gg_gal_BW_fdata_taugg *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] = exp(x[j*ndim+0]) * fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) *
                      sigma_gg_BW__mb( fdata_in.E_gam__GeV, exp(x[j*ndim+0]) );
        }
        return 0;
    }

/* Integrate internal pair-production optical depth over the target radiation spectrum. */
static double tau_gg_gal_BW( double E_gam__GeV, double (*n_phot)(double *, double), double *n_phot_params, double E_phot__GeV_lims[2], double h_pc )
{

    double res = 0.;
    double abserr;

    struct cg_tau_gg_gal_BW_fdata_taugg fdata;

    fdata.E_gam__GeV = E_gam__GeV;
    fdata.n_phot = n_phot;
    fdata.n_phot_params = n_phot_params;

    double xmin[1], xmax[1];

    xmin[0] = log(E_phot__GeV_lims[0]);
    xmax[0] = log(E_phot__GeV_lims[1]);

    hcubature_v( 1, cg_tau_gg_gal_BW_F_taugg, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );

    return h_pc * pc__cm * mb__cm2 * res;
}
