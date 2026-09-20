/**
 * @file hadronic.h
 * @brief Proton spectra, pion production and secondary-particle emission.
 *
 * Defines proton injection normalisation and proton-proton collision kernels,
 * including pion-decay photons, secondary electrons and neutrinos. Energy
 * arguments use GeV; cross-section and production-rate units are indicated
 * by the individual functions. Calorimetry is supplied by the galaxy worker.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */
/* Return the squared proton collision centre-of-mass energy in GeV squared. */
static double s( double T_p ){return 2. * m_p__GeV * ( T_p + 2. * m_p__GeV );}

/* Return the proton collision centre-of-mass Lorentz factor. */
static double gam_CM( double T_p ){return ( T_p + 2. * m_p__GeV )/( sqrt( s( T_p ) ) );}

/* Return the kinematic neutral-pion energy in the centre-of-mass frame, GeV. */
static double E_pi_CM( double T_p ){return ( s( T_p ) - 4. * pow(m_p__GeV, 2) + pow(m_pi0__GeV, 2) )/( 2. * sqrt( s( T_p ) ) );}

/* Return the neutral-pion centre-of-mass momentum in natural GeV units. */
static double P_pi_CM( double T_p ){return sqrt( pow(E_pi_CM( T_p ), 2) - pow(m_pi0__GeV, 2) );}

/* Return the centre-of-mass speed as a fraction of the speed of light. */
static double beta_CM( double T_p ){return sqrt( 1. - pow(gam_CM( T_p ), -2) );}

/* Transform the maximum neutral-pion energy to the laboratory frame, GeV. */
static double E_pi_LAB_max( double T_p ){return gam_CM( T_p ) * ( E_pi_CM( T_p ) + P_pi_CM( T_p ) * beta_CM( T_p ) );}

/* Evaluate the fitted inelastic proton-proton cross section in millibarns. */
static double sigma_inel__mb( double T_p ){
  if (T_p >= T_p_th){
    return ( 30.7 - 0.96 * log( T_p/T_p_th ) + 0.18 * pow(log( T_p/T_p_th ),2) ) * pow( (1. - pow( (T_p_th/T_p), 1.9 ) ), 3 );
    }
  else {return 0.;}
  }

/* Select the energy-dependent amplitude coefficient for pion gamma-ray production. */
static double cg_A_max_b1( double T_p ){if (T_p < 5.){return 9.53;} else {return 9.13;}}

/* Select the power-law coefficient in the pion gamma-ray amplitude fit. */
static double cg_A_max_b2( double T_p ){if (T_p < 5.){return 0.52;} else {return 0.35;}}

/* Select the logarithmic-curvature coefficient in the pion gamma-ray amplitude fit. */
static double cg_A_max_b3( double T_p ){if (T_p < 5.){return 0.054;} else {return 9.7e-3;}}

/* Evaluate the fitted single-pion production cross section in millibarns. */
static double cg_A_max_sigma_1pi_mb( double T_p ){
    double sigma_0 = 7.66e-3;
    double M_res = 1.1883;
    double Gam_res = 0.2264;
    double eta = sqrt( pow( s( T_p ) - pow(m_pi0__GeV, 2) - 4.*pow(m_p__GeV, 2), 2 ) - 16.*pow(m_pi0__GeV, 2)*pow(m_p__GeV, 2) )/( 2. *  m_pi0__GeV * sqrt( s( T_p ) ) );
    double gam_s = sqrt( pow(M_res,2) * (pow(M_res,2) + pow(Gam_res,2)) );
    double K = (sqrt(8.) * M_res * Gam_res * gam_s)/( M_PI * sqrt( pow(M_res,2) + gam_s ) );
    double f_BW = (m_p__GeV * K)/(pow( pow( sqrt( s( T_p ) ) - m_p__GeV, 2) - pow(M_res,2), 2 ) + pow(M_res,2) * pow(Gam_res,2));
    return sigma_0 * pow(eta,1.95) * ( 1. + eta + pow(eta,5) ) * pow(f_BW,1.86);
    }

/* Evaluate the fitted two-pion production cross section in millibarns. */
static double cg_A_max_sigma_2pi_mb( double T_p ){return 5.7 /( 1. + exp( -9.3 * ( T_p - 1.4) ) );}

/* Evaluate the low-energy pion multiplicity fit. */
static double cg_A_max_n_pi_low( double T_p ){
    double Q_p = (T_p - T_p_th)/m_p__GeV;
    return -6e-3 + 0.237 * Q_p - 0.023 * pow( Q_p, 2 );
    }

/* Evaluate the five-parameter pion multiplicity expression for a proton energy. */
static double cg_pion_multiplicity( double T_p, double a1, double a2, double a3, double a4, double a5 ){
      double xi_p = (T_p - 3.)/m_p__GeV;
      return a1 * pow( xi_p, a4 ) * ( 1. + exp( -a2 * pow( xi_p, a5 ) ) ) * ( 1. - exp( -a3 * pow( xi_p, 0.25 ) ) );
      }

/* Select the high-energy pion multiplicity fit for the proton energy range. */
static double cg_A_max_n_pi_high( double T_p ){

    if (T_p >= 5.){
      double a1 = 0.728;
      double a2 = 0.596;
      double a3 = 0.491;
      double a4 = 0.2503;
      double a5 = 0.117;
      return cg_pion_multiplicity( T_p, a1, a2, a3, a4, a5 );
      }
    else if (T_p > 50.){
      double a1 = 0.652;
      double a2 = 0.0016;
      double a3 = 0.488;
      double a4 = 0.1928;
      double a5 = 0.483;
      return cg_pion_multiplicity( T_p, a1, a2, a3, a4, a5 );
      }
    else if (T_p > 100.){
      double a1 = 5.436;
      double a2 = 0.254;
      double a3 = 0.072;
      double a4 = 0.075;
      double a5 = 0.166;
      return cg_pion_multiplicity( T_p, a1, a2, a3, a4, a5 );
      }
    else{
      return 0.;
      }
    }

/* Combine pion production branches into a cross section in millibarns. */
static double cg_A_max_sigma_pi_mb( double T_p ){
   if (T_p < T_p_th){return 0.;}
    else if (T_p < 2.){return cg_A_max_sigma_1pi_mb( T_p ) + cg_A_max_sigma_2pi_mb( T_p );}
    else if (T_p < 5.){return sigma_inel__mb( T_p ) * cg_A_max_n_pi_low( T_p );}
    else {return sigma_inel__mb( T_p ) * cg_A_max_n_pi_high( T_p );}
    }

/* Compute the cross-section enhancement relative to its 1-TeV value. */
static double cg_A_max_G( double T_p ){
    double T_p_0 = 1e3;
    return 1. + log(fmax( 1., sigma_inel__mb( T_p )/sigma_inel__mb( T_p_0 ) ));
    }

/* Calculate the nuclear enhancement factor used by the gamma-ray amplitude. */
static double cg_A_max_eps( double T_p ){
    double eps_c = 1.37;
    double eps_1 = 0.29;
    double eps_2 = 0.1;
    double sigma_pp_R_mb = 31.4;
    return eps_c + (eps_1 + eps_2) * (sigma_pp_R_mb * cg_A_max_G( T_p ))/sigma_inel__mb( T_p );
    }

/* Normalize the pion-decay differential gamma-ray cross section, including enhancement. */
static double A_max( double T_p ){
  double b0 = 5.9;

  double theta_p = T_p/m_p__GeV;

  double Amax;

  if (T_p < T_p_th){Amax = 0.;}
  if (T_p < 1.){Amax = b0 * cg_A_max_sigma_pi_mb( T_p )/E_pi_LAB_max( T_p ) * mb__cm2;}
  else {Amax = cg_A_max_b1( T_p ) * pow(theta_p, -cg_A_max_b2( T_p )) * exp( cg_A_max_b3( T_p ) * pow(log(theta_p), 2) ) * cg_A_max_sigma_pi_mb( T_p )/m_p__GeV * mb__cm2;}

  if ( T_p > m_pi0__GeV/K_pi ){return cg_A_max_eps( T_p ) * Amax;}
  else{return Amax;}
  }

/* Select the alpha shape parameter of the differential gamma-ray cross section. */
static double cg_dsig_dEg_alpha( double T_p ){if ( T_p < T_p_th ){return 0.;} if ( T_p <= 20.){return 1.0;} else {return 0.5;}}

/* Evaluate the kappa parameter controlling the gamma-ray cross-section shape. */
static double cg_dsig_dEg_kappa( double T_p ){return 3.29 - 0.2 * pow(T_p/m_p__GeV, -1.5);}

/* Evaluate the mu parameter used by the low-energy gamma-ray shape fit. */
static double cg_dsig_dEg_mu( double T_p ){double q = (T_p - 1.)/m_p__GeV; return 5./4. * pow(q, 5./4.) * exp( -5./4. * q );}

/* Select the beta shape parameter for the proton kinetic-energy range. */
static double cg_dsig_dEg_beta( double T_p ){
    if ( T_p < T_p_th ){return 0.;}
    if ( T_p <= 1.){return cg_dsig_dEg_kappa( T_p );}
    else if ( T_p <= 4.){return cg_dsig_dEg_mu( T_p ) + 2.45;}
    else if ( T_p <= 20.){return 1.5 * cg_dsig_dEg_mu( T_p ) + 4.95;}
    else if ( T_p <= 100.){return 4.2;}
    else {return 4.9;}
    }

/* Select the gamma shape parameter for the proton kinetic-energy range. */
static double cg_dsig_dEg_gamma( double T_p ){
    if ( T_p <= 1.){return 0.;}
    else if ( T_p <= 4.){return cg_dsig_dEg_mu( T_p ) + 1.45;}
    else if ( T_p <= 20.){return cg_dsig_dEg_mu( T_p ) + 1.5;}
    else {return 1.;}
    }

/* Evaluate the dimensionless gamma-ray spectral shape within its kinematic limits. */
static double cg_dsig_dEg_F( double T_p, double E_gam, double X_gam, double C ){
    if (pow(X_gam, cg_dsig_dEg_alpha( T_p )) <= 1.){
      if ( T_p <= 1.){return pow( (1. - X_gam) , cg_dsig_dEg_beta( T_p ) );}
      else {return pow( (1. - pow(X_gam, cg_dsig_dEg_alpha( T_p ))) , cg_dsig_dEg_beta( T_p ) ) / pow( (1. + X_gam/C) , cg_dsig_dEg_gamma( T_p ));}
      }
    else{return 0.;}
    }

/* Return the differential pion-decay photon cross section for GeV energies. */
static double dsig_dEg( double T_p, double E_gam ){
  double gam_pi0__GeV_LAB =  E_pi_LAB_max( T_p )/m_pi0__GeV;
  double beta_pi_LAB = sqrt( 1. - pow(gam_pi0__GeV_LAB, -2) );
  double E_gam_max = m_pi0__GeV/2. * gam_pi0__GeV_LAB * ( 1. + beta_pi_LAB );

  double Y_gam = E_gam + pow(m_pi0__GeV, 2)/(4. * E_gam);
  double Y_gam_max = E_gam_max + pow(m_pi0__GeV, 2)/(4. * E_gam_max);
  double X_gam = ( Y_gam - m_pi0__GeV )/( Y_gam_max - m_pi0__GeV );

  double lambda = 3.;
  double C = lambda * m_pi0__GeV/Y_gam_max;

  if (T_p < T_p_th){return 0.;}

  return A_max( T_p ) * cg_dsig_dEg_F( T_p, E_gam, X_gam, C );

  }

/* Evaluate the momentum-power-law particle spectrum with an exponential kinetic-energy cutoff. */
static double J( double T, double C, double q, double m, double T_cutoff )
{
    double beta = 1.0;
    double p = sqrt( pow(T,2) + 2.*m*T );

    return C * pow( p, -q ) * exp(-pow(T/T_cutoff, beta)) * p/T;

}

struct cg_C_norm_E_fdata_norm
    {
        double q;
        double m;
        double T_cutoff;
    };

/* Evaluate batched integration samples for injection energy normalization; state is passed explicitly. */
static int cg_C_norm_E_f( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_C_norm_E_fdata_norm fdata_in = *((struct cg_C_norm_E_fdata_norm *)fdata);
        for (j = 0; j < npts; ++j)
        {

            fval[j] = x[j*ndim+0] * J( x[j*ndim+0], 1., fdata_in.q, fdata_in.m, fdata_in.T_cutoff );
        }
        return 0;
    }

/* Integrate the energy-weighted injection shape to obtain its normalization. */
static double C_norm_E( double q, double m, double T_cutoff )
{
    double result;
    double abserr;

    struct cg_C_norm_E_fdata_norm fdata;
    fdata.q = q;
    fdata.m = m;
    fdata.T_cutoff = T_cutoff;

    double xmin[1] = {T_p_norm__GeV[0]};
    double xmax[1] = {T_p_norm__GeV[1]};

    hcubature_v( 1, cg_C_norm_E_f, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &result, &abserr );

    return result;
}

struct cg_eps_pi_fdata_PI
    {
        double E_gam__GeV;
        double n_H__cmm3;
        double T_p_cutoff__GeV;
        double C_p;
        gsl_spline_object_1D gso1D_fcal;
    };

/* Evaluate batched integration samples for calorimetry-weighted pion photon emission; state is
 * passed explicitly.
 */
static int cg_eps_pi_f( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_eps_pi_fdata_PI fdata_in = *((struct cg_eps_pi_fdata_PI *)fdata);
        double beta_p[npts], f_cal[npts];
        for (j = 0; j < npts; ++j)
        {
            beta_p[j] = sqrt( 1. - pow(m_p__GeV,2)/pow( x[j*ndim+0] + m_p__GeV, 2) );
            f_cal[j] = gsl_so1D_eval( fdata_in.gso1D_fcal, x[j*ndim+0] );
            fval[j] = dsig_dEg( x[j*ndim+0], fdata_in.E_gam__GeV ) * J( x[j*ndim+0], fdata_in.C_p, q_p_inject, m_p__GeV, fdata_in.T_p_cutoff__GeV ) *
                      c__cmsm1 * beta_p[j] * f_cal[j] * fdata_in.n_H__cmm3;
        }
        return 0;
    }

/* Integrate pion-decay photon emission with energy-dependent proton calorimetry. */
static double eps_pi( double E_gam__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal )
{
    double res;
    double abserr;

    struct cg_eps_pi_fdata_PI fdata;
    fdata.E_gam__GeV = E_gam__GeV;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.C_p = C_p;
    fdata.T_p_cutoff__GeV = T_p_cutoff__GeV;
    fdata.gso1D_fcal = gso1D_fcal;

    double xmin[1] = { T_CR_lims__GeV[0] };
    double xmax[1] = { T_CR_lims__GeV[1] };

    hcubature_v( 1, cg_eps_pi_f, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );

    return res;
}

struct cg_eps_pi_fcal1_fdata_PI
    {
        double E_gam__GeV;
        double n_H__cmm3;
        double T_p_cutoff__GeV;
        double C_p;
        gsl_spline_object_1D gso1D_fcal;
    };

/* Evaluate batched integration samples for full-calorimetry pion photon emission; state is passed
 * explicitly.
 */
static int cg_eps_pi_fcal1_f_fcal1( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_eps_pi_fcal1_fdata_PI fdata_in = *((struct cg_eps_pi_fcal1_fdata_PI *)fdata);
        double beta_p[npts];
        for (j = 0; j < npts; ++j)
        {
            beta_p[j] = sqrt( 1. - pow(m_p__GeV,2)/pow( x[j*ndim+0] + m_p__GeV, 2) );
            fval[j] = dsig_dEg( x[j*ndim+0], fdata_in.E_gam__GeV ) * J( x[j*ndim+0], fdata_in.C_p, q_p_inject, m_p__GeV, fdata_in.T_p_cutoff__GeV ) *
                      c__cmsm1 * beta_p[j] * fdata_in.n_H__cmm3;
        }
        return 0;
    }

/* Integrate the diagnostic pion-decay spectrum assuming complete proton calorimetry. */
static double eps_pi_fcal1( double E_gam__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal )
{
    double res;
    double abserr;

    struct cg_eps_pi_fcal1_fdata_PI fdata;
    fdata.E_gam__GeV = E_gam__GeV;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.C_p = C_p;
    fdata.T_p_cutoff__GeV = T_p_cutoff__GeV;
    fdata.gso1D_fcal = gso1D_fcal;

    double xmin[1] = { T_CR_lims__GeV[0] };
    double xmax[1] = { T_CR_lims__GeV[1] };

    hcubature_v( 1, cg_eps_pi_fcal1_f_fcal1, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );

    return res;
}

/* Return the unit step, taking the value one at zero. */
static double Heaviside( double arg ){
  if (arg < 0){return 0.;}
  else {return 1.;}
  }

/* Calculate pion production at GeV energy, returning zero outside calorimetry coverage. */
static double q_pi( double E_pi__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal )
{
    double T_p__GeV = E_pi__GeV/K_pi;

    if (T_p__GeV < gso1D_fcal.x_lim[0] || T_p__GeV > gso1D_fcal.x_lim[1])
    {
        return 0.;
    }

    double beta_p = sqrt( 1. - pow( m_p__GeV, 2 )/pow( T_p__GeV + m_p__GeV, 2 ) );
    double f_cal = gsl_so1D_eval( gso1D_fcal, T_p__GeV );
    return n_H__cmm3 * c__cmsm1/K_pi * beta_p * J( T_p__GeV, C_p, q_p_inject, m_p__GeV, T_p_cutoff__GeV ) * sigma_inel__mb( T_p__GeV ) * mb__cm2 * f_cal;
}

/* Evaluate the g_nu_mu branch of the dimensionless neutrino decay-energy distribution. */
static double cg_f_nu_mu2_g_nu_mu( double x )
    {
    double lambda = 1. - pow(m_mu__GeV/m_piC__GeV,2);
    double r = 1. - lambda;
    return (9. * pow(x,2) - 6. * log(x) - 4. * pow(x,3) - 5.) * (3. - 2. * r) / (9. * pow(1.-r,2));
    }

/* Evaluate the h_nu_mu1 branch of the dimensionless neutrino decay-energy distribution. */
static double cg_f_nu_mu2_h_nu_mu1( double x )
    {
    double lambda = 1. - pow(m_mu__GeV/m_piC__GeV,2);
    double r = 1. - lambda;
    return (9. * pow(r,2) - 6. * log(r) - 4. * pow(r,3) - 5.) * (3. - 2. * r) / (9. * pow(1.-r,2));
    }

/* Evaluate the h_nu_mu2 branch of the dimensionless neutrino decay-energy distribution. */
static double cg_f_nu_mu2_h_nu_mu2( double x )
    {
    double lambda = 1. - pow(m_mu__GeV/m_piC__GeV,2);
    double r = 1. - lambda;
    return (9. * (r + x) - 4. * (pow(r,2) + r*x + pow(x,2))) * (1. + 2. * r) * (r - x) / (9. * pow(r,2));
    }

/* Assemble the secondary muon-neutrino decay spectrum from its piecewise branches. */
static double f_nu_mu2( double x )
  {

  double lambda = 1. - pow(m_mu__GeV/m_piC__GeV,2);
  double r = 1. - lambda;
  return cg_f_nu_mu2_g_nu_mu(x) * Heaviside(x-r) + (cg_f_nu_mu2_h_nu_mu1(x) + cg_f_nu_mu2_h_nu_mu2(x)) * Heaviside(r-x);
  }

/* Evaluate the g_nu_e branch of the dimensionless neutrino decay-energy distribution. */
static double cg_f_nu_e_g_nu_e( double x )
    {

    double lambda = 1. - pow(m_mu__GeV/m_piC__GeV,2);
    double r = 1. - lambda;
    return 2. * ( (1.-x) * (6. * pow(1.-x,2) + r * (5. + 5. * x - 4. * pow(x,2))) + 6. * r * log(x))/(3. * pow(1.-r,2));

    }

/* Evaluate the h_nu_e1 branch of the dimensionless neutrino decay-energy distribution. */
static double cg_f_nu_e_h_nu_e1( double x )
    {
    double lambda = 1. - pow(m_mu__GeV/m_piC__GeV,2);
    double r = 1. - lambda;
    return 2. * ((1.-r) * (6. - 7. * r + 11. * pow(r,2) - 4. * pow(r,3)) + 6. * r * log(r))/(3. * pow(1.-r,2));
    }

/* Evaluate the h_nu_e2 branch of the dimensionless neutrino decay-energy distribution. */
static double cg_f_nu_e_h_nu_e2( double x )
    {
    double lambda = 1. - pow(m_mu__GeV/m_piC__GeV,2);
    double r = 1. - lambda;
    return 2. * (r - x) * (7. * pow(r,2) - 4. * pow(r,3) + 7. * x * r - 4. * x * pow(r,2) - 2. * pow(x,2) - 4. * pow(x,2) * r)/(3. * pow(r,2));
    }

/* Assemble the electron-neutrino decay spectrum from its piecewise branches. */
static double f_nu_e( double x )
  {

  double lambda = 1. - pow(m_mu__GeV/m_piC__GeV,2);
  double r = 1. - lambda;
  return cg_f_nu_e_g_nu_e(x) * Heaviside(x-r) + (cg_f_nu_e_h_nu_e1(x) + cg_f_nu_e_h_nu_e2(x)) * Heaviside(r-x);
  }

struct cg_q_nu_fdata_nu
    {
        double E_nu__GeV;
        double n_H__cmm3;
        double T_p_cutoff__GeV;
        double C_p;
        gsl_spline_object_1D gso1D_fcal;
    };

/* Evaluate batched integration samples for muon-decay neutrino production; state is passed explicitly. */
static int cg_q_nu_F_numu2_nue( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_q_nu_fdata_nu fdata_in = *((struct cg_q_nu_fdata_nu *)fdata);
        for (j = 0; j < npts; ++j)
        {
            fval[j] = 2. *  ( f_nu_e(x[j*ndim+0]) + f_nu_mu2(x[j*ndim+0]) ) * q_pi( fdata_in.E_nu__GeV/x[j*ndim+0],
                      fdata_in.n_H__cmm3, fdata_in.C_p, fdata_in.T_p_cutoff__GeV, fdata_in.gso1D_fcal )/x[j*ndim+0];
        }
        return 0;
    }

/* Evaluate batched integration samples for direct pion-decay muon-neutrino production; state is
 * passed explicitly.
 */
static int cg_q_nu_F_numu1( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_q_nu_fdata_nu fdata_in = *((struct cg_q_nu_fdata_nu *)fdata);
        double lambda = 1. - pow( m_mu__GeV/m_piC__GeV, 2 );

        for (j = 0; j < npts; ++j)
        {
            fval[j] = 2./lambda * q_pi( fdata_in.E_nu__GeV/x[j*ndim+0],
                      fdata_in.n_H__cmm3, fdata_in.C_p, fdata_in.T_p_cutoff__GeV, fdata_in.gso1D_fcal )/x[j*ndim+0];
        }
        return 0;
    }

/* Integrate the pion and muon decay contributions to the neutrino production spectrum. */
static double q_nu( double E_nu__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal )
{
    double res_numu2_nue, res_numu1;
    double abserr_numu2_nue, abserr_numu1;

    double xmin[1], xmax[1];

    struct cg_q_nu_fdata_nu fdata;
    fdata.E_nu__GeV = E_nu__GeV;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.C_p = C_p;
    fdata.T_p_cutoff__GeV = T_p_cutoff__GeV;
    fdata.gso1D_fcal = gso1D_fcal;

    xmin[0] = E_nu__GeV/( T_p_norm__GeV[1] * K_pi );
    xmax[0] = 1.;

    if (xmin[0] < xmax[0])
    {
        hcubature_v( 1, cg_q_nu_F_numu2_nue, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res_numu2_nue, &abserr_numu2_nue );
    }
    else
    {
        res_numu2_nue = 0.;
        abserr_numu2_nue = 0.;
    }

    double lambda = 1. - pow( m_mu__GeV/m_piC__GeV, 2 );

    xmin[0] = E_nu__GeV/( T_p_norm__GeV[1] * K_pi );
    xmax[0] = lambda;

    if (xmin[0] < xmax[0])
    {
        hcubature_v( 1, cg_q_nu_F_numu1, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res_numu1, &abserr_numu1 );
    }
    else
    {
        res_numu1 = 0.;
        abserr_numu1 = 0.;
    }

    return res_numu2_nue + res_numu1;
}

struct cg_q_e_fdata_qe
    {
        double E_e__GeV;
        double n_H__cmm3;
        double T_p_cutoff__GeV;
        double C_p;
        gsl_spline_object_1D gso1D_fcal;
    };

/* Evaluate batched integration samples for secondary electron injection; state is passed explicitly. */
static int cg_q_e_F_e( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct cg_q_e_fdata_qe fdata_in = *((struct cg_q_e_fdata_qe *)fdata);
        for (j = 0; j < npts; ++j)
        {
            fval[j] = fmax( 2. * f_nu_mu2(x[j*ndim+0])  * q_pi( fdata_in.E_e__GeV/x[j*ndim+0] , fdata_in.n_H__cmm3, fdata_in.C_p,
                      fdata_in.T_p_cutoff__GeV, fdata_in.gso1D_fcal )/x[j*ndim+0], 0.);
        }
        return 0;
    }

/* Integrate secondary electron injection at the requested kinetic energy in GeV. */
static double q_e( double T_e__GeV, double n_H__cmm3, double C_p, double T_p_cutoff__GeV, gsl_spline_object_1D gso1D_fcal )
{
    double res;
    double abserr;

    struct cg_q_e_fdata_qe fdata;
    fdata.E_e__GeV = T_e__GeV + m_e__GeV;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.C_p = C_p;
    fdata.T_p_cutoff__GeV = T_p_cutoff__GeV;
    fdata.gso1D_fcal = gso1D_fcal;

    double xmin[1] = { fdata.E_e__GeV/( T_p_norm__GeV[1] * K_pi) };
    double xmax[1] = { 1. };

    if (xmin[0] < xmax[0])
    {
        hcubature_v( 1, cg_q_e_F_e, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );
    }
    else
    {
        res = 0.;
        abserr = 0.;
    }
    return res;
}
