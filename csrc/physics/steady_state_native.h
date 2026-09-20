/**
 * @file steady_state_native.h
 * @brief Finite-volume assembly and steady-state electron solves.
 *
 * Defines transition, injection, escape and continuous-loss integrals used
 * to assemble the disc/halo kinetic equations for two electron populations.
 * Uses checked worker-owned allocations and the runtime linear solver;
 * reference energy-grid rounding is retained.
 *
 * Scientific worker kernels extracted from the improvised legacy code.
 * Original author: Matt Roth. GPL-2.0; see repository LICENSE.
 * No lookup-table generation, file I/O or Python callbacks live here.
 */

/* Evaluate transition-rate moments for the lower-boundary electron matrix contribution. */
static int F_Gamma_i0_2D_2_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        if (exp(x[j*ndim+0]) > exp(x[j*ndim+1]))
        {
            fval[j * fdim + 0] = exp(x[j*ndim+0]) * dGammadlogEf_total__logGeVm1sm1( exp(x[j*ndim+0]),
                                 exp(x[j*ndim+1]), fdata_in.n_gso2D, fdata_in.gso_2D_radfield, fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
            fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
        }
        else
        {
            fval[j * fdim + 0] = 0.;
            fval[j * fdim + 1] = 0.;
        }
    }
    return 0;
}

/* Evaluate four log-energy moments of nonlocal electron transitions for matrix assembly. */
static int F_Gamma_2D_4_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    double dummy1;
    for (j = 0; j < npts; ++j)
    {
        if (exp(x[j*ndim+0]) > exp(x[j*ndim+1]))
        {
            dummy1 = dGammadlogEf_total__logGeVm1sm1( exp(x[j*ndim+0]), exp(x[j*ndim+1]),
                                                      fdata_in.n_gso2D, fdata_in.gso_2D_radfield, fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
            fval[j * fdim + 0] = exp(x[j*ndim+0]) * dummy1;
            fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
            fval[j * fdim + 2] = exp(x[j*ndim+1]) * dummy1;
            fval[j * fdim + 3] = x[j*ndim+0] * fval[j * fdim + 2];
        }
        else
        {
            fval[j * fdim + 0] = 0.;
            fval[j * fdim + 1] = 0.;
            fval[j * fdim + 2] = 0.;
            fval[j * fdim + 3] = 0.;
        }
    }
    return 0;
}

/* Evaluate escape and continuous-loss moments for the electron matrix. */
static int F_EdotDE_2_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j * fdim + 0] = exp(x[j*ndim+0]) * gsl_so1D_eval( fdata_in.gso_1D_D__cm2sm1, exp(x[j*ndim+0]) )/pow(fdata_in.h__pc*pc__cm,2) -
                             fdata_in.E_func( exp(x[j*ndim+0]), fdata_in.B__G, fdata_in.n_H__cmm3, fdata_in.h__pc );
        fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
    }
    return 0;
}

/* Evaluate the energy-squared injection term in logarithmic electron-energy coordinates. */
static int F_QE2_i_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j] = pow(exp(x[j*ndim+0]),2) * gsl_so1D_eval( fdata_in.gso_1D_Q, exp(x[j*ndim+0]) );
    }
    return 0;
}

/* Assemble and solve the coupled primary/secondary steady-state electron equations. */
static int CRe_steadystate_solve( int structure, double E_e_lim__GeV[2], int n_E, double n_H__cmm3, double B__G, double h__pc,
    unsigned int n_gso2D, gsl_spline_object_2D * gso_2D_radfields, gsl_spline_object_2D gso2D_BS, gsl_spline_object_1D gso_1D_D__cm2sm1,
    gsl_spline_object_1D gso_1D_Q_inject_1, gsl_spline_object_1D gso_1D_Q_inject_2,
    gsl_spline_object_1D * qe_1_so_1D, gsl_spline_object_1D * qe_2_so_1D )
{

    if(n_E<4 || n_E>500) cg_fail(CG_INVALID);
    double E__GeV[n_E+1];
    logspace_array( n_E+1, E_e_lim__GeV[0], E_e_lim__GeV[1], E__GeV );

    double DeltalogE = log(E_e_lim__GeV[1]/E_e_lim__GeV[0])/n_E;

    int i,j;

    double lnE_i__GeV[n_E];
    for (i = 0; i < n_E; ++i)
    {
        lnE_i__GeV[i] = (log(E__GeV[i])+log(E__GeV[i+1]))/2.;
    }

    double E_out__GeV[n_E+2];
    for (i = 0; i < n_E; ++i)
    {
        E_out__GeV[i+1] = exp((log(E__GeV[i])+log(E__GeV[i+1]))/2.);
    }
    E_out__GeV[0] = E__GeV[0];
    E_out__GeV[n_E+1] = E__GeV[n_E];

    double xmin[1], xmax[1];
    double xmin2D[2], xmax2D[2];
    double res, abserr;
    double res2[2], abserr2[2];
    double res4[4], abserr4[4];

    struct F_int_data fdata;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.B__G = B__G;
    fdata.h__pc = h__pc;
    fdata.n_gso2D = n_gso2D;
    fdata.gso_2D_radfield = gso_2D_radfields;
    fdata.gso2D_BS = gso2D_BS;
    fdata.gso_1D_D__cm2sm1 = gso_1D_D__cm2sm1;

    double Gamma_i0[n_E];
    double Gamma_i0_prime[n_E];
    xmin2D[1] = log(m_e__GeV);
    xmax2D[1] = log(E__GeV[0]);
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_Gamma_i0_2D_2_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        Gamma_i0[i] = res2[0]/DeltalogE;
        Gamma_i0_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_i0[i];
    }

    double Gamma_ii[n_E];
    double Gamma_ii_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        xmin2D[1] = log(E__GeV[i]);
        xmax2D[1] = log(E__GeV[i+1]);
        hcubature_v( 4, F_Gamma_2D_4_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res4, abserr4 );
        Gamma_ii[i] = (res4[0]/DeltalogE) - (res4[2]/DeltalogE) ;
        Gamma_ii_prime[i] = (res4[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * (res4[0]/DeltalogE)) -
                            (res4[3]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * (res4[2]/DeltalogE));
    }

    double **Gamma_ij = malloc(sizeof *Gamma_ij * n_E);
    if (Gamma_ij){for (i = 0; i < n_E; i++){Gamma_ij[i] = malloc(sizeof *Gamma_ij[i] * n_E);}}
    double **Gamma_ij_prime = malloc(sizeof *Gamma_ij_prime * n_E);
    if (Gamma_ij_prime){for (i = 0; i < n_E; i++){Gamma_ij_prime[i] = malloc(sizeof *Gamma_ij_prime[i] * n_E);}}
    double **Gamma_ji_prime = malloc(sizeof *Gamma_ji_prime * n_E);
    if (Gamma_ji_prime){for (i = 0; i < n_E; i++){Gamma_ji_prime[i] = malloc(sizeof *Gamma_ji_prime[i] * n_E);}}
    double **Gamma_ji = malloc(sizeof *Gamma_ji * n_E);
    if (Gamma_ji){for (i = 0; i < n_E; i++){Gamma_ji[i] = malloc(sizeof *Gamma_ji[i] * n_E);}}

    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        for (j = 0; j < i; ++j)
        {
            xmin2D[1] = log(E__GeV[j]);
            xmax2D[1] = log(E__GeV[j+1]);

            hcubature_v( 4, F_Gamma_2D_4_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res4, abserr4 );

            Gamma_ij[i][j] = res4[0]/DeltalogE;
            Gamma_ij_prime[i][j] = res4[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_ij[i][j];
            Gamma_ji[i][j] = res4[2]/DeltalogE;
            Gamma_ji_prime[i][j] = res4[3]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_ji[i][j];
        }
        for (j = i; j < n_E; ++j)
        {
            Gamma_ij[i][j] = 0.;
            Gamma_ij_prime[i][j] = 0.;
            Gamma_ji[i][j] = 0.;
            Gamma_ji_prime[i][j] = 0.;
        }
    }

    double Gamma_i[n_E];
    double Gamma_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        Gamma_i[i] = Gamma_i0[i] + Gamma_ii[i];
        Gamma_i_prime[i] = Gamma_i0_prime[i] + Gamma_ii_prime[i];
        for (j = 0; j < i; ++j)
        {
            Gamma_i[i] += Gamma_ij[i][j];
            Gamma_i_prime[i] += Gamma_ij_prime[i][j];
        }
    }

    double Edot_i[n_E+1];
    if (structure == 1)
    {
        fdata.E_func = dEdtm1_total_disc__GeVsm1;
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_disc__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc );
        }
    }
    else if (structure == 2)
    {
        fdata.E_func = dEdtm1_total_halo__GeVsm1;
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_halo__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc );
        }
    }
    else
    {
        printf("No valid structure specified in CRe steady state solver!");
    }

    double D_i[n_E];
    double D_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_EdotDE_2_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        D_i[i] = res2[0]/DeltalogE;
        D_i_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * D_i[i];
    }

    double **M = malloc(sizeof *M * n_E);
    if (M){for (i = 0; i < n_E; i++){M[i] = malloc(sizeof *M[i] * n_E);}}

    for (i = 0; i < n_E-2; ++i)
    {

        for (j = 0; j < i-1; ++j)
        {
            M[i][j] = 0.;
        }

        if (i > 0)
        {
            M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;
        }

        if (i == 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i]
                      - Gamma_ji_prime[i+1][i]/2. - Edot_i[i]/(2.*DeltalogE) + D_i_prime[i] + Gamma_i_prime[i];
        }
        else if (i > 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] - Gamma_ji_prime[i+1][i]/2.;
        }

        M[i][i+1] = - Edot_i[i+1]/DeltalogE + Edot_i[i]/(4.*DeltalogE) - D_i_prime[i]/2.
                    + Gamma_ji[i+1][i] - Gamma_i_prime[i]/2. - Gamma_ji_prime[i+2][i]/2.;

        if (i == n_E-3)
        {
            M[i][i+2] = - Edot_i[i+1]/(4.*DeltalogE) + Gamma_ji[i+2][i] + Gamma_ji_prime[i+1][i]/2.;
        }
        else if (i < n_E-3)
        {
            M[i][i+2] = - Edot_i[i+1]/(4.*DeltalogE) + Gamma_ji[i+2][i] + Gamma_ji_prime[i+1][i]/2. - Gamma_ji_prime[i+3][i]/2.;
        }

        for (j = i+3; j < n_E; ++j)
        {
            if (j == n_E-1)
            {
                M[i][j] = Gamma_ji[j][i] + Gamma_ji_prime[j-1][i]/2.;
            }
            else if (j < n_E-1)
            {
                M[i][j] = Gamma_ji[j][i] + Gamma_ji_prime[j-1][i]/2. - Gamma_ji_prime[j+1][i]/2.;
            }

        }
    }

    i = n_E-2;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;

    M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] - Gamma_ji_prime[i+1][i]/2.;

    M[i][i+1] = - Edot_i[i+1]/DeltalogE + Edot_i[i]/(4.*DeltalogE) - D_i_prime[i]/2.
                + Gamma_ji[i+1][i] - Gamma_i_prime[i];

    i = n_E-1;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;

    M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i];

    double Q_i_1[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_1;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_1[i] = -1.*res/DeltalogE;
    }

    double Q_i_2[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_2;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_2[i] = -1.*res/DeltalogE;
    }

    double q_e_1[n_E+2];

    if (Q_i_1[n_E-1] == 0.)
    {
        int n_Estar1 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_1[i] == 0.)
            {
                n_Estar1 = n_Estar1 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_1star[n_Estar1];
        for (i = 0; i < n_Estar1; ++i)
        {
            Q_i_1star[i] = Q_i_1[i];
        }
        double *A1 = malloc(sizeof A1 * n_Estar1*n_Estar1);
        for (i = 0; i < n_Estar1; ++i)
        {
            for (j = 0; j < n_Estar1; ++j)
            {
                A1[i*n_Estar1+j] = M[i][j];
            }
        }
        double x_out_1[n_Estar1];
        cg_linear_solve( n_Estar1, A1, Q_i_1star, x_out_1);
        for (i = 0; i < n_Estar1; ++i)
        {
            q_e_1[i+1] = x_out_1[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar1; i < n_E; ++i)
        {
            q_e_1[i+1] = 0.;
        }
        free( A1 );

    }
    else
    {
        double *A1 = malloc(sizeof A1 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A1[i*n_E+j] = M[i][j];
            }
        }
        double x_out_1[n_E];
        cg_linear_solve( n_E, A1, Q_i_1, x_out_1);
        for (i = 0; i < n_E; ++i)
        {
            q_e_1[i+1] = x_out_1[i]/E_out__GeV[i+1];
        }
        free( A1 );
    }

    double q_e_2[n_E+2];

    if (Q_i_2[n_E-1] == 0.)
    {
        int n_Estar2 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_2[i] == 0.)
            {
                n_Estar2 = n_Estar2 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_2star[n_Estar2];
        for (i = 0; i < n_Estar2; ++i)
        {
            Q_i_2star[i] = Q_i_2[i];
        }
        double *A2 = malloc(sizeof A2 * n_Estar2*n_Estar2);
        for (i = 0; i < n_Estar2; ++i)
        {
            for (j = 0; j < n_Estar2; ++j)
            {
                A2[i*n_Estar2+j] = M[i][j];
            }
        }
        double x_out_2[n_Estar2];
        cg_linear_solve( n_Estar2, A2, Q_i_2star, x_out_2);
        for (i = 0; i < n_Estar2; ++i)
        {
            q_e_2[i+1] = x_out_2[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar2; i < n_E; ++i)
        {
            q_e_2[i+1] = 0.;
        }
        free( A2 );
    }
    else
    {
        double *A2 = malloc(sizeof A2 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A2[i*n_E+j] = M[i][j];
            }
        }
        double x_out_2[n_E];
        cg_linear_solve( n_E, A2, Q_i_2, x_out_2);
        for (i = 0; i < n_E; ++i)
        {
            q_e_2[i+1] = x_out_2[i]/E_out__GeV[i+1];
        }
        free( A2 );
    }

    q_e_1[0] = fmax(0.,exp( ((log(q_e_1[2])-log(q_e_1[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_1[1]) ));
    q_e_1[n_E+1] = fmax(0.,exp( ((log(q_e_1[n_E])-log(q_e_1[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_1[n_E]) ));
    q_e_2[0] = fmax(0.,exp( ((log(q_e_2[2])-log(q_e_2[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_2[1]) ));
    q_e_2[n_E+1] = fmax(0.,exp( ((log(q_e_2[n_E])-log(q_e_2[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_2[n_E]) ));

    *qe_1_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_1 );
    *qe_2_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_2 );

    free2D( n_E, Gamma_ij );
    free2D( n_E, Gamma_ij_prime );
    free2D( n_E, Gamma_ji );
    free2D( n_E, Gamma_ji_prime );

    free2D( n_E, M );

    return 0;

}
