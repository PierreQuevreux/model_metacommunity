#include "metacommunity.h"
#include "community.h"
#include "species.h"
#include "functions.h"

#include <gsl/gsl_matrix.h> // to use matrix and vectors
#include <gsl/gsl_linalg.h> // to use gsl_linalg_solve_tridiag
#include <gsl/gsl_multiroots.h> // to use the root finder algorithm
#include <gsl/gsl_eigen.h> // to compute the eigenvalues of the system
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h> // for integration
// For random numbers
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int Metacommunity::get_dim(){
    return(_dim);
}
void Metacommunity::set_object(int seed){
    _species = new Species*[_dim];
        for (int i=0; i<_dim; i++){
            _species[i] = NULL;
        }
    _community = new Community*[_nCommunity];
        for (int i=0; i<_nCommunity; i++){
            _community[i] = NULL;
        }
    _dispersal = new double*[_dim]; // create the table containing the dispersal coefficients
    _asymmetry = new double*[_nCommunity]; // create the table containing the dispersal coefficients
    _perturbation = gsl_matrix_calloc (_dim_prtb, _dim_prtb); // create the covariance matrix of the perturbations
    //_response_press = gsl_matrix_calloc (_dim, _dim);
    _covariance = gsl_vector_calloc (_dim*_dim);
    _correlation = gsl_vector_calloc (_dim*_dim);
    // Random objects
    const gsl_rng_type* T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    _rng = gsl_rng_alloc(T); // instances the random  umber generator
    gsl_rng_set(_rng, seed); // change the seed of the generator to have different values after each run

    // Biomass of each species
    _t_extinction = new double[_dim]; // value at t
    _biomass_species = new double[_dim]; // value at t
    _biomass_species_equilibrium = new double[_dim]; // sum of values
    _biomass_species_mean = new double[_dim]; // sum of values
    _biomass_species_SQ = new double[_dim]; // sum of square
    _biomass_species_CV = new double[_dim]; // coefficient of variation
    _balance = new double[_dim]; // ratio of dispersal process to demographic process
    _balance_plus = new double[_dim]; // ratio of dispersal process to demographic process (positive terms)
    _balance_minus = new double[_dim]; // ratio of dispersal process to demographic process (negative terms)
    for (int i=0; i<_dim; i++){
        _t_extinction[i] = 0;
        _biomass_species[i] = 0;
        _biomass_species_equilibrium[i] = 0;
        _biomass_species_mean[i] = 0;
        _biomass_species_SQ[i] = 0;
        _biomass_species_CV[i] = 0;
        _balance[i] = 0;
        _balance_plus[i] = 0;
        _balance_minus[i] = 0;
        }
    // Biomass of each species across communities
    _biomass_species_total = new double[_nSpecies]; // value at t
    _biomass_species_total_mean = new double[_nSpecies]; // sum of values
    _biomass_species_total_SQ = new double[_nSpecies]; // sum of square
    _biomass_species_total_CV = new double[_nSpecies]; // coefficient of variation
    for (int i=0; i<_nSpecies; i++){
        _biomass_species_total[i] = 0;
        _biomass_species_total_mean[i] = 0;
        _biomass_species_total_SQ[i] = 0;
        _biomass_species_total_CV[i] = 0;
        }
    // Total biomass of a community
    _biomass_community = new double[_nCommunity]; // value at t
    _biomass_community_mean = new double[_nCommunity]; // sum of values
    _biomass_community_SQ = new double[_nCommunity]; // sum of square
    _biomass_community_CV = new double[_nCommunity]; // coefficient of variation
    for (int i=0; i<_nCommunity; i++){
        _biomass_community[i] = 0;
        _biomass_community_mean[i] = 0;
        _biomass_community_SQ[i] = 0;
        _biomass_community_CV[i] = 0;
    }
    _time_series = new double*[_n_step_record];
    create_table(_time_series, _n_step_record, _dim+1);
}
void Metacommunity::set_dispersal_matrix(const string file_matrix){
    create_table(_dispersal, _dim, _dim); // create the dispersal matrix (empty)
    load_table(file_matrix, _dim, _dispersal, false); // load the dispersal matrix
}
void Metacommunity::set_perturbation_matrix(const string file_matrix){
    double** perturbation = new double* [_dim_prtb];
    create_table(perturbation,_dim_prtb,_dim_prtb);
    load_table(file_matrix, _dim_prtb, perturbation, false); // load the dispersal matrix
    for (int i=0; i<_dim_prtb; i++){
        for (int j=0; j<_dim_prtb; j++){
            gsl_matrix_set(_perturbation, i, j, perturbation[i][j]); // converts the table into a matrix
        }
    }
    delete_table(perturbation, _dim_prtb);
}
void Metacommunity::set_asymmetry_matrix(const string file_matrix, const int n_params_asymmetry){
    create_table(_asymmetry, _nCommunity, n_params_asymmetry);
    load_table(file_matrix, _nCommunity, _asymmetry, false); // load the dispersal matrix
}
void Metacommunity::set_community(const double m, const double g, const double r, const double D, const double e, const double a, const double FR, const double h){
    for (int i=0; i<_nCommunity; i++){
        _community[i] = new Community(_nSpecies, i); // create the communities
        _community[i]->set_object();
        _community[i]->set_species(m, g, r, D, e, a, FR, h,
                                   _asymmetry[i][0], // m
                                   _asymmetry[i][1], // g
                                   _asymmetry[i][2], // r
                                   _asymmetry[i][3], // D
                                   _asymmetry[i][4], // e
                                   _asymmetry[i][5], // a
                                   _asymmetry[i][6], // FR
                                   _asymmetry[i][7]); // h
        _community[i]->set_interaction();
        for (int j=0; j<_nSpecies; j++){
            _species[i*_nSpecies+j] = _community[i]->get_species_address(j); //
        }
    }
}
// INTEGRATION
// Jacobian matrix
void Metacommunity::Derivative(const double* y, double* f){
    double prey(0); // effects of prey on growth
    double pred(0); // effects of predators on growth
    double disp(0); // effect of dispersal
    for (int i=0; i<_dim; i++){
        prey = 0; pred=0; disp=0; // reinitialisation
        for (int j=0; j<_species[i]->get_nPrey(); j++){
            prey += _species[i]->get_e() * _species[i]->get_a() * pow(y[_species[i]->get_preyID(j)],_species[i]->get_FR()) /
            (1 + _species[i]->get_a() * _species[i]->get_h() * pow(y[_species[i]->get_preyID(j)],_species[i]->get_FR()) );
        }
        for (int j=0; j<_species[i]->get_nPred(); j++){
            pred -= _species[i]->get_m() * _species[_species[i]->get_predID(j)]->get_a() * y[_species[i]->get_predID(j)] * pow(y[i],_species[_species[i]->get_predID(j)]->get_FR()-1) /
            (1 + _species[_species[i]->get_predID(j)]->get_a() * _species[_species[i]->get_predID(j)]->get_h() * pow(y[i],_species[_species[i]->get_predID(j)]->get_FR()) );
        }
        for (int j=0; j<_dim; j++){
            disp += _dispersal[i][j] * y[j]; // effect of dispersal
            //disp += _dispersal[i][j] * pow(y[j],2); // effect of dispersal with the square of biomass
        }
        f[i] = pow(_species[i]->get_m(),_species[i]->get_IDc()) * _species[i]->get_D() * (y[i] * (_species[i]->get_g() / _species[i]->get_D() - y[i] + prey + pred) + disp);
        //cout << (_species[i]->get_g() / _species[i]->get_D() - y[i] + prey + pred) << " - ";
        //cout << pow(_species[i]->get_m(),_species[i]->get_IDc()) << " - ";
        //f[i] = pow(_species[i]->get_m(),i) * _species[i]->get_D() * y[i] * (_species[i]->get_g() / _species[i]->get_D() - y[i] + prey + pred); // no dispersal
        //f[i] = _species[i]->get_D() / y[i] * (_species[i]->get_g() / _species[i]->get_D() - y[i] + prey + pred + disp); // rescaled equation
    }
    //cout << endl;
}
void Metacommunity::Balance(const double* y){
    for (int i=0; i<_dim; i++){
        _balance[i] = 0;
        _balance_plus[i] = 0;
        _balance_minus[i] = 0;
    }
    double disp_plus(0); // effect of immigration
    double disp_minus(0); // effect of emigration
    for (int i=0; i<_dim; i++){
        disp_plus=0; disp_minus=0; // reinitialisation
        for (int j=0; j<_species[i]->get_nPrey(); j++){ // effect of prey
            _balance_plus[i] += abs(_species[i]->get_e() * _species[i]->get_a() * y[_species[i]->get_preyID(j)] * y[i]);
        }
        for (int j=0; j<_species[i]->get_nPred(); j++){ // effect of pred
            _balance_minus[i] += abs(_species[i]->get_m() * _species[i]->get_a() * y[_species[i]->get_predID(j)] * y[i]);
        }
        _balance_minus[i] += abs(_species[i]->get_g()/_species[i]->get_D() * y[i]); // basal growth/mortality rate
        _balance_minus[i] += pow(y[i],2); // self-regulation
        for (int j=0; j<_dim; j++){
            if (j==i){
                disp_minus = abs(_dispersal[i][j] * y[j]); // effect of dispersal
            }
            else {
                disp_plus += abs(_dispersal[i][j] * y[j]); // effect of dispersal
            }
        }
        _balance[i] = (disp_plus + disp_minus)/(disp_plus + disp_minus + _balance_plus[i] + _balance_minus[i]);
        _balance_plus[i] = disp_plus / (disp_plus + _balance_plus[i]);
        _balance_minus[i] = disp_minus / (disp_minus + _balance_minus[i]);
    }
}
void Metacommunity::Jacobian(const gsl_vector* x, gsl_matrix* jacobian){
    double J(0); // element of the Jacobian matrix
    gsl_matrix_set_zero(jacobian); // reset the matrix to zero
    for (int i=0; i<_dim; i++){
        J = _species[i]->get_g()/_species[i]->get_D() -2*gsl_vector_get(x, i);
        for (int j=0; j<_species[i]->get_nPrey(); j++){
            J += _species[i]->get_e() * _species[i]->get_a() * gsl_vector_get(x, _species[i]->get_preyID(j));
        }
        for (int j=0; j<_species[i]->get_nPred(); j++){
            J += -_species[i]->get_m() * _species[i]->get_a() * gsl_vector_get(x, _species[i]->get_predID(j));
        }
        gsl_matrix_set(jacobian, i, i, J);
        for (int j=0; j<_species[i]->get_nPrey(); j++){
            J = _species[i]->get_e() * _species[i]->get_a() * gsl_vector_get(x, i); // effects of prey on i
            gsl_matrix_set(jacobian, i, _species[i]->get_preyID(j), J);
        }
        for (int j=0; j<_species[i]->get_nPred(); j++){
            J = -_species[i]->get_m() * _species[i]->get_a() * gsl_vector_get(x, i); // effects of pred on i
            gsl_matrix_set(jacobian, i, _species[i]->get_predID(j), J);
        }
        for (int j=0; j<_dim; j++){
            J = pow(_species[i]->get_m(),_species[i]->get_IDc()) * _species[i]->get_D() * (gsl_matrix_get(jacobian, i, j) + _dispersal[i][j]); // add the effect of dispersal
            //J = pow(_species[i]->get_m(),_species[i]->get_IDc()) * _species[i]->get_D() * (gsl_matrix_get(jacobian, i, j) + 2 * _dispersal[i][j] * pow(gsl_vector_get(x, j),2)); // jacobian with square mass effects of dispersal
            //J = _species[i]->get_D() / gsl_vector_get(x,i) * (gsl_matrix_get(jacobian, i, j) + _dispersal[i][j]); // rescaled jacobian
            gsl_matrix_set(jacobian, i, j, J);
        }
    }
}
void Metacommunity::Interaction(gsl_matrix* jacobian){
    double J(0); // element of the Jacobian matrix
    gsl_matrix_set_zero(jacobian); // reset the matrix to zero
    for (int i=0; i<_dim; i++){
        J = -1;
        gsl_matrix_set(jacobian, i, i, J);
        for (int j=0; j<_species[i]->get_nPrey(); j++){
            J = _species[i]->get_e() * _species[i]->get_a(); // effects of prey on i
            gsl_matrix_set(jacobian, i, _species[i]->get_preyID(j), J);
        }
        for (int j=0; j<_species[i]->get_nPred(); j++){
            J = -_species[i]->get_m() * _species[i]->get_a(); // effects of pred on i
            gsl_matrix_set(jacobian, i, _species[i]->get_predID(j), J);
        }
    }
}
gsl_complex Metacommunity::Equilibrium(){
    // INITIALISATION
    // The equilibria without dispersal can be analytically computed by resolving Ax=b
    // Creation of vectors
    gsl_vector * diag = gsl_vector_calloc (_dim); // vector with the diagonal elements of the matrix A
    gsl_vector * abovediag = gsl_vector_calloc (_dim-1); // vector with the above sub-diagonal elements of the matrix A
    gsl_vector * belowdiag = gsl_vector_calloc (_dim-1); // vector with the below sub-diagonal elements of the matrix A
    gsl_vector * b = gsl_vector_calloc (_dim); // vector with the elements of b
    gsl_vector * x = gsl_vector_calloc (_dim); // vector with the biomass at equilibrium
    // Initialisation of vectors
    gsl_vector_set_all(diag, -1);
    int index(0);
    for (int i=0; i<_nSpecies-1; i++){
        for (int j=0; j<_nCommunity; j++){
            index = i+j*_nSpecies;
            gsl_vector_set(abovediag, index, -_species[index]->get_m() * _species[index]->get_a());
            gsl_vector_set(belowdiag, index, _species[index+1]->get_e() * _species[index+1]->get_a());
        }
    }
    for (int i=0; i<_dim; i++){
        gsl_vector_set(b, i, -_species[i]->get_g() / _species[i]->get_D());
    }
    // Compute the value of x (biomasses at equilibrium without dispersal)
    gsl_linalg_solve_tridiag(diag, abovediag, belowdiag, b, x);

    // ROOT-FINDER ALGORITHM
    size_t dimension = _dim;
    gsl_multiroot_function_fdf multiroot_functions = {&multiroot_f, &multiroot_df, &multiroot_fdf, dimension, this};

    const gsl_multiroot_fdfsolver_type * T = gsl_multiroot_fdfsolver_newton;
    gsl_multiroot_fdfsolver * s = gsl_multiroot_fdfsolver_alloc (T, dimension);

    gsl_multiroot_fdfsolver_set(s, &multiroot_functions, x);

    // ROOT FINDING
    int status(GSL_CONTINUE);
    int nIteration(0); // number of iterations
    while (status == GSL_CONTINUE && nIteration < 1000){
        nIteration++;
        status = gsl_multiroot_fdfsolver_iterate (s);
        if (status != GSL_SUCCESS)
        break;
        status = gsl_multiroot_test_residual (s->f, 1e-7);
    }
    for (int i=0; i<_dim; i++){
        _biomass_species_equilibrium[i] = gsl_vector_get(x, i); // stores the values
    }
    Balance(_biomass_species_equilibrium); // compute the ratio of dispersal process to demographic process

    // STABLE EQUILIBRIUM ?
    gsl_matrix * jacobian = gsl_matrix_calloc (_dim, _dim); // initialise the Jacobian matrix
    Jacobian(x, jacobian); // compute the final Jacobian matrix
    // Calculation of eigenvalues
    gsl_vector_complex *eval = gsl_vector_complex_alloc (_dim);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (_dim, _dim);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (_dim);
    gsl_eigen_nonsymmv (jacobian, eval, evec, w);
    double real[_dim];
    for (int i=0; i<_dim; i++){
        real[i] = GSL_REAL(gsl_vector_complex_get (eval, i)); // real part of the eigen values
    }
    //int index(0); // get the dominant eigen value
    for (int i=0; i<_dim; i++){
        if (real[i]>real[index]){
            index = i;
        }
    }
    gsl_complex eval_dom = gsl_vector_complex_get (eval, index);

    //FREE THE MEMORY
    gsl_vector_free (diag);
    gsl_vector_free (abovediag);
    gsl_vector_free (belowdiag);
    gsl_vector_free (b);
    //gsl_vector_free (x_0);
    gsl_vector_free (x);
    gsl_matrix_free (jacobian);
    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_eigen_nonsymmv_free (w);

    return(eval_dom);
}
void Metacommunity::Perturbation_equilibrium(){
    // RESPONSE TO A PRESS PERTURBATION
    gsl_vector* x = gsl_vector_calloc (_dim); // vector with the biomass at equilibrium
    gsl_matrix* jacobian = gsl_matrix_calloc (_dim, _dim); // Jacobian matrix
    for (int i=0; i<_dim; i++){
        gsl_vector_set(x, i, _biomass_species_equilibrium[i]);
    }
    //Interaction(jacobian); // compute the final Jacobian matrix
    //string monfichier = "results/interaction.txt";
    //save_matrix(jacobian, _dim, _dim, monfichier);
    Jacobian(x, jacobian); // compute the final Jacobian matrix
    //string monfichier = "results/jacobienne.txt";
    //save_matrix(jacobian, _dim, _dim, monfichier);
    _response_press = invert_matrix(jacobian, _dim); // invert the Jacobian matrix
    gsl_matrix_scale(_response_press, -1);

    // RESPONSE TO A STOCHASTIC PERTURBATION
    // defines T
    gsl_matrix* T = gsl_matrix_calloc (_dim, _dim_prtb); // matrix T
    for (int i=0; i<_dim; i++){
        gsl_matrix_set(T, i, i, pow(gsl_vector_get(x, i), 0)); // exogenous stochastic perturbation
        gsl_matrix_set(T, i, _dim+i, pow(gsl_vector_get(x, i), 0.5)); // demographic stochastic perturbation
        gsl_matrix_set(T, i, 2*_dim+i, pow(gsl_vector_get(x, i), 1)); // environmental stochastic perturbation
    }
    // effect of the white noise
    gsl_matrix* TV = gsl_matrix_calloc (_dim, _dim_prtb); // matrix TV
    gsl_matrix* TVT = gsl_matrix_calloc (_dim, _dim); // matrix TVT^T
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,T,_perturbation,0.0,TV); // compute TV
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,TV,T,0.0,TVT); // compute TVT^T
    gsl_vector* TVT_stack = column_stack(TVT, _dim, _dim); // turn TT into a vector stacking its columns

    // Kronecker product on the Jacobian matrix
    gsl_matrix* I = gsl_matrix_calloc (_dim, _dim); // identity matrix
    gsl_matrix_set_identity(I);
    gsl_matrix* kro = kronecker_product(jacobian, I, _dim, _dim, _dim, _dim); // Kronecker product JxI
    gsl_matrix* kro_bis = kronecker_product(I, jacobian, _dim, _dim, _dim, _dim); // Kronecker product IxJ
    gsl_matrix_add(kro, kro_bis); // JxI + IxJ
    //string monfichier = "results/kronecker.txt";
    //save_matrix(kro, _dim*_dim, _dim*_dim, monfichier);
    gsl_matrix* kro_inv = invert_matrix(kro, _dim*_dim); // invert the (JxI + IxJ) matrix
    // Final operation to compute the output cov matrix
    gsl_blas_dgemv(CblasNoTrans, -1.0, kro_inv, TVT_stack, 0.0, _covariance); // solve the Lyapunov equation
    // compute the correlation matrix
    double corr(0); // element of the correlation matrix
    for (int j=0; j<_dim; j++){
        for (int i=0; i<_dim; i++){
            corr = gsl_vector_get(_covariance, i+j*_dim) / pow(gsl_vector_get(_covariance, i+i*_dim) * gsl_vector_get(_covariance, j+j*_dim),0.5);
            gsl_vector_set(_correlation, i+j*_dim, corr);
        }
    }
    // FREE THE MEMORY
    gsl_vector_free (x);
    gsl_matrix_free (jacobian);
    gsl_matrix_free (T);
    gsl_matrix_free (TV);
    gsl_matrix_free (TVT);
    gsl_vector_free (TVT_stack);
    gsl_matrix_free (I);
    gsl_matrix_free (kro);
    gsl_matrix_free (kro_bis);
    gsl_matrix_free (kro_inv);
}
void Metacommunity::Dynamic_Euler_Maruyama(double t, const double t_step, const int n_step_trans, const int n_step_record, double t_step_record, const double extinction_threshold){
    double** noise = new double*[n_step_trans+n_step_record];
    create_table(noise, n_step_trans+n_step_record, _dim_prtb);

    double* y = new double[_dim]; // initial value
    double* f = new double[_dim]; // df/dt
    for (int i=0; i<_dim; i++){
        y[i] = _biomass_species_equilibrium[i];
        f[i] = 0;
    }
    for (int i=0; i<n_step_record; i++){
        for (int j=0; j<_dim+1; j++){
            _time_series[i][j] = 0;
        }
    }
    for (int i=0; i<n_step_trans+n_step_record; i++){
        for (int j=0; j<_dim_prtb; j++){
            noise[i][j] = sqrt(gsl_matrix_get(_perturbation, j, j)) * gsl_ran_gaussian(_rng, 1) * sqrt(t_step); // create the noise function
        }
    }
    double t_record(0); // time between two long recorded time steps
    _nPoint = 0; // initialisation of the number of recorded points counter for the calculation of CVs
    _nPoint_TS = 0; // initialisation of the number of recorded points counter

/////////////////
// INTEGRATION //
/////////////////

    for (int k=0; k<n_step_trans+n_step_record; k++){
        Derivative(y, f); // compute df/dt
        t = t + t_step;
        for (int i=0; i<_dim; i++){
            if (y[i]!=0){
                y[i] = y[i] + f[i] * t_step + noise[k][i] + pow(y[i], 0.5) * noise[k][_dim+i] + y[i] * noise[k][2*_dim+i];
                //y[i] = y[i] + f[i] * t_step ;
            }

            // Extinction
            if (y[i]<extinction_threshold && y[i]!=0){
                y[i]=0;
                _biomass_species[i] = 0;
                _biomass_species_mean[i] = 0;
                _biomass_species_SQ[i] = 0; // no more biomass
                _t_extinction[i] = t;
            }
        }
        if (k>=n_step_trans){

            _biomass_total = 0;
            // Biomass of each species
            for (int i=0; i<_dim; i++){
                _biomass_species[i] = y[i];
                _biomass_species_mean[i] += y[i];
                _biomass_species_SQ[i] += pow(y[i],2);
                // Total biomass
                _biomass_total += y[i];
            }
            // Biomass of each species across communities
            for (int i=0; i<_nSpecies; i++){
                _biomass_species_total[i] = 0;
                for (int j=0; j<_nCommunity; j++){
                    _biomass_species_total[i] += y[j*_nSpecies+i];
                }
                _biomass_species_total_mean[i] += _biomass_species_total[i];
                _biomass_species_total_SQ[i] += pow(_biomass_species_total[i],2);
            }
            // Total biomass of a community
            for (int i=0; i<_nCommunity; i++){
                _biomass_community[i] = 0;
                for (int j=0; j<_nSpecies; j++){
                    _biomass_community[i] += y[i*_nSpecies+j];
                }
                _biomass_community_mean[i] += _biomass_community[i];
                _biomass_community_SQ[i] += pow(_biomass_community[i],2);
            }
            // Total biomass
            _biomass_total_mean += _biomass_total;
            _biomass_total_SQ += pow(_biomass_total,2);

            // Time series
            t_record += t_step;
            if (t_record>t_step_record){
                _time_series[_nPoint_TS][0] = t;
                for (int j=0; j<_dim; j++){
                    _time_series[_nPoint_TS][j+1] = y[j];
                }
                _nPoint_TS++;
                t_record = 0; // reset t_record
            }
            _nPoint++;
        }
    }

////////////////////////
// VARIABLE COMPUTING //
////////////////////////

    for (int i=0; i<_dim; i++){
        _biomass_species_mean[i] /= _nPoint;
        if (_biomass_species_mean[i] > pow(10,-15)){
            _biomass_species_CV[i] = CV(_biomass_species_mean[i],_biomass_species_SQ[i],_nPoint);
        }
    }
    for (int i=0; i<_nSpecies; i++){
        _biomass_species_total_mean[i] /= _nPoint;
        if (_biomass_species_total_mean[i] > pow(10,-15)){
            _biomass_species_total_CV[i] = CV(_biomass_species_total_mean[i],_biomass_species_total_SQ[i],_nPoint);
        }
    }
    for (int i=0; i<_nCommunity; i++){
        _biomass_community_mean[i] /= _nPoint;
        if (_biomass_community_mean[i] > pow(10,-15)){
            _biomass_community_CV[i] = CV(_biomass_community_mean[i],_biomass_community_SQ[i],_nPoint);
        }
    }
    _biomass_total_mean /= _nPoint;
    if (_biomass_total_mean > pow(10,-15)){
        _biomass_total_CV = CV(_biomass_total_mean,_biomass_total_SQ,_nPoint); // coefficient of variation of primary production
    }

    delete_table(noise, n_step_trans+n_step_record);
    delete[] y;
    delete[] f;
}
// OUTPUT
void Metacommunity::output_analytic(const string parameter_values, ofstream& file_biomass_species_equilibrium, ofstream& file_press,
                                    ofstream& file_covariance, ofstream& file_correlation,
                                    ofstream& file_balance, ofstream& file_balance_plus, ofstream& file_balance_minus){
    // EQUILIBRIUM BIOMASS
    file_biomass_species_equilibrium << parameter_values;
    for (int i=0; i<_dim; i++){
        file_biomass_species_equilibrium << ";" << _biomass_species_equilibrium[i];
    }
    file_biomass_species_equilibrium << endl;

    // RESPONSE TO PRESS PERTURBATION
    file_press << parameter_values;
    for (int j=0; j<_dim; j++){
        for (int i=0; i<_dim; i++){
            file_press << ";" << gsl_matrix_get(_response_press, i, j) ; // element ij of the matrix
        }
    }
    file_press << endl;
    gsl_matrix_free(_response_press);

    // VARIANCE-COVARIANCE MATRIX
    file_covariance << parameter_values;
    for (int i=0; i<_dim*_dim; i++){
        file_covariance << ";" << gsl_vector_get(_covariance, i);
    }
    file_covariance << endl;

    // CORRELATION MATRIX
    file_correlation << parameter_values;
    for (int i=0; i<_dim*_dim; i++){
        file_correlation << ";" << gsl_vector_get(_correlation, i);
    }
    file_correlation << endl;

    // BALANCE BETWEEN DISPERSION PROCESS AND DEMOGRAPHIC PROCESS
    file_balance << parameter_values;
    for (int i=0; i<_dim; i++){
        file_balance << ";" << _balance[i];
    }
    file_balance << endl;

    file_balance_plus << parameter_values;
    for (int i=0; i<_dim; i++){
        file_balance_plus << ";" << _balance_plus[i];
    }
    file_balance_plus << endl;

    file_balance_minus << parameter_values;
    for (int i=0; i<_dim; i++){
        file_balance_minus << ";" << _balance_minus[i];
    }
    file_balance_minus << endl;
}
void Metacommunity::output_dynamic(const string parameter_values,
                                   ofstream& file_biomass_species_mean, ofstream& file_biomass_species_CV,
                                   ofstream& file_biomass_species_total_mean, ofstream& file_biomass_species_total_CV,
                                   ofstream& file_biomass_community_mean, ofstream& file_biomass_community_CV,
                                   ofstream& file_metacommunity, ofstream& file_time_series){
// FILE WITH SPECIES MEAN BIOMASS AND CV IN EACH COMMUNITY (PATCH)
    file_biomass_species_mean << parameter_values;
    file_biomass_species_CV << parameter_values;
    for (int i=0; i<_dim; i++){
        file_biomass_species_mean << ";" << _biomass_species_mean[i];
        file_biomass_species_CV << ";" << _biomass_species_CV[i];
    }
    file_biomass_species_mean << endl;
    file_biomass_species_CV << endl;

// FILE WITH SPECIES MEAN BIOMASS AND CV ACROSS ALL COMMUNITIES
    file_biomass_species_total_mean << parameter_values;
    file_biomass_species_total_CV << parameter_values;
    for (int i=0; i<_nSpecies; i++){
        file_biomass_species_total_mean << ";" << _biomass_species_total_mean[i];
        file_biomass_species_total_CV << ";" << _biomass_species_total_CV[i];
    }
    file_biomass_species_total_mean << endl;
    file_biomass_species_total_CV << endl;

// FILE WITH COMMUNITY MEAN BIOMASS AND CV
    file_biomass_community_mean << parameter_values;
    file_biomass_community_CV << parameter_values;
    for (int i=0; i<_nCommunity; i++){
        file_biomass_community_mean << ";" << _biomass_community_mean[i];
        file_biomass_community_CV << ";" << _biomass_community_CV[i];
    }
    file_biomass_community_mean << endl;
    file_biomass_community_CV << endl;

// METACOMMUNITY LEVELS
    file_metacommunity << ";" << parameter_values;
    file_metacommunity << ";" << _biomass_total_mean;
    file_metacommunity << _biomass_total_CV << endl;

// TIME SERIES
    for (int i=0; i<_nPoint_TS; i++){
        file_time_series << parameter_values;
        for (int j=0; j<_dim+1; j++){
            file_time_series << ";" << _time_series[i][j];
        }
        file_time_series << endl;
    }
}

Metacommunity::Metacommunity(int nSpecies,
                             int nCommunity,
                             int dim,
                             int n_step_record)
                             :_nSpecies(nSpecies)
                             ,_nCommunity(nCommunity)
                             ,_dim(dim)
                             ,_dim_prtb(3*dim)
                             ,_n_step_record(n_step_record)
                             ,_biomass_total(0)
                             ,_biomass_total_mean(0)
                             ,_biomass_total_SQ(0)
                             ,_biomass_total_CV(0)
                             ,_nPoint(0)
                             ,_nPoint_TS(0)
{}

Metacommunity::~Metacommunity(){
    delete[] _species;
    for (int i=0; i<_nCommunity; i++){
        delete _community[i];
    }
    delete[] _community;
    delete_table(_dispersal, _dim);
    delete_table(_asymmetry, _nCommunity);
    gsl_matrix_free(_perturbation);
    gsl_vector_free(_covariance);
    gsl_vector_free(_correlation);
    gsl_rng_free (_rng);
    delete[] _t_extinction;
    delete[] _biomass_species;
    delete[] _biomass_species_equilibrium;
    delete[] _biomass_species_mean;
    delete[] _biomass_species_SQ;
    delete[] _biomass_species_CV;
    delete[] _balance;
    delete[] _balance_plus;
    delete[] _balance_minus;
    delete[] _biomass_species_total;
    delete[] _biomass_species_total_mean;
    delete[] _biomass_species_total_SQ;
    delete[] _biomass_species_total_CV;
    delete[] _biomass_community;
    delete[] _biomass_community_mean;
    delete[] _biomass_community_SQ;
    delete[] _biomass_community_CV;
    delete_table(_time_series, _n_step_record);
}
