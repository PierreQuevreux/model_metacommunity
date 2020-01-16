#include "metacommunity.h"
#include "community.h"
#include "species.h"
#include "functions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "string"
#include "math.h"

#include <gsl/gsl_complex.h> // for complex numbers
#include <stdio.h> // for random numbers
#include <gsl/gsl_rng.h> // for random numbers

#include <omp.h> // OpenMP library for parallel calculation
// add the option -fopenmp to the compiler

using namespace std;

int main()
{
    int nSim(0); // number written on the file (starts at zero)
    string N; // string version of the int
    N = to_string(nSim); // string version of the int

    ////////////////////////////////
    // PARALLELISATION PARAMETERS //
    ////////////////////////////////
    string path("results/"); // path to the folder where results are saved
    string file_path(path + "parameters_data_" + N + ".txt"); // file containing the parameters
    bool header(true); // is there a header in the parameters file ?
    int n_start(0); // number of the first simulation (useful if the loop was interrupted)
    int n_simu(0); // number of simulations
    n_simu = load_parameters(file_path,1);
    int n_params(0); // number of parameters
    n_params = load_parameters(file_path,2);
    int nSpecies(0); // number of different species
    nSpecies = load_parameters(file_path,3);
    int nCommunity(0); // number of communities (i.e. patches)
    nCommunity = load_parameters(file_path,4);
    int dim(nSpecies*nCommunity); // dimension of the system
    int n_threads(0); // number of threads
    n_threads = load_parameters(file_path,5);
    omp_set_num_threads(n_threads); // number of threads used
    int n_params_asymmetry(8); // number of parameters subject to asymmetry between patches

    double **params; // table containing the parameters feeding the simulation
    params = new double*[n_simu]; // create the table
    create_table(params,n_simu,n_params); // create the table
    file_path = path + string("parameters_") + N + string(".txt"); // file containing the parameters
    string parameter_names = load_table(file_path,n_simu,params,header); // read the file containing the parameters

    /////////////////////////////
    // PROGRAM FOR SIMULATIONS //
    /////////////////////////////

    // PARAMETERS
    //double q(0); // scaling constant of the metabolic rate
    //double m(0.65); // predator/prey metabolic rate ratio
    //double g(1); // scaling constant of primary producers growth rate
    //double r(0); // scaling constant of consumers death rate
    //double D(1); // scaling constant of self-regulation
    //double e(0.65); // assimilation efficiency
    //double a(1/0.65); // scaling constant of consumption rate
    //double v(0); // exponent splitting the dependency of consumption rate between predator and prey metabolic rates

    // Integration parameters
    double t = 0.0, t_final = 100; // time span of integration // 350
    double t_record = 0; // time from witch recording begins // 50
    double t_step = 0.001; // time step
    double t_step_record = 0.1; // larger recorded time step (to reduce the number of points)
    //double t_fixed_step = 0.0001; // fixed time step for the ODE solving
    int n_step_trans(static_cast<int>((t_record-t)/t_step)+1); // total number of steps of the transitory dynamics
    int n_step_record(static_cast<int>((t_final-t_record)/t_step)+1); // number of recorded steps
    //cout << n_step_step << " - " << n_step_step*t_fixed_step << endl;
    //cout << n_tot_step << " - " << n_step << endl;
    //double h = 1e-6; // absolute accuracy
    double extinction_threshold (pow(10,-30)); // extinction biomass threshold

    // Seeding the noise random generator
    int seed[n_simu];
    const gsl_rng_type * T;
    gsl_rng * rng;
    T = gsl_rng_default;
    rng = gsl_rng_alloc (T);
    gsl_rng_env_setup();
    for (int i=0; i<n_simu; i++){
        seed[i] = gsl_rng_uniform_int (rng,n_simu*100);
    }
    gsl_rng_free (rng);

// GENERAL INFORMATION SFOR OUTPUT FILES
    string monfichier;
    //parameter_names= "simu_ID;dispersal;perturbation;stochastic_type;a;m;e;ma;ea"; // name of the parameters of the simulation
    string parameter_values; // values of parameters for output files

// SPECIES AND COMMUNITY LEVEL
    // FILE WITH SPECIES BIOMASS AT EQUILIBRIUM (FIXED POINTS)
    monfichier = path + "biomass_species_equilibrium_" + N + ".txt";
    ofstream file_biomass_species_equilibrium (monfichier.c_str()); // creation of the file
    file_biomass_species_equilibrium << parameter_names; // write the header
    for (int i=1; i<=nCommunity; i++){
        for (int j=1; j<=nSpecies; j++){
            file_biomass_species_equilibrium << ";x_" << j << "_" << i ; // each species of each patch
        }
    }
    file_biomass_species_equilibrium << endl;
    // FILE WITH SPECIES MEAN BIOMASS IN EACH COMMUNITY (PATCH)
    monfichier = path + "biomass_species_mean_" + N + ".txt";
    ofstream file_biomass_species_mean (monfichier.c_str()); // creation of the file
    file_biomass_species_mean << parameter_names; // write the header
    for (int i=1; i<=nCommunity; i++){
        for (int j=1; j<=nSpecies; j++){
            file_biomass_species_mean << ";x_" << j << "_" << i ; // each species of each patch
        }
    }
    file_biomass_species_mean << endl;

    // FILE WITH SPECIES BIOMASS CV IN EACH COMMUNITY (PATCH)
    monfichier = path + "biomass_species_CV_" + N + ".txt";
    ofstream file_biomass_species_CV (monfichier.c_str()); // creation of the file
    file_biomass_species_CV << parameter_names; // write the header
    for (int i=1; i<=nCommunity; i++){
        for (int j=1; j<=nSpecies; j++){
            file_biomass_species_CV << ";x_" << j << "_" << i ; // each species of each patch
        }
    }
    file_biomass_species_CV << endl;

// SPECIES LEVEL
    // FILE WITH SPECIES MEAN BIOMASS ACROSS ALL COMMUNITIES
    monfichier = path + "biomass_species_total_mean_" + N + ".txt";
    ofstream file_biomass_species_total_mean (monfichier.c_str()); // creation of the file
    file_biomass_species_total_mean << parameter_names; // write the header
    for (int i=1; i<=nSpecies; i++){
        file_biomass_species_total_mean << ";x" << i ; // each species
    }
    file_biomass_species_total_mean << endl;

    // FILE WITH SPECIES BIOMASS CV ACROSS ALL COMMUNITIES
    monfichier = path + "biomass_species_total_CV_" + N + ".txt";
    ofstream file_biomass_species_total_CV (monfichier.c_str()); // creation of the file
    file_biomass_species_total_CV << parameter_names; // write the header
    for (int i=1; i<=nSpecies; i++){
        file_biomass_species_total_CV << ";x" << i ; // each species
    }
    file_biomass_species_total_CV << endl;

// COMMUNITY LEVEL
    // FILE WITH COMMUNITY MEAN BIOMASS
    monfichier = path + "biomass_community_mean_" + N + ".txt";
    ofstream file_biomass_community_mean (monfichier.c_str()); // creation of the file
    file_biomass_community_mean << parameter_names; // write the header
    for (int i=1; i<=nCommunity; i++){
        file_biomass_community_mean << ";x" << i ; // each community
    }
    file_biomass_community_mean << endl;

    // FILE WITH COMMUNITY BIOMASS CV
    monfichier = path + "biomass_community_CV_" + N + ".txt";
    ofstream file_biomass_community_CV (monfichier.c_str()); // creation of the file
    file_biomass_community_CV << parameter_names; // write the header
    for (int i=1; i<=nCommunity; i++){
        file_biomass_community_CV << ";x" << i ; // each community
    }
    file_biomass_community_CV << endl;

// METACOMMUNITY LEVEL
    monfichier = path + "metacommunity_" + N + ".txt";
    ofstream file_metacommunity (monfichier.c_str()); // creation of the file
    file_metacommunity << parameter_names << ";biomass_total_mean;biomass_total_CV_"; // write the header
    file_metacommunity << endl;

// ANALYTICAL OUTPUTS
    // FILE WITH THE RESPONSE TO PRESS PERTURBATIONS
    monfichier = path + "press_" + N + ".txt";
    ofstream file_press (monfichier.c_str()); // creation of the file
    file_press << parameter_names; // write the header
    for (int i=1; i<=dim; i++){
        for (int j=1; j<=dim; j++){
            file_press << ";x_" << j << "_" << i ; // element ij of the press response matrix
        }
    }
    file_press << endl;

    // FILE WITH THE VARIANCE-COVARIANCE MATRIX
    monfichier = path + "covariance_" + N + ".txt";
    ofstream file_covariance (monfichier.c_str()); // creation of the file
    file_covariance << parameter_names; // write the header
    for (int i=1; i<=dim; i++){
        for (int j=1; j<=dim; j++){
            file_covariance << ";x_" << j << "_" << i ; // element ij of the matrix
        }
    }
    file_covariance << endl;

    // FILE WITH THE CORRELATION MATRIX
    monfichier = path + "correlation_" + N + ".txt";
    ofstream file_correlation (monfichier.c_str()); // creation of the file
    file_correlation << parameter_names; // write the header
    for (int i=1; i<=dim; i++){
        for (int j=1; j<=dim; j++){
            file_correlation << ";x_" << j << "_" << i ; // element ij of the matrix
        }
    }
    file_correlation << endl;

// FILE WITH TIME SERIES
    monfichier = path + "time_series_" + N + ".txt";
    ofstream file_time_series (monfichier.c_str()); // creation of the file
    file_time_series << parameter_names << ";t"; // write the header
    for (int i=1; i<=nCommunity; i++){
        for (int j=1; j<=nSpecies; j++){
            file_time_series << ";x_" << j << "_" << i ; // each species of each patch
        }
    }
    file_time_series << endl;

// FILE WITH BALANCE BETWEEN dispersal PROCESS AND DEMOGRAPHIC PROCESS
    monfichier = path + "balance_" + N + ".txt";
    ofstream file_balance (monfichier.c_str()); // creation of the file
    file_balance << parameter_names; // write the header
    for (int i=1; i<=nCommunity; i++){
        for (int j=1; j<=nSpecies; j++){
            file_balance << ";x_" << j << "_" << i ; // each species of each patch
        }
    }
    file_balance << endl;

    monfichier = path + "balance_plus_" + N + ".txt";
    ofstream file_balance_plus (monfichier.c_str()); // creation of the file
    file_balance_plus << parameter_names; // write the header
    for (int i=1; i<=nCommunity; i++){
        for (int j=1; j<=nSpecies; j++){
            file_balance_plus << ";x_" << j << "_" << i ; // each species of each patch
        }
    }
    file_balance_plus << endl;

    monfichier = path + "balance_minus_" + N + ".txt";
    ofstream file_balance_minus (monfichier.c_str()); // creation of the file
    file_balance_minus << parameter_names; // write the header
    for (int i=1; i<=nCommunity; i++){
        for (int j=1; j<=nSpecies; j++){
            file_balance_minus << ";x_" << j << "_" << i ; // each species of each patch
        }
    }
    file_balance_minus << endl;

// OUTPUT STRUCTURES

    /////////////////////////
    // LOOP TO PARALLELISE //
    /////////////////////////

    int count_simu(0); // counter of simulations
    #pragma omp parallel for // instruction to parallelise the following for loop
    for(int i=n_start; i<n_simu; i++){
        // Simulation
        Metacommunity MC(nSpecies,
                         nCommunity,
                         dim,
                         n_step_record); // stochastic type
        MC.set_object(seed[i]); // create the tables, matrix and vectors
        string N_matrix ; // string version of the int
        N_matrix = to_string((int)params[i][9]);
        string file_path_perturbation = path + string("matrix_perturbation_") + N_matrix + string(".txt"); // file containing the dispersal matrix
        MC.set_perturbation_matrix(file_path_perturbation); // load the perturbation matrix
        N_matrix = to_string((int)params[i][10]);
        string file_path_dispersal = path + string("matrix_dispersal_") + N_matrix + string(".txt"); // file containing the dispersal matrix
        MC.set_dispersal_matrix(file_path_dispersal); // load the dispersal matrix
        N_matrix = to_string((int)params[i][11]);
        string file_path_asymmetry = path + string("matrix_asymmetry_") + N_matrix + string(".txt"); // file containing the asymmetry between patches matrix
        MC.set_asymmetry_matrix(file_path_asymmetry, n_params_asymmetry); // load the asymmetry matrix
        MC.set_community(params[i][1], // m --- create the communities and the species
                         params[i][2], // g
                         params[i][3], // r
                         params[i][4], // D
                         params[i][5], // e
                         params[i][6], // a
                         params[i][7], // FR
                         params[i][8]); // h
        gsl_complex eigen; // largest eigenvalue of the system
        eigen = MC.Equilibrium(); // compute the equilibrium without dispersal
        //if (GSL_REAL(eigen)<0){
        MC.Perturbation_equilibrium(); // compute the response to a press and a stochastic perturbation
        //}
        //else {
            //MC.Dynamic(t, t_final, t_record, t_step, h, extinction_threshold, t_fixed_step); // solve numerically the SDE system
            //MC.Dynamic_fixed_step(t, t_final, t_record, t_step, t_fixed_step, h, extinction_threshold); // solve numerically the SDE system with a fixed time step
            //MC.Dynamic_Euler_Maruyama(t, t_step, n_step_trans, n_step_record, t_step_record, extinction_threshold); // solve the SDE with the Euler-Maruyama method
        //}

        #pragma omp critical
        {
            file_path = path + string("parameters_") + N + string(".txt"); // file containing the parameters
            parameter_values = output_parameter(file_path,i+1,header); // read ith line of the parameter file
            //if (GSL_REAL(eigen)<0){
                MC.output_analytic(parameter_values, file_biomass_species_equilibrium, file_press,
                                   file_covariance, file_correlation,
                                   file_balance, file_balance_plus, file_balance_minus);
            //}
            //else {
//                MC.output_dynamic(parameter_values,
//                                  file_biomass_species_mean, file_biomass_species_CV,
//                                  file_biomass_species_total_mean, file_biomass_species_total_CV,
//                                  file_biomass_community_mean, file_biomass_community_CV,
//                                  file_metacommunity, file_time_series);
            //}
            count_simu ++;
            cout << "simulation " << count_simu << "/" << n_simu << " ID=" << i+1 << " in thread " << omp_get_thread_num() << " - eigenvalue = " << GSL_REAL(eigen) << " + i" << GSL_IMAG(eigen) <<  endl;
        }
    }

    // Delete the table of parameters
    delete_table(params,n_simu);

    return 0;
}
