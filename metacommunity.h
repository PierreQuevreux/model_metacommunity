#ifndef METACOMMUNITY_H
#define METACOMMUNITY_H

#include "community.h"
#include "species.h"

// INTEGRATION
#include <stdio.h>
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration, matrix and vectors
#include <gsl/gsl_odeiv2.h> // for integration
#include <gsl/gsl_complex.h> // for complex numbers

#include <iostream>
#include <fstream>
#include <sstream>

#include <gsl/gsl_rng.h> // random number

class Metacommunity
{
private:
    // PARAMETERS
    int _nSpecies; // total number of different species
    int _nCommunity; // number of communities (=patches)
    int _dim; // total dimension of the system
    int _dim_prtb; // size of the perturbation matrix (exogenous, demographic and environmental)
    int _n_step_record; // maximum number of time steps
    //STRUCTURES
    Species** _species; // array of species
    Community** _community; // vector with the pointers of each community
    double** _dispersal; // dispersal matrix between patches
    double** _asymmetry; // asymmetry matrix between patches
    gsl_matrix* _perturbation; // effect of white noise on species dynamics
    gsl_matrix* _response_press; // invert of the Jacobian matrix
    gsl_vector* _covariance; // vector with the stacked columns of the variance-covariance matrix
    gsl_vector* _correlation; // vector with the stacked columns of the correlation matrix

    gsl_rng* _rng; // random number generator
    int _nPoint; // number of recorded points counter for the calculation of CVs
    int _nPoint_TS; // number of recorded points for the time series table

    // RECORDED OUT PUT
    double* _t_extinction; // time at which species get extinct
    // Biomass of each species
    double* _biomass_species; // value at t
    double* _biomass_species_equilibrium; // species biomass equilibrium computed analytically
    double* _biomass_species_mean; // sum of values
    double* _biomass_species_SQ; // sum of square
    double* _biomass_species_CV; // coefficient of variation
    double* _balance; // ratio of dispersal process to demographic process
    double* _balance_plus; // ratio of dispersal process to demographic process (positive terms)
    double* _balance_minus; // ratio of dispersal process to demographic process (negative terms)
    // Biomass of each species across communities
    double* _biomass_species_total; // value at t
    double* _biomass_species_total_mean; // sum of values
    double* _biomass_species_total_SQ; // sum of square
    double* _biomass_species_total_CV; // coefficient of variation
    // Total biomass of a community
    double* _biomass_community; // value at t
    double* _biomass_community_mean; // sum of values
    double* _biomass_community_SQ; // sum of square
    double* _biomass_community_CV; // coefficient of variation
    // Total biomass
    double _biomass_total; // value at t
    double _biomass_total_mean; // sum of values
    double _biomass_total_SQ; // sum of square
    double _biomass_total_CV; // coefficient of variation
    // Time series
    double** _time_series;

public:
    int get_dim(); // return the dimension of the system
    void set_object(int seed); // create vectors
    void set_dispersal_matrix(const string file_matrix); // create and load the dispersal matrix
    void set_perturbation_matrix(const string file_matrix); // load the covariance matrix of the perturbation
    void set_asymmetry_matrix(const string file_matrix, const int n_params_asymmetry); // load the asymmetry matrix between patches
    void set_community(const double m, const double g, const double r, const double D, const double e, const double a, const double FR, const double h); // create the communities
    // INTEGRATION
    void Derivative(const double* y, double* f); // compute the derivative of biomass dynamics
    void Balance(const double* y); // compute the ratio of dispersal process to demographic process
    void Jacobian(const gsl_vector* x, gsl_matrix* jacobian); // compute the Jacobian matrix
    void Interaction(gsl_matrix* jacobian); // compute the interaction matrix
    gsl_complex Equilibrium(); // compute the equilibrium without dispersal and return the value of the largest eigenvalue
    void Perturbation_equilibrium(); // compute the response to a press and a stochastic perturbation
    void Dynamic_Euler_Maruyama(double t, const double t_step, const int n_step_trans, const int n_step_record, double t_step_record, const double extinction_threshold); // SDE solver
    // RECORDS THE OUTPUTS
    void output_analytic(const string parameter_values, ofstream& file_biomass_species_mean, ofstream& file_press,
                         ofstream& file_covariance, ofstream& file_correlation,
                         ofstream& file_balance, ofstream& file_balance_plus, ofstream& file_balance_minus);
    void output_dynamic(const string parameter_values,
                        ofstream& file_biomass_species_mean, ofstream& file_biomass_species_CV,
                        ofstream& file_biomass_species_total_mean, ofstream& file_biomass_species_total_CV,
                        ofstream& file_biomass_community_mean, ofstream& file_biomass_community_CV,
                        ofstream& file_metacommunity, ofstream& file_time_series);

    Metacommunity(int nSpecies,
                  int nCommunity,
                  int dim,
                  int n_step);
    virtual ~Metacommunity();
};

#endif // METACOMMUNITY_H
