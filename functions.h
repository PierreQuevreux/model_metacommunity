#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "string"
using namespace std;
#include <gsl/gsl_matrix.h> // for integration

// GENERAL FUNCTIONS
void create_table(double** table, const int nrow, const int ncol); // create a two dimensions table
void delete_table(double** table, const int nrow); // delete a two dimension table

// PARALLELISATION
int load_parameters(const string file_path, int line_num);
string load_table(const string file_path, const int n_simu, double** params, const bool header); // load the file containing the parameters
string output_parameter(const string file_path, const int row, const bool header); // return a string version of parameters values
//string output_parameter(double** params, const int row, const int n_col); // return a string version of parameters values

// OTHER FUNCTIONS
double CV(const double mean, const double SQ, const int nPoint); // calculate the coefficient of variation
double SD(const double mean, const double SQ, const int nPoint); // calculate the standard deviation
double det_matrix(const gsl_matrix* matrix, const int n); // matrix determinant
gsl_matrix* invert_matrix(const gsl_matrix* matrix, const int n); // matrix inversion
void save_matrix(gsl_matrix* matrix, int n_row, int n_col, string file_path); // save the matrix as .txt file
gsl_vector* column_stack(const gsl_matrix* A, const int n_row, const int n_col); // turns the matrix A into a vector stacking A's columns
gsl_matrix* kronecker_product(const gsl_matrix* A, const gsl_matrix* B, const int n_row_A, const int n_col_A, const int n_row_B, const int n_col_B); // Kronecker product

//INTEGRATION
int multiroot_f (const gsl_vector* x, void* params, gsl_vector* f);
int multiroot_df (const gsl_vector* x, void* params, gsl_matrix* jacobian);
int multiroot_fdf (const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* jacobian);
// System of ODE, function friend of all classes
int func(double t, const double y[], double f[], void *params); // system of ODE

#endif // FUNCTIONS_H_INCLUDED
