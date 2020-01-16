#include "functions.h"
#include "metacommunity.h"
#include "community.h"
#include "species.h"

#include <iostream>
#include <fstream>
#include <sstream> // to use istringstream and convert string to double

#include <string>
using namespace std;

#include <time.h>
#include "math.h"
#include <gsl/gsl_errno.h> // for integration
#include <gsl/gsl_matrix.h> // for integration
#include <gsl/gsl_odeiv.h> // for integration
#include <gsl/gsl_randist.h> // for random distributions
#include <gsl/gsl_linalg.h> // linear algebra with matrix and vectors

// GENERAL FUNCTION
void create_table(double** table, const int nrow, const int ncol){
    for(int i=0; i<nrow; i++){
        table[i] = new double[ncol];
        for(int j=0; j<ncol; j++){
            table[i][j] = 0;
        }
    }
}
void delete_table(double** table, const int nrow){
    for (int i=0; i<nrow; i++){
        delete[] table[i];
        table[i] = NULL;
    }
    delete[] table;
}

// PARALLELISATION
int load_parameters(const string file_path, int line_num){
    int out;
    // OPENING OF THE FILE
    ifstream file(file_path.c_str()); // identify the file to open
    if (file){} // control the opening of the file
    else {cout << "ERROR: unable to open file at path" << endl;}
    // DECLARATION OF THE VARIABLEs
    string line; // line of the file
    istringstream iss; // variable for the conversion from string to double

    for(int i=1; i<line_num+1; i++){
        getline(file, line); // read the line of 'file' and store it in 'line'
    }
    iss.str( line.c_str() ); // convert 'value' from string to double
    iss >> out; // write ,the value of the parameter in the table
    return (out);
}
string load_table(const string file_path, const int nRow, double** table, const bool header){
    // OPENING OF THE FILE
    ifstream file(file_path.c_str()); // identify the file to open
    if (file){} // control the opening of the file
    else {cout << "ERROR: unable to open " << file_path << endl;}
    // DECLARATION OF THE VARIABLEs
    string header_names; // names of the variables in the header
    string line; // line of the file
    string letter; // character in the line
    string value; // value of the parameter (concatenation of the values of 'letter')
    istringstream iss; // variable for the conversion from string to double
    int k(0); // counter of columns
    // HEADER
    if (header){
        getline(file, header_names); // skip the header
    }
    //

    for(int i=0; i<nRow; i++){
        getline(file, line); // read the line of 'file' and store it in 'line'
        value.clear(); // clear value
        k = 0;
        for (unsigned int j=0; j<line.size(); j++){
            letter = line[j];
            if (letter!=";"){
                value += line[j]; // write the data letter by letter
            }
            else {
                iss.str( value.c_str() ); // convert 'value' from string to double
                iss >> table[i][k];  // write ,the value of the parameter in the table
                iss.clear(); // clear iss
                value.clear(); // clear value
                k++; // next variable
            }
        }
        iss.str( value.c_str() ); // convert 'value' from string to double
        iss >> table[i][k]; // write ,the value of the parameter in the table
        iss.clear(); // clear iss
        value.clear(); // clear value
    }
    return (header_names); // return the header
}
string output_parameter(const string file_path, const int row, const bool header){
    ifstream file(file_path.c_str()); // identify the file to open
    if (file){} // control the opening of the file
    else {cout << "ERROR: unable to open " << file_path << endl;}
    // DECLARATION OF THE VARIABLEs
    string line; // line of the file
    // HEADER
    if (header){
        getline(file, line); // skip the header
    }
    //
    for(int i=0; i<row; i++){
        getline(file, line); // read the line of 'file' and store it in 'line'
    }
    return(line);
}
//string output_parameter(double** params, const int row, const int n_col){
//    string output_parameter;
//    string value; // string version of the int
//    for (int i=0; i<n_col-1; i++){
//        value = to_string(params[row][i]);
//        //cout << params[row][i] << ";";
//        output_parameter += value + ";";
//    }
//    //cout << endl;
//    value = to_string(params[row][n_col-1]);
//        output_parameter += value;
//        //cout << output_parameter << endl;
//    return(output_parameter);
//}

// OTHER FUNCTIONS
double CV(const double mean, const double SQ, const int nPoint){
    return(pow(SQ/nPoint - pow(mean,2),0.5)/mean);
}
double SD(const double mean, const double SQ, const int nPoint){
    return(pow(SQ/nPoint - pow(mean,2),0.5));
}
double det_matrix(const gsl_matrix* matrix, const int n){
    gsl_permutation* p = gsl_permutation_alloc(n);
    gsl_matrix* LU = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(LU, matrix);
    int s; // sign of the permutation
    gsl_linalg_LU_decomp(LU, p, &s); // compute the LU decomposition of this matrix
    double det = gsl_linalg_LU_det(LU, s); // invert the matrix
    gsl_permutation_free(p);
    gsl_matrix_free(LU);
    return det;
}
gsl_matrix* invert_matrix(const gsl_matrix* matrix, const int n){
    gsl_permutation* p = gsl_permutation_alloc(n);
    gsl_matrix* LU = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(LU, matrix);
    int s; // sign of the permutation
    gsl_linalg_LU_decomp(LU, p, &s); // compute the LU decomposition of this matrix
    gsl_matrix* inv = gsl_matrix_alloc(n, n); // compute the  inverse of the LU decomposition
    gsl_linalg_LU_invert(LU, p, inv); // invert the matrix
    gsl_permutation_free(p);
    gsl_matrix_free(LU);
    return inv;
}
void save_matrix(gsl_matrix* matrix, int n_row, int n_col, string file_path){
    ofstream file_matrix (file_path.c_str()); // creation of the file
    for (int i=0; i<n_row; i++){
        for (int j=0; j<n_col-1; j++){
            file_matrix << gsl_matrix_get(matrix, i, j) << ";"; // each species of each patch
        }
    file_matrix << gsl_matrix_get(matrix, i, n_col-1) << endl;
    }
}
gsl_vector* column_stack(const gsl_matrix* A, const int n_row, const int n_col){
    gsl_vector* V = gsl_vector_alloc(n_row*n_col);
    for (int j=0; j<n_col; j++){
        for (int i=0; i<n_row; i++){
            gsl_vector_set(V, j*n_row+i, gsl_matrix_get(A, i, j));
        }
    }
    return V;
}
gsl_matrix* kronecker_product(const gsl_matrix* A, const gsl_matrix* B, const int n_row_A, const int n_col_A, const int n_row_B, const int n_col_B){
    gsl_matrix *kro = gsl_matrix_alloc(n_row_A*n_row_B, n_col_A*n_col_B); // matrix for the Kronecker product
    for (int row_A=0; row_A<n_row_A; row_A++){
        for (int col_A=0; col_A<n_col_A; col_A++){
            for (int row_B=0; row_B<n_row_B; row_B++){
                for (int col_B=0; col_B<n_col_B; col_B++){
                    gsl_matrix_set(kro, row_A*n_row_B+row_B, col_A*n_col_B+col_B, gsl_matrix_get(A, row_A, col_A)*gsl_matrix_get(B, row_B, col_B));
                }
            }
        }
    }
    return kro;
}

//INTEGRATION
// Root finder
int multiroot_f (const gsl_vector* x, void* params, gsl_vector* f){
    Metacommunity* MC = (Metacommunity*)params;
    double* array_x =  new double[MC->get_dim()];
    double* array_f =  new double[MC->get_dim()];
    for (int i=0; i<MC->get_dim();i++){
        array_x[i] = gsl_vector_get(x,i);
        array_f[i] = gsl_vector_get(f,i);
    }
    MC->Derivative(array_x, array_f); // direct access to the array included in vectors
    for (int i=0; i<MC->get_dim();i++){
        gsl_vector_set(f,i,array_f[i]);
    }
    delete[] array_x;
    delete[] array_f;
    return GSL_SUCCESS;
}
int multiroot_df (const gsl_vector* x, void* params, gsl_matrix* jacobian){
    Metacommunity* MC = (Metacommunity*)params;
    MC->Jacobian(x, jacobian); // compute the Jacobian matrix
   return GSL_SUCCESS;
}
int multiroot_fdf (const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* jacobian){
    multiroot_f (x, params, f);
    multiroot_df (x, params, jacobian);
   return GSL_SUCCESS;
}
