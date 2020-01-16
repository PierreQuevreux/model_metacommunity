# model_metacommunity

Model used in the article "Synchrony and perturbation transmission in trophic metacommunities"

Program coded in C++ and paralellised with openMP.

# generator_matrix.R

R code producing the following matrices :
- the perturbation matrices saved as matrix_perturbation_N.txt containing the the variance of each independent perturbation
- the dispersal matricies saved as matrix_dispersal_N.txt containg the dispersal rates
- the asymmetry matrices saved as matrix_asymmetry_N.txt containing the asymmetry of parameters between patches
All the informations on the matrices (ID, parameters...) are in the parameters_matrix.txt file.

# models.txt

File selecting the matrices that will be used during the simulations.

# generator_parameter_table.R

R code generating the file containing all the parameters and the ID of matrices used to feed the simulations. The parameters and the informatiopns on parameters are contained in the files parameters_N.txt and parameters_data_N.txt.

# aggregation_data.R

R code merging the results of the different simulations. The output data are saved in the "results" folder containing a README.txt detailing the recorded variables and parameters.

# Figures.R

R code to produce the figures of the article
