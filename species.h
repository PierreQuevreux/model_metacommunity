#ifndef SPECIES_H_INCLUDED
#define SPECIES_H_INCLUDED

#include <iostream>
#include "string"
#include <cstdlib> // to use calloc
using namespace std;

class Species
{
protected :
    // GENERAL PARAMETERS
    double _m; // metabolic rate
    double _g; // density independent demographic process
    double _D; // self regulation
    double _e; // assimilation efficiency
    double _a; // attack rate
    double _FR; // Hill coefficient
    double _h; // handling time
    string _type; // type of the species
    int _IDc; // rank of the species within a community
    int _ID; // rank of the species within a metacommunity

    int _diversity; // number of species
    int _nResource; // number of resources
    int _nba; // number of possible interactions

    // PREYS
    int _nPrey; // number of prey
    int* _preyID; // array containing the ID (or rank in the array) of prey
    Species** _prey; // array containing the pointers toward prey
    double* _aPrey; // interaction strength with prey
    int _countPrey; // counter to fill in the array of prey pointers

    // PREDATORS
    int _nPred; // number of predators
    int* _predID; // array containing the number of the predators
    Species** _pred; // array containing the pointers toward predators
    double* _aPred; // interaction strength with predators
    int _countPred; // counter to fill in the array of predators pointers

public :

    // GENERAL FUNCTIONS
    Species* get_address(); // return the pointer this
    double get_mass(); // function returning the body mass
    double get_m(); // function returning the metabolic rate
    double get_g(); // function returning the growth/mortality rate
    double get_D(); // function returning the self regulation coefficient
    double get_e(); // function returning the assimilation efficiency
    double get_a(); // function returning the attack rate
    double get_FR(); // function returning the Hill coefficient
    double get_h(); // function returning the handling time
    string get_type(); // function returning the type
    int get_IDc(); // function returning the rank within the community
    int get_ID(); // function returning the rank within the metacommunity

    // FUNCTIONS RELATIVE TO PREDATORS
    int get_nPred(); // function returning the number of predators
    void change_nPred(int nPred); // change the total number of predators
    void set_pred(); // create the vector of predators'pointers
    void add_pred(int nPred, Species* pred, double a); // add the pointer of the predator to the predators list
    int get_predID(int i); // return the index of the ith pred

    // FUNCTIONS RELATIVE TO PREYS
    int get_nPrey(); // return the number of preys
    void change_nPrey(int nPrey); // change the number of prey
    void set_prey(); // create the array containing the pointers of preys
    void add_prey(int nPrey, Species* prey, double a); // add the pointer of a prey to the diet
    int get_preyID(int i); // return the index of the ith prey

    //constructor
    Species(double m, double g, double D, double e, double a, double FR, double h, string type, int IDc, int ID);
    // destructor
    ~Species();
};

#endif // SPECIES_H_INCLUDED
