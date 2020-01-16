#include "species.h"
#include "functions.h"

#include <iostream>
#include "string"
#include "math.h"
#include <cstdlib> // to use calloc
using namespace std;

// GENERAL FUNCTIONS
Species* Species::get_address(){
    return (this);
}
double Species::get_m(){
    return(_m);
}
double Species::get_g(){
    return(_g);
}
double Species::get_D(){
    return(_D);
}
double Species::get_e(){
    return(_e);
}
double Species::get_a(){
    return(_a);
}
double Species::get_FR(){
    return(_FR);
}
double Species::get_h(){
    return(_h);
}
string Species::get_type(){
    return(_type);
}
int Species::get_IDc(){
    return(_IDc);
}
int Species::get_ID(){
    return(_ID);
}

// FUNCTIONS RELATIVE TO PREDATORS
int Species::get_nPred(){
    return(_nPred);
}
void Species::change_nPred(int nPred){
    _nPred = nPred;
}
void Species::set_pred(){
    _pred = new Species*[_nPred];
        for (int i=0; i<_nPred; i++){
            _pred[i] = NULL;
        }
    _predID = new int[_nPred]; // create the vector of predators ID (or rank)
        for (int i=0; i<_nPred; i++){
            _predID[i] = 0;
        }
    _aPred = new double[_nPred]; // create the vector of interaction strength with predators
        for (int i=0; i<_nPred; i++){
            _aPred[i] = 0;
        }
}
void Species::add_pred(int predID, Species* pred, double a){
    _predID[_countPred] = predID; // position of the predator
    _pred[_countPred] = pred; // pointer toward the predator
    _aPred[_countPred] = a; // interaction strength with predators
    _countPred++; // counter of predators
}
int Species::get_predID(int i){
    return (_predID[i]);
}

// FUNCTIONS RELATIVE TO PREY
int Species::get_nPrey(){
    return(_nPrey);
}
void Species::change_nPrey(int nPrey){
    _nPrey = nPrey;
}
void Species::set_prey(){
    _prey = new Species*[_nPrey];
        for (int i=0; i<_nPrey; i++){
            _prey[i] = NULL;
        }
    _preyID = new int[_nPrey]; // create the vector of predators ID (or rank)
        for (int i=0; i<_nPrey; i++){
            _preyID[i] = 0;
        }
    _aPrey = new double[_nPrey]; // create the vector of interaction strength with predators
        for (int i=0; i<_nPrey; i++){
            _aPrey[i] = 0;
        }
}
void Species::add_prey(int preyID, Species* prey, double a){
    _preyID[_countPrey] = preyID; // position of the prey in the vector
    _prey[_countPrey] = prey; // pointer toward the prey
    _aPrey[_countPrey] = a; // interaction strength with prey
}
int Species::get_preyID(int i){
    return (_preyID[i]);
}

// constructor
Species::Species(double m,
                 double g,
                 double D,
                 double e,
                 double a,
                 double FR,
                 double h,
                 string type,
                 int IDc,
                 int ID)
                 :_m(m)
                 ,_g(g)
                 ,_D(D)
                 ,_e(e)
                 ,_a(a)
                 ,_FR(FR)
                 ,_h(h)
                 ,_type(type)
                 ,_IDc(IDc)
                 ,_ID(ID)
                 ,_nPrey(0)
                 ,_countPrey(0)
                 ,_nPred(0)
                 ,_countPred(0){}
// destructor
Species::~Species(){
    delete[] _pred;
    delete[] _predID;
    delete[] _aPred;
    delete[] _prey;
    delete[] _preyID;
    delete[] _aPrey;
}
