#include "community.h"
#include "species.h"
#include "functions.h"

#include "string"
#include "math.h"

using namespace std;




Species* Community::get_species_address(int i){
    return (_species[i]->get_address());
}
int Community::get_ID(){
    return (_ID);
}
// CREATION OF THE FOOD-WEB
void Community::set_object(){
    _species = new Species*[_nSpecies];
        for (int i=0; i<_nSpecies; i++){
            _species[i] = NULL;
        }
}
void Community::set_species(const double m, const double g, const double r, const double D, const double e, const double a, const double FR, const double h,
                            const double as_m, const double as_g, const double as_r, const double as_D, const double as_e, const double as_a, const double as_FR, const double as_h){
    _species[0] = new Species(as_m*m, as_g*g, as_D*D, as_e*e, as_a*a, 0, 0,
                              "Primary_producer", 0, 0 +_ID*_nSpecies); // create the primary producer at the bottom of the food chain
    for (int i=1; i<_nSpecies; i++){
        _species[i] = new Species(as_m*m, -r*as_r, as_D*D, as_e*e, as_a*a, as_FR*FR, as_h*h,
                                  "Consumer", i, i+_ID*_nSpecies); // creation the consumers
    }
}
void Community::set_interaction(){
    double feeding_rate(0);
    _species[0]->set_prey(); // _prey needs to be initialised for the primary producer
    for (int i=1; i<_nSpecies; i++){
        _species[i]->change_nPrey(1);
        _species[i-1]->change_nPred(1);
        _species[i]->set_prey();
        _species[i-1]->set_pred();
        feeding_rate = 0;//set_feeding_rate(a, _species[i-1]->get_metabolic(), _species[i]->get_metabolic(), v) // calculate the feeding rate
        _species[i]->add_prey(_species[i-1]->get_ID(), _species[i-1]->get_address(), feeding_rate); // add species i-1 in the diet of species i
        _species[i-1]->add_pred(_species[i]->get_ID(), _species[i]->get_address(), feeding_rate); // add species i-1 in the diet of species i
    }
    _species[_nSpecies-1]->set_pred(); // _pred needs to be initialised for the top predator
}
// constructor

Community::Community(int nSpecies,
                     int ID)
    :_nSpecies(nSpecies)
    ,_ID(ID){}
// destructor
Community::~Community(){
    for (int i=0; i<_nSpecies; i++){
        delete _species[i];
    }
    delete[] _species;
}
