#ifndef COMMUNITY_H_INCLUDED
#define COMMUNITY_H_INCLUDED

#include "species.h"
#include "string"
using namespace std;

class Community
{
private :
    //PARAMETERS
    int _nSpecies; // number of species
    int _ID; // number of the community

    //STRUCTURES
    Species** _species; // array of species

public :
    Species* get_species_address(int i); // return the address of species i
    int get_ID(); // return the ID of the community
    // CREATION OF THE FOODWEB
    void set_object(); // create the vectors of populations and parameters
    void set_species(const double m, const double g, const double r, const double D, const double e, const double a, const double FR, const double h,
                     const double as_m, const double as_g, const double as_r, const double as_D, const double as_e, const double as_a, const double as_FR, const double as_h); // create the primary producers
    void set_interaction(); // defines the interactions between species

    //constructor
    Community(int diversity, int ID);
    // destructor
    ~Community();
};

#endif // COMMUNITY_H_INCLUDED
