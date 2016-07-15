
#include <iostream>
#include "Cluster.h"
#include "Darksusy.h"
#include "Functions.h"
main(){
	Cluster fish;
	std::cout << fish.rconst(1.4 *kpc2cm )/kpc2cm<< " " << fish.rconst(1500. *kpc2cm )/mpc2cm<<std::endl;
}

/*
 g++ -o main Constants.cpp Cluster.cpp dist.cpp pIC.cpp emissivity.cpp main.cpp -lgsl -lgslcblas

*/

/*
 need to do diffusion, so incorporate darksusy, also get elec_spect. issue in emissivity with misusing p.mx?? 
 don't get it, hopefully will fix itself later. can restructure a bit tomorrow and keep more or less same organization as in dark_matter, 
 just with everything through ssyn, sic as Cluster member functions. 
 */