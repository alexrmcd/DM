
#include <iostream>
#include "Run.h"
main(){
	dsinit_();
//	Cluster::createGLUT();
	Cluster::SD = 0;
	runExCurve(25, 0.64);
	runSED(40);
	runFlux(40);
}

/*
 g++ -o main Constants.cpp Cluster.cpp greens.cpp diffusion.cpp dist.cpp psyn.cpp pIC.cpp emissivity.cpp flux.cpp calc_sv.cpp run.cpp main.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
*/

/*
 need to do diffusion, so incorporate darksusy, also get elec_spect. issue in emissivity with misusing p.mx?? 
 don't get it, hopefully will fix itself later. can restructure a bit tomorrow and keep more or less same organization as in dark_matter, 
 just with everything through ssyn, sic as Cluster member functions. 
 */