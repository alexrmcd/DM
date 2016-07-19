
#include <iostream>
#include "Run.h"
main(){
	dsinit_();
	Cluster::createGLUT();
	//Cluster::SD = 0;
	//runExCurve(25, 0.64);
	//runSED(40);
	runFlux(40);
}

/*
 g++ -o main Constants.cpp Cluster.cpp greens.cpp diffusion.cpp dist.cpp psyn.cpp pIC.cpp emissivity.cpp flux.cpp calc_sv.cpp run.cpp main.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
*/

/*
why are runSED() and runFlux() so much slower
	slows even more at 48/50 ( runSED() )
need to make bmu ave of bfield
write funcs, runGreens_rdv() runGreens_r()
 */