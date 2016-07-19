
#include <iostream>
#include "Run.h"
main(){
	dsinit_();
	Cluster c;

c.name = "e1_im200x";
c.z = 1e-6;
c.rh = 1.4;
c.rcore = 0.1;

c.B0 = 1.0;


c.DM = 1;
c.rhos_Ein = 6.6;	
c.rs_Ein = 	0.15;
c.alpha = 0.3;
	

c.bsynch = 0.02549;
c.bIC = 0.2499;
c.bbrem = 0;
c.bcoul = 0;

c.DD = 0;
c.gamma = 0.7;
c.db = 1;
c.D0 =3e26;

c.n_r = 101; 
c.imNum = 250;
std::cout << c.greens(1*kpc2cm, 0.5*kpc2cm) << std::endl;

c.imNum = 25;
std::cout << c.greens(1*kpc2cm, 0.5*kpc2cm) << std::endl;
c.createGLUT();
runExCurve(25, 0.001);
c.SD = 0;
runExCurve(25, 0.001);



/*


///////////////////////////

c.name = "e2";
c.z = 1e-6;
c.rh = 1.4;
c.rcore = 0.1;

c.B0 = 0.6;


c.DM = 1;
c.rhos_Ein = 7.0;		
c.rs_Ein = 	0.25;	
c.alpha = 0.3;
	

c.bsynch = 0.02549;
c.bIC = 0.2499;
c.bbrem = 0;
c.bcoul = 0;

c.DD = 0;
c.gamma = 0.7;
c.db = 1;
c.D0 = 3e26;		


	runGreens_rdv(100);
	std::cout << "chill" <<std::endl;
	runGreens_r(10);
	std::cout << "chill" <<std::endl;
	runGreens_rdv(20);
	std::cout << "chill" <<std::endl;
	runGreens_r(20);
	std::cout << "chill" <<std::endl;
	runGreens_rdv(200);
	std::cout << "chill" <<std::endl;
	runGreens_r(15);
	std::cout << "chill" <<std::endl;

////////////////

c.name = "e3";
c.z = 1e-6;
c.rh = 1.4;
c.rcore = 0.1;

c.B0 = 1.0;

c.DM = 1;
c.rhos_Ein = 6.6;	
c.rs_Ein = 	0.15;
c.alpha = 0.3;


c.bsynch = 0.02549;
c.bIC = 0.2499;
c.bbrem = 0;
c.bcoul = 0;

c.DD = 0;
c.gamma = 0.7;
c.db = 1;
c.D0 = 3e27;		


	runGreens_rdv(100);
	std::cout << "chill" <<std::endl;
	runGreens_r(10);
	std::cout << "chill" <<std::endl;
	runGreens_rdv(20);
	std::cout << "chill" <<std::endl;
	runGreens_r(20);
	std::cout << "chill" <<std::endl;
	runGreens_rdv(200);
	std::cout << "chill" <<std::endl;
	runGreens_r(15);
	std::cout << "chill" <<std::endl;

	*/
}

/*
 g++ -o main Constants.cpp Cluster.cpp greens.cpp diffusion.cpp dist.cpp psyn.cpp pIC.cpp emissivity.cpp flux.cpp calc_sv.cpp run.cpp main.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran
*/

/*
why are runSED() and runFlux() so much slower
	slows even more at 48/50 ( runSED() )
write funcs, runGreens_rdv() runGreens_r()
 */