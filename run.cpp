//*****run.cpp*****//

#include <ctime>
#include <fstream>
#include <sstream>

#include "Run.h"

std::string channel(){
	std::string channel;

	if(Cluster::p.ch == 13){
		channel = "WW";
	}
	else if(Cluster::p.ch == 15){
		channel = "ee";
	}
	else if(Cluster::p.ch == 17){
		channel = "mumu";
	}
	else if(Cluster::p.ch == 19){
		channel = "tt";
	}
	else if(Cluster::p.ch == 25){
		channel = "bb";
	};

	return channel;
}

void runExCurve(int ch, double flux_obs){
	

	Cluster::p.ch = ch;							
	//Cluster::createGLUT();

	double mx_min; 
	
	if(Cluster::p.ch == 13)
		mx_min = 81;
	else
		mx_min = 5;
	
	
	double mx_max = 1000;
	
	double data;

	std::ostringstream makefilename;
	makefilename <<"NULL"<<Cluster::name << "_" << channel() << "_gamma_"<<Cluster::gamma <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());


	int n_mx = 50 ;//number of mx values used

	for (int i = 0 ; i < n_mx + 1 ; ++i){

		// iteration timer start
		std::clock_t start;
		double duration;
		start = std::clock();
		///////before algorithm

			Cluster::p.mx = mx_min * ( exp(    (log(mx_max) - log(mx_min))/ n_mx * i));

			if(Cluster::SD==1)
				Cluster::create_vLUT();

			if (flux_obs == 0)
				data = calc_sv();
			else
				data = calc_sv(flux_obs);

			file << Cluster::p.mx << "\t" <<  data <<std::endl;
			std::cout << "sv( " << Cluster::p.mx << " ) = " << data << std::endl;

		////////after algorithm
		duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
		std::cout << i << "/"<< n_mx << " ";
		std::cout << Cluster::p.ch << ", " << Cluster::name << ", time = " << duration <<std::endl;

	};

}

void runSED(double mx){ 

	Cluster::p.mx = mx;
	int n_nu = 50;
	if(Cluster::SD == 1) 
		Cluster::create_vLUT();

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	///////before algorithm

	std::ostringstream makefilename;
	makefilename << "IC_SED_JUNK" << Cluster::p.mx << "Gev_coma.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n_nu + 1; ++i  ){
		//total time timer start
		std::clock_t iter_start;
		double iter_duration;
		iter_start = std::clock();
		///////before algorithm
		double nu_min = 1e12;
		double nu_max = 1e23;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n_nu * i));
		double data = Gev2erg*Cluster::sIC( nu )*nu ;
		file << log10(nu) << "\t" << data <<std::endl;

		iter_duration = (std::clock()  -  iter_start)/(double) CLOCKS_PER_SEC;
		std::cout << i<< "/" << n_nu << " " << log10(nu) << "\t" << data << ", duration: " << iter_duration << std::endl;

	};
	
	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}


void runFlux(double mx){ 
	
	Cluster::p.mx = mx;
	int n = 50;
	if(Cluster::SD == 1) 
		Cluster::create_vLUT();

	//total time timer start
	std::clock_t start;
	double duration;
	start = std::clock();
	///////before algorithm

	std::ostringstream makefilename;
	makefilename <<Cluster::name << "_SynchFlux_" << Cluster::p.mx << "Gev_coma.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int i = 0 ; i < n + 1; ++i  ){
		std::clock_t Sstart;
		double Stime;
		Sstart = std::clock();
		///////before algorithm

		double nu_min = 1e6;
		double nu_max = 1e11;
	
		double nu = nu_min * ( exp(    (log(nu_max) - log(nu_min))/ n * i));
		double data = GeVJy*Cluster::ssyn(nu);
		file << nu << "\t" <<  data <<std::endl;

		////////after algorithm
		Stime = (std::clock()  -  Sstart)/(double) CLOCKS_PER_SEC;
		std::cout << "\n" <<Cluster::name << ", "<< Cluster::p.mx <<" : "<< nu << "\t" << data <<" duration = "<< Stime << std::endl;

	};
	
	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	std::cout << "Total time:  " << duration <<std::endl;
}

