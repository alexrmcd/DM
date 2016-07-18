//***** calc_sv.cpp *****//

/*
g++ -c calc_sv.cpp
*/

#include "Run.h"

double min_flux(double r){

	double thetaB = 25.0; 					// beam size in arcsec
	double frms   = 1e-5; 					//noise per beam in Jy

	double dist_z = Cluster::dist() / (1.0 + Cluster::z);

	double thetaH = r/dist_z * 180.0/pi * 3600.0;

	double min_flux = 4.0 * log(2.0) * frms * pow(thetaH/thetaB, 2.0); 

	return min_flux;
}


double calc_sv(){ 
	double rcm = Cluster::rh * kpc2cm ; 
	double r = Cluster::rconst(rcm);	

	double flux_dm  = Cluster::ssyn()/Cluster::p.sv * GeVJy ; 
	double flux_obs = min_flux(r);

	double sv = flux_obs/flux_dm;

	return sv ; 
}

double calc_sv(double flux_obs ){

	double flux_dm  = Cluster::ssyn()/Cluster::p.sv * GeVJy ; 

	double sv =  flux_obs/flux_dm;

	return sv ; 
}