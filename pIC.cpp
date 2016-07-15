
//***** pIC.cpp *****//

/*
g++ -c pIC.cpp 
*/

#include <iostream>

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include "Constants.h"
#include "Cluster.h"

double Cluster::G(double eps, double E_gamma, double boost){

	double GAMMA = 4 * eps * boost/me;
	double q = E_gamma / (GAMMA * (boost * me  - E_gamma) );

	double G = 2 * q* log(q) + (1+2*q)*( 1 - q ) + pow(GAMMA*q , 2)*(1-q) / ( 2*(1+GAMMA*q ) );//
	return G;
}

double Cluster::IC_cross_Section( double eps, double E_gamma, double E){

	double boost = E/me;

	double sigma = 3 * sigma_thompson/(4* eps * pow(boost ,2))* G(eps, E_gamma, boost);

	return sigma;
}

double Cluster::CMB_bbSpectrum(double eps){

	double T = 2.73;

	double nu = eps/(hplanck*J2Gev);

	double CMB_bbSpectrum = 8*pi* pow(nu, 2.0)/pow(clight, 3.0)/(	exp(hplanck * nu /(kb * T )) - 1	);

	return CMB_bbSpectrum;
}

double Cluster::dpIC(double eps, void * params ){

	std::vector<double> pICParams = *(std::vector<double> *)params;
	double E_gamma = pICParams[0] ;
	double E = pICParams[1] ;
	
	double q = pow(me,2)*E_gamma / ( 4 * eps * E * (E  - E_gamma) );
	
	double dpIC;
	if( q>1 or q < 1/(4* pow(E/me , 2) ) or eps/hplanck/J2Gev > 1e13){
		dpIC = 0;
	}
	else{
		dpIC = CMB_bbSpectrum(eps) *  IC_cross_Section(eps, E_gamma, E); 
	};
	return dpIC;
}


double Cluster::pIC(double nu, double E){			//int over eps

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (5000);

	double result, error, pIC;
	size_t size;

	std::vector<double> pICParams (2);

	double E_gamma = hplanck * nu * J2Gev;
	double boost = E/me;
	double eps_max = E_gamma*E/(E - E_gamma );//1.24e-12;//
	double eps_min = E_gamma/(4 * boost/me*( E - E_gamma)); //1.24e-15;//
		//std::cout << eps_max <<", E_gamma = " << E_gamma <<std::endl;
	if (eps_max > 1e13*hplanck*J2Gev){
		eps_max = 1e13*hplanck*J2Gev;

	}
	pICParams[0] = E_gamma;
	pICParams[1] = E;

	gsl_function F;
	F.function = &dpIC;
	F.params = &pICParams;
	gsl_set_error_handler_off();
	gsl_integration_qng (&F,  eps_min, eps_max,  0, 1e-3, &result, &error, &size); 

	gsl_integration_workspace_free (w);

	if(E<E_gamma)
		pIC = 0;
	else
		pIC = clight* E_gamma*result; 

	return pIC;
}
