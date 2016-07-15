//***** emissivity.cpp *****//

#include <ctime>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Cluster.h"

double Cluster::djIC(double E , void * params){

	std::vector<double> jICParams = *(std::vector<double> *)params;
	double nu = jICParams[0];
	double r = jICParams[1];

	double djIC = 2* pIC(nu, E)* elec_spect(E , r);

	return djIC;
}


double Cluster::jIC(double nu, double r){ 				// int over E

	///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	///////////


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;
	std::vector<double> jICParams (2);

	jICParams[0] = nu;
	jICParams[1] = r;

	gsl_function F;
	F.function = &djIC;
	F.params = &jICParams;

	gsl_integration_qags (&F, me , p.mx, 0, 1e-2, 1000,
	                    w, &result, &error); 


	gsl_integration_workspace_free (w);
	
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	if (duration > 1)
		std::cout << "jsyn( r = " << r/kpc2cm <<" ) = "<< result <<" duration: " << duration <<std::endl;  //      ~30s-60s

	return result;
}

double Cluster::djsyn(double E , void * params){

	double r = *(double *)params;

	double djsyn = 2* psyn(E, r)* elec_spect(E , r);

	return djsyn;
}


double Cluster::jsyn(double r){ 				// int over E

	///////////
	std::clock_t start;
	double duration;
	start = std::clock();
	///////////

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &djsyn;
	F.params = &r;

	gsl_integration_qags (&F, me, p.mx, 0, 1e-3, 1000,
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;
	if (duration > 1)
		std::cout << "jsyn( r = " << r/kpc2cm <<" ) = "<< result <<" duration: " << duration <<std::endl;  //      ~30s-60s

	return result;
}
