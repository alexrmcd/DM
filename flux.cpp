
//***** flux.cpp ******//

#include <iostream>

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Cluster.h"

double Cluster::dsIC( double r, void * params ){

	double nu = *(double *)params;

	double dist_z = dist() / (1+z);

	double sICIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(DM_profile(r) , 2)*jIC(nu, r);

	return sICIntegrand;
}



double Cluster::sIC( double nu, double r ){				

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &dsIC;
	F.params = &nu;

	gsl_integration_qags (&F, 1e-16, r, 0, 1e-3, 1000,
	                w, &result, &error); 

	result *= p.sv/(8* pi*pow( p.mx , 2.0 ));

	return result;

}


double Cluster::dssyn( double r, void * params ){

	double dist_z = dist() / (1+z);

	double ssynIntegrand = 4 *pi /pow(dist_z , 2) *pow(r,2)  *  pow(DM_profile(r) , 2)*jsyn(r);	
	
	return ssynIntegrand;
}


double Cluster::ssyn( double r ){				// int over r

	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &dssyn;

	gsl_integration_qags (&F, 1e-16 ,  r, 0, 1e-3, 1000,
	                w, &result, &error); 
	
	result *= p.sv/(8* pi*pow( p.mx , 2.0 ));

	return result;

}








