
//***** psyn.cpp *****//

/*
g++ -c psyn.cpp 
*/

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Functions.h"
//Synchrotron emmission spectral function form Cola2006.
double fff(double x){

	double fff = 1.25 * pow( x , 1.0/3.0) * exp( -x )* pow((648 + x*x) , 1.0/12.0);

	return fff;
}

double dpsyn(double theta, void * params ){

	std::vector<double> psynParams = *(std::vector<double> *)params;
	double E = psynParams[0];
	double r = psynParams[1];

	double nu = 1.4;					//Ghz
	double psyn0 = 1.46323e-25 ; 		//Gev/s/Hz
	double x0 = 62.1881 ;				//dimensionless constant
	double nu_em = ( 1 + c.z )* nu; 	//(observing freq)*(1+z)


	double x = x0 *nu_em / ( c.bfield_model( r ) * pow( E, 2)	 );
	//std::cout << "x = " << x  << " B = " << c.bfield_model(psynParams[1]) << std::endl;
	double dpsyn = psyn0 * c.bfield_model(r) * 0.5 * pow(  sin(theta) , 2)* fff( x  /sin(theta) ); 
	
	

	return dpsyn;

}

double psyn(  double E, double r){			//int over theta

		
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> psynParams (2);

	psynParams[0] = E;
	psynParams[1] = r;


	gsl_function F;
	F.function = &dpsyn;
	F.params = &psynParams;

	gsl_integration_qags (&F, 1e-16, pi, 0, 1e-3, 1000, //x?
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}
