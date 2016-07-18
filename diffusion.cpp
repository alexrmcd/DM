
//***** diffusion.cpp *****//

#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include "Constants.h"
#include "Cluster.h"
#include "Darksusy.h"

double darksusy (double Ep){

	int yieldk = 151;
	int istat;

	double ds = dshayield_(&Cluster::p.mx, &Ep, &Cluster::p.ch, &yieldk, &istat);

	return ds;
}

double Cluster::ddiffusion(double Ep, void * params){
	std::vector<double> diffusionParams = *(std::vector<double> *)params;

	double E = diffusionParams[0];
	double vE = diffusionParams[1];
	double r = diffusionParams[2];


	double ddiffusion;


	if(SD == 1){
		double Ep_scaled = (int)(Ep/vscale) ;
		double rootdv = sqrt( std::abs(vE - vlookup[Ep_scaled]) ); 

		double r_scale = rh/n_r;
		double rootdv_scale = rootdv_max*kpc2cm/n_rootdv;

		double r_int = (int)(( 1 + r/kpc2cm )/r_scale);
		double rootdv_int = (int)(rootdv/rootdv_scale);

		if(r_int >= n_r)
			r_int = n_r - 1;

		ddiffusion = GLUT[r_int][rootdv_int] * darksusy(Ep);
		
		if(isnan(GLUT[r_int][rootdv_int]) == 1)
			std::cout << "GLUT " << GLUT[r_int][rootdv_int] << ", r = "<<r/kpc2cm << " r_int = "<< r_int <<", rootdv = " << rootdv/kpc2cm<<std::endl;
	}

	else{
		ddiffusion = darksusy(Ep);
	}

	return ddiffusion;

}

double Cluster::diffusion( double E, double r){			// int over Ep
	double E_scaled = (int)(E/vscale);

	double vE = vlookup[E_scaled];
	
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> diffusionParams (3);

	diffusionParams[0] = E;
	diffusionParams[1] = vE;
	diffusionParams[2] = r;

	gsl_function F;
	F.function = &ddiffusion;
	F.params = &diffusionParams; 								//pass Ep to rootdv(), pass r from dndeeq as well, 
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, E, p.mx, 0, 1e-2, 1000, 
	                    w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}

double Cluster::elec_spect(double E, double r ){

	double elec_spect = (1 / bloss(E,r))*diffusion(E, r);	

	if ( isnan(elec_spect) == 1)
		std::cout << "elec_spect(E = "<< E  << " , r = " << r/kpc2cm << " ) = "<< elec_spect << std::endl;
	
	return elec_spect;

}
