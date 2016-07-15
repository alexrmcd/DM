//***** Cluster.h *****//

#include <ctime>
#include <iostream>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Particle.h"

class Cluster{  

	public:

	std::string name;

	//Cluster Parameters
	static double z ;					//redshift
	static double rh; 					//halo radius in kpc
	static double rcore ; 				//core radius in kpc

	//Bfield parameters
	static double B0 ; 				//microGauss
	static double Bmu;
	static double beta;
	static double eta;

	//DM density parameters
	static int DM;						//selects DM profile, 0 -> NFW, 1-> Einasto
	static double rhos_NFW;			//characteristic density for NFW in GeV/cm^3
	static double rs_NFW;				//scale radius for NFW in kpc
	static double rs_NFWa;
	static double rhos_Ein;			//characteristic density for Einasto in GeV/cm^3
	static double rs_Ein;				//scale radius for Einasto in kpc
	static double alpha;

	//energy loss coeff. in 1e-16 Gev/s 
	static double bsynch;				
	static double bIC;
	static double bbrem;
	static double bcoul;
	
	// Diffusion Parameters
	int SD;						//switch to turn on diffusion, 0 -> NSD, 1-> SD
	static double gamma;				// D(E) ~ E^gamma
	static double db;					// minimum scale of uniformity for mag field, from Colafrancesco. 
	static double D0;					//in cm^2/s

	// v lookup table parameters
	int vsize;
	double vscale;
	std::vector<double> vlookup;

	// Greens lookupt table parameters 
	int imNum;
	int n_r;
	int n_rootdv;
	int rootdv_max;
	std::vector< std::vector<double> > GLUT;

	Particle p;


	Cluster() ;
	static double ddist(double z  , void * params);
	static double dist();
	static double rconst(double rcm);

	double bfield_model (double r) ;

	static double DM_profile(double r);

	double bloss(double E , double r);
	static double bloss(double E );

	static double D(double E);
	static double dv(double E , void * params);
	double v( double E );
	double root_dv(double Ep,  double vE);
	static double dgreens_term(double rp, void * params );
	double greens_term(double ri , double r, double root_dv, double rh);
	double greens (double r, double root_dv);
	void   create_vLUT();
	void   createGLUT();
	
	static double G(double eps, double E_gamma, double boost);
	static double IC_cross_Section( double eps, double E_gamma, double E);
	static double CMB_bbSpectrum(double eps);
	static double dpIC(double eps, void * params );
	static double pIC(double nu, double E);

	static double djIC(double E , void * params);
 double jIC(double nu, double r);
	static double djsyn(double E , void * params);
 double jsyn(double r);

	static double dsIC( double r, void * params );
	double sIC( double nu);
	static double dssyn( double r, void * params );
	double ssyn();

};



