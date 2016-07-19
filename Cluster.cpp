//***** Cluster.cpp *****//

#include "Constants.h"
#include "Cluster.h"
std::string Cluster::name = "default";
double Cluster::z = 0.0232; 						
double Cluster::rh = 1000.0;						
double Cluster::rcore = 291.0;

double Cluster::B0 = 5.0;
double Cluster::Bavg = Cluster::Bmu();
double Cluster::beta = 0.75;
double Cluster::eta = 0.5;

int    Cluster::DM = 0;
double Cluster::rhos_NFW = 0.039974;		
double Cluster::rs_NFW = 404.0;	
double Cluster::rhos_Ein = 0.08296;
double Cluster::rs_Ein = 0.28;
double Cluster::alpha = 0.17;

double Cluster::bsynch = 0.0253;
double Cluster::bIC = 0.265;
double Cluster::bbrem = 1.51;
double Cluster::bcoul = 6.13;

int    Cluster::SD = 1;
int    Cluster::DD = 1;
double Cluster::gamma = 0.3;
double Cluster::db = 20;
double Cluster::D0 = 3e28;	

int    Cluster::vsize = 1000000;					
std::vector<double> Cluster::vlookup (vsize);
double Cluster::vscale = 1;					//calculated in create_vLUT()

int    Cluster::imNum = 7;
int    Cluster::n_r = 1001;
int    Cluster::n_rootdv = 1001;
int    Cluster::rootdv_max = 100;			//calculated in createGLUT()
std::vector< std::vector<double> >  Cluster::GLUT( n_r , std::vector<double>(n_rootdv) );

Particle Cluster::p;


Cluster::Cluster()

	{
		std::cout << "creating cluster... " << std::endl;
	}

double Cluster::bfield_model (double r) {

	double rc = rcore * kpc2cm;

	double B_field = B0 * pow(( 1 + r*r/(rc*rc)),(-1.5*beta*eta));		// Storm et al 2013 

	return B_field;

}

double Cluster::bfield_model (double r, void * params) {

	double rc = rcore * kpc2cm;

	double B_field = B0 * pow(( 1 + r*r/(rc*rc)),(-1.5*beta*eta));		// Storm et al 2013 

	return B_field;

}

double Cluster::Bmu(){
	
	double Rh = rh*kpc2cm;

	gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &bfield_model;

	gsl_integration_qags (&F, 1e-16, Rh, 0, 1e-8, 1000, 
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= 1.0/(Rh);
	
	return result;
}

double Cluster::bloss(double E , double r){

	double bloss = bsynch*pow(bfield_model(r), 2)*E*E 					//bsyn bfield_model(r)
					+ bIC * pow(1 + z, 4 )*E*E  					//bIC
					+ bbrem*nele*(0.36 + log(E/me/nele) )			//+ 1.51*n*(0.36 + log(E/me) )*E						//bbrem Note this is most likely incorrect, no factor of E, and should be E/me/nin log
					+ bcoul*nele*( 1 + log(E/me/nele)/75); 					//bcoul

	bloss *=1e-16;	

	return bloss;
}


double Cluster::bloss(double E ){ 												//overload so that we can just have bloss(E), could also have same as b(E, r) but set r=0 ??  kinda sloppy

	double bloss = bsynch*pow(Bavg, 2.0)*E*E 								//bsyn bfield_model(r)
					+ bIC /** pow(1 + z, 4 )*/*E*E  					//bIC
					+ bbrem*nele*(0.36 + log(E/me/nele) )						//brem , in Emma's code  has + 1.51*n*(0.36 + log(E/me) )*E
					+ bcoul*nele*( 1 + log(E/me/nele)/75); 					//bcoul

	//bloss *=1e-16;	

	return bloss;
}

double Cluster::DM_profile(double r){
	double rho;
	if (DM == 0){
		//NFW
		double rs = rs_NFW * kpc2cm;
		rho = rhos_NFW / ( r/rs * pow(1 + r/rs , 2 )) ; 
	
	}
	else if (DM == 1){
		// Einasto
		double rs = rs_Ein * kpc2cm;
		rho = (rhos_Ein) * exp(-2.0/alpha * ( pow(r/rs, alpha) - 1 ));
	}

	return rho;

}

