#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//***** Functions.h *****//

/*


#include "Cluster.h"
#include "Particle.h"

double ddist(double z  , void * params);
double dist();
double rconst(double rcm);

double v( double E );
double root_dv(double Ep,  double vE);
double dgreens_term(double rp, void * params );
double greens_term(double ri , double r, double root_dv, double rh);
double greens (double r, double root_dv);
void   create_vLUT();
void   createGLUT();

extern"C" {	 							//interface with fortran code to initialize darksusy

void dsinit_();

}	

extern"C" {			
										// dshayield gives injection spectrum
double dshayield_(double *mwimp, double *emuthr, int *ch,  int *yieldk, int *istat);

}

double darksusy(double Ep);

double ddiffusion(double Ep, void * params);
double diffusion( double E, double r);
double elec_spect(double E, double r ); */

/////uses c.z and c.bfield_model
double fff(double x);
double dpsyn(double theta, void * params );
double psyn(  double E, double r);


/////no use of Cluster or particle class
double G(double eps, double E_gamma, double boost);
double IC_cross_Section( double eps, double E_gamma, double E);
double CMB_bbSpectrum(double eps);
double dpIC(double eps, void * params );
double pIC(double nu, double E);

/////uses p.mx
double djIC(double E , void * params);
double jIC(double nu, double r);
double djsyn(double E , void * params);
double jsyn(double r);

//////uses c.z c.DM_profile p.sv p.mx
double dsIC( double r, void * params );
double sIC( double nu, double r );
double dssyn( double r, void * params );
double ssyn( double r );

double min_flux(double r);
double calc_sv();
double calc_sv(double Sout);

void   runExCurve(int ch, double flux);

#endif