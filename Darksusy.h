//****** Darksusy.h *****//
extern"C" {	 							//interface with fortran code to initialize darksusy

void dsinit_();

}	

extern"C" {			
										// dshayield gives injection spectrum
double dshayield_(double *mwimp, double *emuthr, int *ch,  int *yieldk, int *istat);

}

double darksusy(double Ep);