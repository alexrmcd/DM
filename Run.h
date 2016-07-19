//*****Run.h*****//
#include "Cluster.h"


void runGreens_rdv(double r);
void runGreens_r(double rdv);

double min_flux(double r);
double calc_sv();
double calc_sv(double flux_obs );

void runExCurve(int ch, double flux_obs);
void runSED(double mx);
void runFlux(double mx);

