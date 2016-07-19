//*****greens.cpp*****//

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "Constants.h"
#include "Cluster.h"

double Cluster::D(double E){
	
	double D;
	
	if(DD == 0)
		D = D0 * pow(E, gamma);
	else if (DD == 1)
		D = D0 *pow( db , 2.0/3.0 )* pow(E, gamma)/pow(Bavg, 1.0/3.0);

	return D;
}

double Cluster::dv(double E , void * params){

	double dv = D(E)/bloss(E);

	return dv;
}

double Cluster::v( double E ){

	gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &dv;

	gsl_integration_qags (&F, E, p.mx, 0, 1e-8, 1000, 
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	result *= 1e16;
	return result;
}



double Cluster::root_dv(double Ep,  double vE){

	double root_dv  =   sqrt( ( vE ) - v(Ep)   ) ;

	return root_dv;
}

double Cluster::dgreens_term(double rp, void * params ){ 

	std::vector<double> greenParam = *(std::vector<double> *)params;
	double ri = greenParam[0];
	double r = greenParam[1];
	double root_dv = greenParam[2]; // 

	double dgreens_term = pow(root_dv , -1) * rp/ri * (exp( - pow( (rp-ri)/(2*root_dv) , 2)) 
		- exp( - pow( ( rp + ri)/(2*root_dv) , 2)) ) * pow( DM_profile(rp),2)/pow( DM_profile(r),2);
	

	return dgreens_term;

}

double Cluster::greens_term(double ri , double r, double root_dv, double Rh){


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> greenParam (3);


	greenParam[0] = ri;
	greenParam[1] = r;
	greenParam[2] = root_dv;

	gsl_function F;
	F.function = &dgreens_term;
	F.params = &greenParam;
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, 1e-16, Rh, 0, 1e-3, 1000, //x?
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}

double Cluster::greens (double r, double root_dv) {  //called by ddsyn

	double Rh = rh * kpc2cm ;
	//int imNum = 7; //number of image charges = 2*imNum + 1
	double Gsum = 0 ;
	for (int i = - imNum; i < imNum + 1; ++i ){

		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*Rh);

		Gsum += pow(-1, i) * greens_term(ri, r, root_dv, Rh);
	}

	double Greens = pow(4*pi , -1.0/2.0)*Gsum ;

	return Greens;

}

void Cluster::create_vLUT(){
	// iteration timer start
	std::clock_t vstart;
	double vduration;
	vstart = std::clock();
	///////before algorithm

	std::cout << "creating LUT..." <<std::endl; 
	
	for (int j = 0 ; j < vsize ; ++j ){

		vscale = p.mx/vsize;
		
		vlookup[j] = v(j*vscale);
		//std::cout << j<<"  " << sqrt(vlookup[j])/kpc2cm <<std::endl;
	}

	////////after algorithm
	vduration = (std::clock()  -  vstart)/(double) CLOCKS_PER_SEC;
	std::cout << "vlookup time = " << vduration <<std::endl;
}

void Cluster::createGLUT(){
	// iteration timer start
	std::clock_t Gstart;
	double Greensduration;
	Gstart = std::clock();
	///////before algorithm

	std::cout << "creating GLUT..." << std::endl; 


	double Rh = rh*kpc2cm;
	double r_scale = Rh/n_r;

	double rdv = sqrt(v(me));
	rootdv_max = (int)(rdv/kpc2cm + 1);
	std::cout << rootdv_max<<std::endl;

	double rootdv_scale = rootdv_max*kpc2cm/n_rootdv;

	GLUT[0][0] = 1;
	
	for (int i = 1 ; i < n_r ; ++i ){
																																																																																								
		for(int j = 1; j < n_rootdv ; ++j){


			GLUT[i][j] = greens(i*r_scale, j*rootdv_scale);

			GLUT[0][j] = 0;
			GLUT[i][0] = 1;
			GLUT[n_r - 1][j] = 0;

		}

		int a = (int)(  (kpc2cm*0.5)/rootdv_scale  );
		std::cout << i << "/" << n_r - 1 << " " << GLUT[i][a] << std::endl;
	};
	////////after algorithm
	Greensduration = (std::clock()  -  Gstart)/(double) CLOCKS_PER_SEC;
	std::cout << "GLUT time = " << Greensduration <<std::endl;

}

