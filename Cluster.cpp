//***** Cluster.cpp *****//

#include "Constants.h"
#include "Cluster.h"

double Cluster::z = 0.0232; 						
double Cluster::rh = 1000.0;						
double Cluster::rcore = 291.0;

double Cluster::B0 = 5.0;
double Cluster::Bmu = 2.0;
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

double Cluster::gamma = 0.3;
double Cluster::db = 20;
double Cluster::D0 = 3e28;	


Cluster::Cluster() :
				name("Default"),
				SD(1),			
				vsize(1000000),					
				vlookup(vsize),
				vscale(1),					//calculated in create_vLUT()

				imNum(7),
				n_r((int)(rh) + 1),
				n_rootdv(1000 + 1),
				rootdv_max(100),			//calculated in createGLUT()
				GLUT( n_r , std::vector<double>(n_rootdv) )

	{	//everything labelled or Coma, should work in user options
		std::cout << "creating cluster... " << std::endl;
	}

double Cluster::bfield_model (double r) {

	double rc = Cluster::rcore * kpc2cm;

	double B_field = B0 * pow(( 1 + r*r/(rc*rc)),(-1.5*beta*eta));		// Storm et al 2013 

	return B_field;

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

	double bloss = bsynch*pow(Bmu, 2.0)*E*E 								//bsyn bfield_model(r)
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

double Cluster::D(double E){

	double D = D0 *pow( db , 2.0/3.0 )* pow(E, gamma)/pow(Bmu, 1.0/3.0);

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
	int imNum = 7; //number of image charges = 2*imNum + 1
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

		//std::cout << "rdv " << root_dv/kpc2cm << std::endl;
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

		int a = (int)(  (kpc2cm*0.1)/rootdv_scale  );
		std::cout << i << "/" << n_r - 1 << " " << GLUT[i][a] << std::endl;
	};
	////////after algorithm
	Greensduration = (std::clock()  -  Gstart)/(double) CLOCKS_PER_SEC;
	std::cout << "GLUT time = " << Greensduration <<std::endl;

}

