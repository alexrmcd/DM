#ifndef PARTICLE_H
#define PARTICLE_H

//***** Particle.h *****//

class Particle {
	public:

	int ch;
	double mx;
	double sv;

	Particle() : ch(25) , mx(1000), sv(4.7e-25){}
};

#endif