#ifndef VEL_MODULE_H_
#define VEL_MODULE_H_

#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<sys/time.h>
#include <vector>
#include <fstream>      // std::ofstream
#include <iostream>
#include <iomanip>


#include "real_type.h"
#include "usekokkos.hpp"
#include "params.h"
#include "utils.hpp"


class vel_module {
public:
	int nx;
	int ny;
	real_t dt;
	real_t dx;
	char* filename;
	real_t* vel;
	real_t* courant;
	real_t* courant2;

	vel_module(Params& params):
	nx(params.NX), 
	ny(params.NY), 
	dt(params.dt), 
	dx(params.dx){

		vel  = new real_t[nx*ny];
		for (int i = 0; i < nx * ny; i++) {
    		vel[i] = 0; // set 0 for initiation
		}
		courant  = new real_t[nx*ny];
		courant2  = new real_t[nx*ny];
		for (int i = 0; i < nx * ny; i++) {
			courant[i] = 0; // set 0 for initiation
    		courant2[i] = 0; // set 0 for initiation
		}
	};

	void readbin_vel(char* filename);
	void default_vel();
	void set_courant();


	virtual ~vel_module() {
		delete[] vel;
		delete[] courant;
		delete[] courant2;
	};
};	// class vel_module


void vel_module::default_vel(){
	for (int i = 0; i < nx * ny; i++) {
    		vel[i] = 2500.0; // set uniform 6000 m/s
		}
}

//compute courant
void vel_module::set_courant(){
	for (int i = 0; i < nx * ny; i++) {
		courant[i] = vel[i]*dt/dx;
		courant2[i] = pow(courant[i],2);

		if (courant[i]>=1) {
			printf("Error for CFL > 1\n Exit");
			exit (EXIT_FAILURE);
		}
	}
	#ifdef DEBUG
 	printf("[vel_module] CFL MIN=%10.5e \n[vel_module] CFL MAX=%10.5e\n",
 		sqrt(FindMin<real_t>(courant,nx*ny)),sqrt(FindMax<real_t>(courant,nx*ny)));
 	#endif
}




// reading velocity structure file
void vel_module::readbin_vel(char* filename){

	#ifdef DEBUG
	printf("[vel_module] NX=%d NY=%d\n",nx,ny);
	#endif

	// std::ofstream outfile;
 // 	outfile.open(filename);
 // 	outfile.write((char*) vel, nx*ny*sizeof(real_t));  // write to file outfile
 // 	outfile.close();

	if (filename != nullptr) {
		#ifdef DEBUG
		printf("[vel_module] input file is %s\n",filename);
		#endif

 		std::ifstream infile;
 		infile.open(filename);
 		infile.read((char*) vel, nx*ny*sizeof(real_t));  // read from file infile
 		infile.close();
 	} else {
 		#ifdef DEBUG
		printf("[vel_module] no input file; use default vel structure\n");
		#endif
 		default_vel();
 	}

 	#ifdef DEBUG
 	printf("[vel_module] VEL MIN=%10.5e m/s\n[vel_module] VEL MAX=%10.5e m/s\n",
 		FindMin<real_t>(vel,nx*ny),FindMax<real_t>(vel,nx*ny));
 	#endif

 	// compute CFL at each grid
 	set_courant();
} // read file for vel structure



#endif
