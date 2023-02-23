#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<sys/time.h>
#include <vector>
// #include <fstream>
// #include <iostream>
// #include<iomanip>

#include "real_type.h"
#include "params.h"
#include "usekokkos.hpp"
#include "useserial.hpp"
// vel_module
#include "vel_module.hpp"
#include "src_module.hpp"
// utils
#include "utils.hpp"

#ifdef KOKKOS_ENABLE_CUDA
#include "CudaTimer.h"
using Timer = CudaTimer;
#elif defined( KOKKOS_ENABLE_HIP)
#include "HipTimer.h"
using Timer = HipTimer;
#elif defined(KOKKOS_ENABLE_OPENMP)
#include "OpenMPTimer.h"
using Timer = OpenMPTimer;
#else
#include "SimpleTimer.h"
using Timer = SimpleTimer;
#endif
using DataArray2d = Kokkos::View<real_t**, Device>;
using DataArray1d = Kokkos::View<real_t*, Device>;
real_t serial(Params par, DataContext dataserial);
real_t kokkos_0(Params par, DataContext dataserial);
real_t kokkos_1(Params par, DataContext dataserial);
real_t kokkos_2(Params par, DataContext dataserial);
real_t kokkos_3(Params par, DataContext dataserial);
real_t kokkos_4(Params par, DataContext dataserial);
real_t kokkos_5(Params par, DataContext dataserial);

int main(int argc, char* argv[]) {
	int dim[150];
	int temp = 50;
	int timeStamp = 5; 
	int row = 6; 
	for(int i = 0; i < 150; i++)
	{
		dim[i] = temp;
		temp+=500; 
	}
	Kokkos::initialize(argc,argv);
	// real_t time_serial[timeStamp];
	// real_t time_kokkos_0[timeStamp];  //implemtxxxxx
	// real_t time_kokkos_1[timeStamp];


	real_t haha[row][timeStamp];

	for (int id=0; id<timeStamp; id++){ // diffrent dimension
		Params par(dim[id],dim[id], 1,0.005,40,40);
		DataContext       dataserial(par);
		//serial time
		// haha[0][id]=serial(par, dataserial);
		// std::cout << "serial at " << id << " is : " << haha[0][id] << std::endl; 

		//kokkos - verison 1 time
		haha[1][id]=kokkos_0(par, dataserial);
		std::cout << "Kokkos_0 at " << id << " is : " << haha[0][id] << std::endl;

		//kokkos - verison 1 time
		// haha[2][id]=kokkos_1(par, dataserial);
		// std::cout << "Kokkos_1 at " << id << " is : " << haha[2][id] << std::endl; 
		
		// // kokkos - version 2
		// haha[3][id]=kokkos_2(par, dataserial);
		//std::cout << "Kokkos_2 at " << id << " is : " << haha[3][id] << std::endl; 

		// haha[4][id]=kokkos_3(par, dataserial);
		// std::cout << "Kokkos_3 at " << id << " is : " << haha[4][id] << std::endl; 

		//kokkos - verison 4 time
		// haha[5][id]=kokkos_4(par, dataserial);
		// std::cout << "Kokkos_4 at " << id << " is : " << haha[5][id] << std::endl;
	}
	Kokkos::finalize();

 
	std::ofstream ofs ("plot_test_speed.py", std::ofstream::out);
 	ofs << "import numpy as np\n";
  	ofs << "import matplotlib.pyplot as plt\n";
 	ofs << "from matplotlib import rc\n";
  	ofs << "#rc('text', usetex=True)\n\n";
	ofs << "size=np.array([";
  	// output size array
	for (int i=0; i<timeStamp; ++i){
		if(i<timeStamp-1)
			ofs << dim[i] << "," ;
		else
			ofs << dim[i];
	}
	ofs << "])\n\n";
	for (int iv=0; iv<row; ++iv) {
		ofs << "v" << iv << "=np.array([";
		for (int i=0; i<timeStamp; ++i)
			i<timeStamp-1 ? ofs << haha[iv][i]<< "," : ofs << haha[iv][i];
		ofs << "])\n\n";
  	}
  ofs << "plt.plot(size,v0, label='serial')\n";
  ofs << "plt.plot(size,v1, label='kokkod with struct')\n";
  ofs << "plt.plot(size,v2, label='kokkod half flat half struct')\n";
  ofs << "plt.plot(size,v3, label='kokkod using struct with 2d range')\n";
  ofs << "plt.plot(size,v4, label='kokkod flat without 2d range')\n";
  ofs << "plt.plot(size,v5, label='kokkod flat with 2d range')\n";
  ofs << "plt.grid(True)\n";

  ofs << "plt.title('Size vs Time')\n";
  ofs << "plt.xlabel('Domain size')\n";
  ofs << "plt.ylabel(r'Time')\n";
  
  ofs << "plt.legend()\n";
  ofs << "plt.show()\n";
  
  ofs.close();
}

// real_t serial(Params par, DataContext dataserial){
// 	char sprintfBuffer[500];
// 	int fileIndex = 0; 
// 	Timer timer;

// 	int BC1[2]={1,2}; // BC: 0 (rigid); 1 (free surface); 2 (absorbing);
//     int BC2[2]={2,2};
//     vel_module vm(par);
//     vm.readbin_vel((char*)"in/randomvel.bin");
  
//     //Generate source 
//     src_module src(par,1);

//     for (int i=0; i<par.NX*par.NY; ++i){
//     	dataserial.u0[i] = 0.;
//     	dataserial.u1[i] = 0.;
//     	dataserial.u2[i] = 0.;
//     } 

// 	timer.start();
// 	for (int it=0; it < par.NT; ++it){
//     	// computing time
//     	par.Time = it * par.dt;
// 		for (int i=1; i<par.NX-1; ++i){
//     		for (int j=1; j<par.NY-1; ++j){
//     			const int ij = INDEX(i,j,par.NX,par.NY);
// 				const int ip1j = INDEX(i+1,j,par.NX,par.NY);
// 				const int im1j = INDEX(i-1,j,par.NX,par.NY);
// 				const int ijp1 = INDEX(i,j+1,par.NX,par.NY);
// 				const int ijm1 = INDEX(i,j-1,par.NX,par.NY);

// 				dataserial.u2[ij] = (2-4*vm.courant2[ij])*dataserial.u1[ij]+
// 				vm.courant2[ij]*
// 				(dataserial.u1[ip1j]+dataserial.u1[im1j]+
// 				 dataserial.u1[ijp1]+dataserial.u1[ijm1])-dataserial.u0[ij];
//     		}
//     	}
// 		src.add_src(par.Time);
// 		for (int i=0; i<par.NS; ++i){
//     		if (par.Time < src.tlen[i]) {
//     			int ind = INDEX(src.isx[i],src.isy[i],par.NX,par.NY);
//     			int ind2 = INDEX(src.isx[i],src.isy[i],par.NX,par.NY);
//     			dataserial.u2[ind2] = src.stf[i];
//     			// printf("%d %d %d %d\n", src.isx[i],src.isy[i],ind,ind2);
//     		}	
//     	}
// 		for (int i=0; i<par.NX; ++i){
//     		for (int j=0; j<par.NY; ++j){

//     			if (BC1[0]==2 and i==0) {
//     				const int ij = INDEX(i,j,par.NX,par.NY);
//     				const int i2=INDEX(1,j,par.NX,par.NY);
//     				const real_t c = vm.courant[i2];
//     				dataserial.u2[ij] = c*dataserial.u1[i2]+(1-c)*dataserial.u1[ij];
//     			}


//     			if (BC1[1]==2 and i==par.NX-1) {
//     				const int ij = INDEX(i,j,par.NX,par.NY);
//     				const int i2=INDEX(par.NX-2,j,par.NX,par.NY);
//     				const real_t c = vm.courant[i2];
//     				dataserial.u2[ij] = c*dataserial.u1[i2]+(1-c)*dataserial.u1[ij];
//     			}

//     			if (BC2[0]==2 and j==0) {
//     				const int ij = INDEX(i,j,par.NX,par.NY);
//     				const int i2=INDEX(i,1,par.NX,par.NY);
//     				const real_t c = vm.courant[i2];
//     				dataserial.u2[ij] = c*dataserial.u1[i2]+(1-c)*dataserial.u1[ij];
//     			}


//     			if (BC2[1]==2 and j==par.NY-1) {
//     				const int ij = INDEX(i,j,par.NX,par.NY);
//     				const int i2=INDEX(i,par.NY-2,par.NX,par.NY);
//     				const real_t c = vm.courant[i2];
//     				dataserial.u2[ij] = c*dataserial.u1[i2]+(1-c)*dataserial.u1[ij];
//     			}


//     			if (BC1[0]==1 and i==0) {
//     				const int ij = INDEX(i,j,par.NX,par.NY);
//     				const int i2=INDEX(1,j,par.NX,par.NY);
//     				dataserial.u2[ij] = dataserial.u2[i2];
//     			}


//     			if (BC2[0]==1 and i==par.NX-1) {
//     				const int ij = INDEX(i,j,par.NX,par.NY);
//     				const int i2=INDEX(par.NX-2,j,par.NX,par.NY);
//     				dataserial.u2[ij] = dataserial.u2[i2];
//     			}

//     			if (BC1[1]==1 and j==0) {
//     				const int ij = INDEX(i,j,par.NX,par.NY);
//     				const int i2=INDEX(i,1,par.NX,par.NY);
//     				dataserial.u2[ij] = dataserial.u2[i2];
//     			}


//     			if (BC2[1]==1 and j==par.NY-1) {
//     				const int ij = INDEX(i,j,par.NX,par.NY);
//     				const int i2=INDEX(i,par.NY-2,par.NX,par.NY);
//     				dataserial.u2[ij] = dataserial.u2[i2];
//     			}


//     		}
// 		}
// 		for (int i=0; i<par.NX*par.NY; ++i){
//     		dataserial.u0[i] = dataserial.u1[i];
//     		dataserial.u1[i] = dataserial.u2[i];
//     	}
// 		// sprintf(sprintfBuffer, "out/Serial_%03u.csv", fileIndex);
//       	// FILE* file2 = fopen(sprintfBuffer, "w");
// 		// int ind=0;
//         // for (unsigned int i = 0; i < par.NX; i += 1) {
// 		// 	ind = INDEX(i,0,par.NX,par.NY);
//         // 	fprintf(file2, "%4.2f", dataserial.u1[ind]);
//     	// 	for (unsigned int j = 1; j < par.NY; j += 1) {
//     	// 		ind = INDEX(i,j,par.NX,par.NY);
//         // 		fprintf(file2, ", %4.2f", dataserial.u1[ind]);
//         // 	}
//         // 		fprintf(file2, "\n");
//       	// 	}
// 		// fclose(file2);
// 		// ++fileIndex; 

// 	}
// 	timer.stop();
// 	double time_seconds = timer.elapsed();
// 	return time_seconds;
// }
real_t kokkos_0(Params par, DataContext dataserial)
{
	char sprintfBuffer[500];
	int fileIndex = 0; 
	Timer timer;
	int BC1[2]={1,2}; //BC: 0 (rigid); 1 (free surface); 2 (absorbing);
    int BC2[2]={2,2};
    vel_module vm(par);
    vm.readbin_vel((char*)"in/randomvel.bin");
	src_module src(par,1);
	DataContextKokkos datakokkos(par);
	for (int i=0; i<par.NX*par.NY; ++i){
    	datakokkos.c2_h(i) = vm.courant2[i];
    }
    Kokkos::deep_copy(datakokkos.c2_d,datakokkos.c2_h);
	for (int i=0; i<par.NX*par.NY; ++i){
    	datakokkos.c_h(i) = vm.courant[i];
    }
    Kokkos::deep_copy(datakokkos.c_d,datakokkos.c_h);
	Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) {
    	datakokkos.u0_d(index) = 0.;
    	datakokkos.u1_d(index) = 0.;
    	datakokkos.u2_d(index) = 0.;
    });

    Kokkos::fence();


    for (int i=0; i<par.NX*par.NY; ++i){
    	datakokkos.u_h(i) = 0.;
	}
	fileIndex = 0;
	timer.start();
    for (int it=0; it < par.NT; ++it){
		// timer.start();
    	// computing time
    	par.Time = it * par.dt;
		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) 
		{  
    		int i,j;
			index2coord(index,i,j,par.NX,par.NY);

			const int ij = INDEX(i,j,par.NX,par.NY);
			const int ip1j = INDEX(i+1,j,par.NX,par.NY);
			const int im1j = INDEX(i-1,j,par.NX,par.NY);
			const int ijp1 = INDEX(i,j+1,par.NX,par.NY);
			const int ijm1 = INDEX(i,j-1,par.NX,par.NY);

			// if (index == 100) printf("%03d %03d %010d vs ",i,j,ij);
			if (i>0 and i<par.NX-1 and j>0 and j<par.NY-1 ){
				datakokkos.u2_d(ij) = (2-4*datakokkos.c2_d(ij))*datakokkos.u1_d(ij)+
				datakokkos.c2_d(ij)*
				(datakokkos.u1_d(ip1j)+datakokkos.u1_d(im1j)+
				 datakokkos.u1_d(ijp1)+datakokkos.u1_d(ijm1))-datakokkos.u0_d(ij);
			}

    	});
		src.add_src(par.Time);
    	//deep copy n source from Host to Device
    	for (int i=0; i<par.NS; ++i){
    		datakokkos.src_h(i) = src.stf[i];
    		// if (it<3) printf("%d %f\n",it,src.stf[i]);
    	}

    	Kokkos::deep_copy(datakokkos.src_d,datakokkos.src_h);
    	Kokkos::fence();

    	for (int i=0; i<par.NS; ++i){
    		if (par.Time < src.tlen[i]) {
    			int ind = INDEX(src.isx[i],src.isy[i],par.NX,par.NY);

    			Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
    				if(index ==  ind) datakokkos.u2_d(index)=datakokkos.src_d(i);
    			});
    			// printf("%d %d %d %d\n", src.isx[i],src.isy[i],ind,ind2);
    		}	
    	}
		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
    		int i,j,i2;
			index2coord(index,i,j,par.NX,par.NY);

			// // Absorbing boundary on the top edge (0-x)
			if (BC1[0]==2 and i==0) {
				i2=INDEX(1,j,par.NX,par.NY);
				const real_t c = datakokkos.c_d(i2);
				datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+
										(1-c)*datakokkos.u1_d(index);
			}

			// // Absorbing boundary on the bottom edge (x-1)
			if (BC2[0]==2 and i==par.NX-1) {
				i2=INDEX(par.NX-2,j,par.NX,par.NY);
				const real_t c = datakokkos.c_d(i2);
				datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+
										(1-c)*datakokkos.u1_d(index);
			}


			// // Absorbing boundary on the left edge (0-y)
			if (BC1[1]==2 and j==0) {
				i2=INDEX(i,1,par.NX,par.NY);
				const real_t c = datakokkos.c_d(i2);
				datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+
										(1-c)*datakokkos.u1_d(index);
			}

			// // Absorbing boundary on the right edge (y-1)
			if (BC2[1]==2 and j==par.NY-1) {
				i2=INDEX(i,par.NY-2,par.NX,par.NY);
				const real_t c = datakokkos.c_d(i2);
				datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+
										(1-c)*datakokkos.u1_d(index);
			}


			// below are free surface BC
			if (BC1[0]==1 and i==0) {
				i2=INDEX(1,j,par.NX,par.NY);
				datakokkos.u2_d(index) = datakokkos.u2_d(i2);
			}


			if (BC2[0]==1 and i==par.NX-1) {
				i2=INDEX(par.NX-2,j,par.NX,par.NY);
				datakokkos.u2_d(index) = datakokkos.u2_d(i2);
			}


			if (BC1[1]==1 and j==0) {
				i2=INDEX(i,1,par.NX,par.NY);
				datakokkos.u2_d(index) = datakokkos.u2_d(i2);
			}


			if (BC2[1]==1 and j==par.NY-1) {
				i2=INDEX(i,par.NY-2,par.NX,par.NY);
				datakokkos.u2_d(index) = datakokkos.u2_d(i2);
			}
    	});
		Kokkos::deep_copy(datakokkos.u0_d,datakokkos.u1_d);
    	Kokkos::deep_copy(datakokkos.u1_d,datakokkos.u2_d);
		// Kokkos::deep_copy(datakokkos.u_h, datakokkos.u1_d);

		// sprintf(sprintfBuffer, "out/Kokkos_%03u.csv", fileIndex);
      	// FILE* file = fopen(sprintfBuffer, "w");
      	// if (!file){
    	// 	std::cout << "Couldn't write the file!\nPlease 'mkdir out/' and rerun" << std::endl;
		// 	return 0;
		// }
		// int ind=0;
      	// for (unsigned int i = 0; i < par.NX; i += 1) {
      	// 	ind = INDEX(i,0,par.NX,par.NY);
		// 	fprintf(file, "%4.2f", datakokkos.u_h(ind));
        // 	for (unsigned int j = 1; j < par.NY; j += 1) {
        // 		ind = INDEX(i,j,par.NX,par.NY);
        // 		fprintf(file, ", %4.2f", datakokkos.u_h(ind));
        // 	}
        // 	fprintf(file, "\n");
      	// }

      	// fclose(file);
      	// ++fileIndex;
	}
	timer.stop();
	#ifdef DEBUG
		std::cout << "debug begin" << std::endl; 
    	for (int index=0; index<par.NX*par.NY; ++index){
    		if (abs(datakokkos.u_h(index)-dataserial.u1[index])>1e-5) 
					printf("U1host= %f serial= %f\n",datakokkos.u_h(index),dataserial.u1[index]);
		}
		std::cout << "debug end" << std::endl; 
    #endif
	double time_seconds = timer.elapsed();
	return time_seconds;
}

//2d Range
// real_t kokkos_1(Params par, DataContext dataserial)
// {
// 	char sprintfBuffer[500];
// 	int fileIndex = 0; 
// 	Timer timer;
//     DataArray1d BC1("BC1",2);
//   	DataArray1d BC2("BC2",2);
// 	BC1(0) = 1;
// 	BC1(1) = 2;
// 	BC2(0) = 2; 
// 	BC2(1) = 2; 

//     vel_module vm(par);
//     vm.readbin_vel((char*)"in/randomvel.bin");
// 	src_module src(par,1);
// 	DataContextKokkos datakokkos(par);
// 	for (int i=0; i<par.NX*par.NY; ++i){
//     	datakokkos.c2_h(i) = vm.courant2[i];
//     }
//     Kokkos::deep_copy(datakokkos.c2_d,datakokkos.c2_h);
// 	for (int i=0; i<par.NX*par.NY; ++i){
//     	datakokkos.c_h(i) = vm.courant[i];
//     }
//     Kokkos::deep_copy(datakokkos.c_d,datakokkos.c_h);
// 	Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) {
//     	datakokkos.u0_d(index) = 0.;
//     	datakokkos.u1_d(index) = 0.;
//     	datakokkos.u2_d(index) = 0.;
//     });

//     Kokkos::fence();


//     for (int i=0; i<par.NX*par.NY; ++i){
//     	datakokkos.u_h(i) = 0.;
// 	}
// 	fileIndex = 0;
// 	timer.start();
//     for (int it=0; it < par.NT; ++it){
// 		// timer.start();
//     	// computing time
//     	par.Time = it * par.dt;
// 		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) 
// 		{  
//     		int i,j;
// 			index2coord(index,i,j,par.NX,par.NY);

// 			const int ij = INDEX(i,j,par.NX,par.NY);
// 			const int ip1j = INDEX(i+1,j,par.NX,par.NY);
// 			const int im1j = INDEX(i-1,j,par.NX,par.NY);
// 			const int ijp1 = INDEX(i,j+1,par.NX,par.NY);
// 			const int ijm1 = INDEX(i,j-1,par.NX,par.NY);

// 			// if (index == 100) printf("%03d %03d %010d vs ",i,j,ij);
// 			if (i>0 and i<par.NX-1 and j>0 and j<par.NY-1 ){
// 				datakokkos.u2_d(ij) = (2-4*datakokkos.c2_d(ij))*datakokkos.u1_d(ij)+
// 				datakokkos.c2_d(ij)*
// 				(datakokkos.u1_d(ip1j)+datakokkos.u1_d(im1j)+
// 				 datakokkos.u1_d(ijp1)+datakokkos.u1_d(ijm1))-datakokkos.u0_d(ij);
// 			}

//     	});
// 		src.add_src(par.Time);
//     	//deep copy n source from Host to Device
//     	for (int i=0; i<par.NS; ++i){
//     		datakokkos.src_h(i) = src.stf[i];
//     		// if (it<3) printf("%d %f\n",it,src.stf[i]);
//     	}

//     	Kokkos::deep_copy(datakokkos.src_d,datakokkos.src_h);
//     	Kokkos::fence();

//     	for (int i=0; i<par.NS; ++i){
//     		if (par.Time < src.tlen[i]) {
//     			int ind = INDEX(src.isx[i],src.isy[i],par.NX,par.NY);

//     			Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     				if(index ==  ind) datakokkos.u2_d(index)=datakokkos.src_d(i);
//     			});
//     			// printf("%d %d %d %d\n", src.isx[i],src.isy[i],ind,ind2);
//     		}	
//     	}
// 		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     		int i,j,i2;
// 			index2coord(index,i,j,par.NX,par.NY);

// 			// // Absorbing boundary on the top edge (0-x)
// 			if (BC1(0)==2 and i==0) {
// 				i2=INDEX(1,j,par.NX,par.NY);
// 				const real_t c = datakokkos.c_d(i2);
// 				datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+
// 										(1-c)*datakokkos.u1_d(index);
// 			}

// 			// // Absorbing boundary on the bottom edge (x-1)
// 			if (BC2(0)==2 and i==par.NX-1) {
// 				i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 				const real_t c = datakokkos.c_d(i2);
// 				datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+
// 										(1-c)*datakokkos.u1_d(index);
// 			}


// 			// // Absorbing boundary on the left edge (0-y)
// 			if (BC1(1)==2 and j==0) {
// 				i2=INDEX(i,1,par.NX,par.NY);
// 				const real_t c = datakokkos.c_d(i2);
// 				datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+
// 										(1-c)*datakokkos.u1_d(index);
// 			}

// 			// // Absorbing boundary on the right edge (y-1)
// 			if (BC2(1)==2 and j==par.NY-1) {
// 				i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 				const real_t c = datakokkos.c_d(i2);
// 				datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+
// 										(1-c)*datakokkos.u1_d(index);
// 			}


// 			// below are free surface BC
// 			if (BC1(0)==1 and i==0) {
// 				i2=INDEX(1,j,par.NX,par.NY);
// 				datakokkos.u2_d(index) = datakokkos.u2_d(i2);
// 			}


// 			if (BC2(0)==1 and i==par.NX-1) {
// 				i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 				datakokkos.u2_d(index) = datakokkos.u2_d(i2);
// 			}


// 			if (BC1(1)==1 and j==0) {
// 				i2=INDEX(i,1,par.NX,par.NY);
// 				datakokkos.u2_d(index) = datakokkos.u2_d(i2);
// 			}


// 			if (BC2(1)==1 and j==par.NY-1) {
// 				i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 				datakokkos.u2_d(index) = datakokkos.u2_d(i2);
// 			}
//     	});

// 		Kokkos::deep_copy(datakokkos.u0_d,datakokkos.u1_d);
//     	Kokkos::deep_copy(datakokkos.u1_d,datakokkos.u2_d);
// 		//timer.elapsed();
// 		// Kokkos::deep_copy(datakokkos.u_h, datakokkos.u1_d);
// 		// sprintf(sprintfBuffer, "out/Kokkos1_%03u.csv", fileIndex);
//       	// FILE* file = fopen(sprintfBuffer, "w");
//       	// if (!file){
//     	// 	std::cout << "Couldn't write the file!\nPlease 'mkdir out/' and rerun" << std::endl;
// 		// 	return 0;
// 		// }
// 		// int ind=0;
//       	// for (unsigned int i = 0; i < par.NX; i += 1) {
//       	// 	ind = INDEX(i,0,par.NX,par.NY);
// 		// 	fprintf(file, "%4.2f", datakokkos.u_h(ind));
//         // 	for (unsigned int j = 1; j < par.NY; j += 1) {
//         // 		ind = INDEX(i,j,par.NX,par.NY);
//         // 		fprintf(file, ", %4.2f", datakokkos.u_h(ind));
//         // 	}
//         // 	fprintf(file, "\n");
//       	// }

//       	// fclose(file);
//       	// ++fileIndex;
// 	}
// 	timer.stop();
// 	#ifdef DEBUG
// 		std::cout << "debug begin" << std::endl; 
//     	for (int index=0; index<par.NX*par.NY; ++index){
//     		if (abs(datakokkos.u_h(index)-dataserial.u1[index])>1e-5) 
// 					printf("U1host= %f serial= %f\n",datakokkos.u_h(index),dataserial.u1[index]);
// 		}
// 		std::cout << "debug end" << std::endl; 
//     #endif
// 	double time_seconds = timer.elapsed();
// 	return time_seconds;
// }
// // kokkos_2d range without view 
// real_t kokkos_2(Params par, DataContext dataserial){
// 	char sprintfBuffer[500];
// 	int fileIndex = 0; 
// 	Timer timer;
// 	int BC1[2]={1,2}; //BC: 0 (rigid); 1 (free surface); 2 (absorbing);
// 	int BC2[2]={2,2};
// 	vel_module vm(par);
// 	vm.readbin_vel((char*)"in/randomvel.bin");
// 	src_module src(par,1);
// 	DataContextKokkos datakokkos(par);
// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		datakokkos.c2_h(i) = vm.courant2[i];
// 	}
// 	Kokkos::deep_copy(datakokkos.c2_d,datakokkos.c2_h);
// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		datakokkos.c_h(i) = vm.courant[i];
// 	}
// 	Kokkos::deep_copy(datakokkos.c_d,datakokkos.c_h);
// 	Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) {
// 		datakokkos.u0_d(index) = 0.;
// 		datakokkos.u1_d(index) = 0.;
// 		datakokkos.u2_d(index) = 0.;
// 	});

// 	using Range2D = typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<2> >;
// 	Range2D range2d( {{0,0}}, {{par.NX,par.NY}} );
// 	Kokkos::fence();
// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		datakokkos.u_h(i) = 0.;
// 	}
// 	fileIndex = 0;
// 	timer.start();
// 	for (int it=0; it < par.NT; ++it){
// 		// timer.start();
// 		// computing time
// 		par.Time = it * par.dt;
// 	using Range2D = typename Kokkos::Experimental::MDRangePolicy
// 											< Kokkos::Experimental::Rank<2> >;
// 	Range2D range2d( {{0,0}}, {{par.NX,par.NY}} );
// 	Kokkos::parallel_for( range2d, KOKKOS_LAMBDA(const int& i, const int& j) 
// 	{  
// 		const int ij = INDEX(i,j,par.NX,par.NY);
// 		const int ip1j = INDEX(i+1,j,par.NX,par.NY);
// 		const int im1j = INDEX(i-1,j,par.NX,par.NY);
// 		const int ijp1 = INDEX(i,j+1,par.NX,par.NY);
// 		const int ijm1 = INDEX(i,j-1,par.NX,par.NY);

// 		// if (index == 100) printf("%03d %03d %010d vs ",i,j,ij);
// 		if (i>0 and i<par.NX-1 and j>0 and j<par.NY-1 ){
// 			datakokkos.u2_d(ij) = (2-4*datakokkos.c2_d(ij))*datakokkos.u1_d(ij)+
// 			datakokkos.c2_d(ij)*
// 			(datakokkos.u1_d(ip1j)+datakokkos.u1_d(im1j)+
// 			datakokkos.u1_d(ijp1)+datakokkos.u1_d(ijm1))-datakokkos.u0_d(ij);
// 		}

// 	});
// 	src.add_src(par.Time);
// 	//deep copy n source from Host to Device
// 	for (int i=0; i<par.NS; ++i){
// 		datakokkos.src_h(i) = src.stf[i];
// 		// if (it<3) printf("%d %f\n",it,src.stf[i]);
// 	}

// 	Kokkos::deep_copy(datakokkos.src_d,datakokkos.src_h);
// 	Kokkos::fence();

// 	for (int i=0; i<par.NS; ++i){
// 		if (par.Time < src.tlen[i]) {
// 		int ind = INDEX(src.isx[i],src.isy[i],par.NX,par.NY);

// 		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
// 			if(index ==  ind) datakokkos.u2_d(index)=datakokkos.src_d(i);
// 		});
// 		// printf("%d %d %d %d\n", src.isx[i],src.isy[i],ind,ind2);
// 		} 
// 	}
// 	Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     int i,j,i2;
//    	index2coord(index,i,j,par.NX,par.NY);

// 	// // Absorbing boundary on the top edge (0-x)
// 	if (BC1[0]==2 and i==0) {
// 		i2=INDEX(1,j,par.NX,par.NY);
// 		const real_t c = datakokkos.c_d(i2);
// 		datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+(1-c)*datakokkos.u1_d(index);
// 	}

// 	// // Absorbing boundary on the bottom edge (x-1)
// 	if (BC2[0]==2 and i==par.NX-1) {
// 		i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 		const real_t c = datakokkos.c_d(i2);
// 		datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+(1-c)*datakokkos.u1_d(index);
// 	}


// 	// // Absorbing boundary on the left edge (0-y)
// 	if (BC1[1]==2 and j==0) {
// 		i2=INDEX(i,1,par.NX,par.NY);
// 		const real_t c = datakokkos.c_d(i2);
// 		datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+(1-c)*datakokkos.u1_d(index);
// 	}

// 	// // Absorbing boundary on the right edge (y-1)
// 	if (BC2[1]==2 and j==par.NY-1) {
// 		i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 		const real_t c = datakokkos.c_d(i2);
// 		datakokkos.u2_d(index) = c*datakokkos.u1_d(i2)+(1-c)*datakokkos.u1_d(index);
// 	}


// 	// below are free surface BC
// 	if (BC1[0]==1 and i==0) {
// 		i2=INDEX(1,j,par.NX,par.NY);
// 		datakokkos.u2_d(index) = datakokkos.u2_d(i2);
// 	}


// 	if (BC2[0]==1 and i==par.NX-1) {
// 		i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 		datakokkos.u2_d(index) = datakokkos.u2_d(i2);
// 	}


// 	if (BC1[1]==1 and j==0) {
// 		i2=INDEX(i,1,par.NX,par.NY);
// 		datakokkos.u2_d(index) = datakokkos.u2_d(i2);
// 	}


// 	if (BC2[1]==1 and j==par.NY-1) {
// 		i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 		datakokkos.u2_d(index) = datakokkos.u2_d(i2);
// 	}
// 		});
// 	// timer.elapsed();
// 	Kokkos::deep_copy(datakokkos.u0_d,datakokkos.u1_d);
// 	Kokkos::deep_copy(datakokkos.u1_d,datakokkos.u2_d);
// 	// Kokkos::deep_copy(datakokkos.u_h, datakokkos.u1_d);
// 	// sprintf(sprintfBuffer, "out/Kokkos2_%03u.csv", fileIndex);
// 	// FILE* file = fopen(sprintfBuffer, "w");
// 	// if (!file){
// 	// 	std::cout << "Couldn't write the file!\nPlease 'mkdir out/' and rerun" << std::endl;
// 	// 	return 0;
// 	// }
// 	// int ind=0;
// 	// for (unsigned int i = 0; i < par.NX; i += 1) {
// 	// ind = INDEX(i,0,par.NX,par.NY);
// 	// fprintf(file, "%4.2f", datakokkos.u_h(ind));
// 	// for (unsigned int j = 1; j < par.NY; j += 1) {
// 	// 	ind = INDEX(i,j,par.NX,par.NY);
// 	// 	fprintf(file, ", %4.2f", datakokkos.u_h(ind));
// 	// }
// 	// fprintf(file, "\n");
// 	// }
// 	// fclose(file);
// 	// ++fileIndex;
// 	}
// 	timer.stop();
// 	#ifdef DEBUG
// 	std::cout << "debug begin" << std::endl; 
// 		for (int index=0; index<par.NX*par.NY; ++index){
// 		if (abs(datakokkos.u_h(index)-dataserial.u1[index])>1e-5) 
// 		printf("U1host= %f serial= %f\n",datakokkos.u_h(index),dataserial.u1[index]);
// 	}
// 	std::cout << "debug end" << std::endl; 
// 		#endif
// 	double time_seconds = timer.elapsed();
// 	return time_seconds;
// }
// // kokkos_1d array without 2d range
// real_t kokkos_3(Params par, DataContext dataserial)
// {
// 	char sprintfBuffer[500];
// 	int fileIndex = 0; 
// 	Timer timer;
// 	int BC1[2]={1,2}; //BC: 0 (rigid); 1 (free surface); 2 (absorbing);
//     int BC2[2]={2,2};
//     vel_module vm(par);
//     vm.readbin_vel((char*)"in/randomvel.bin");
// 	src_module src(par,1);
// 	DataArray1d c2h("c2h", par.NX*par.NY);
// 	DataArray1d c2d("c2d", par.NX*par.NY);
// 	DataArray1d ch("ch", par.NX*par.NY);
// 	DataArray1d cd("cd", par.NX*par.NY);

// 	DataArray1d srch("srch", par.NS);
// 	DataArray1d srcd("srcd", par.NS);

// 	DataArray1d u0d("u0d", par.NX *par.NY);
// 	DataArray1d u1d("u1d", par.NX*par.NY);
// 	DataArrayHost uh=Kokkos::create_mirror_view(u1d);
// 	DataArray1d u2d("u2d", par.NX*par.NY);

// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		c2h[i] = vm.courant2[i];
//     }
//     Kokkos::deep_copy(c2d, c2h);

// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		ch[i]= vm.courant[i];
//     }
// 	Kokkos::deep_copy(cd,ch);

// 	Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) {
// 		u0d[index] = 0.;
// 		u1d[index] = 0.;
// 		u2d[index] = 0.; 
//     });

//     Kokkos::fence();


//     for (int i=0; i<par.NX*par.NY; ++i){
// 		uh[i] = 0; 
// 	}
// 	fileIndex = 0;
// 	timer.start();
//     for (int it=0; it < par.NT; ++it){
// 		//timer.start();
//     	// computing time
//     	par.Time = it * par.dt;
// 		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) 
// 		{  
//     		int i,j;
// 			index2coord(index,i,j,par.NX,par.NY);

// 			const int ij = INDEX(i,j,par.NX,par.NY);
// 			const int ip1j = INDEX(i+1,j,par.NX,par.NY);
// 			const int im1j = INDEX(i-1,j,par.NX,par.NY);
// 			const int ijp1 = INDEX(i,j+1,par.NX,par.NY);
// 			const int ijm1 = INDEX(i,j-1,par.NX,par.NY);

// 			// if (index == 100) printf("%03d %03d %010d vs ",i,j,ij);
// 			if (i>0 and i<par.NX-1 and j>0 and j<par.NY-1 ){
// 				u2d(ij) = (2-4*c2d(ij))*u1d(ij)+
// 				c2d(ij)*
// 				(u1d(ip1j)+u1d(im1j)+
// 				 u1d(ijp1)+u1d(ijm1))-u0d(ij);
// 			}

//     	});
// 		src.add_src(par.Time);
//     	//deep copy n source from Host to Device
//     	for (int i=0; i<par.NS; ++i){
//     		srch(i) = src.stf[i];
//     		// if (it<3) printf("%d %f\n",it,src.stf[i]);
//     	}

//     	Kokkos::deep_copy(srcd, srch);
//     	Kokkos::fence();

//     	for (int i=0; i<par.NS; ++i){
//     		if (par.Time < src.tlen[i]) {
//     			int ind = INDEX(src.isx[i],src.isy[i],par.NX,par.NY);

//     			Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     				if(index ==  ind) u2d(index)=srcd(i);
//     			});
//     			// printf("%d %d %d %d\n", src.isx[i],src.isy[i],ind,ind2);
//     		}	
//     	}
// 		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     		int i,j,i2;
// 			index2coord(index,i,j,par.NX,par.NY);

// 			// // Absorbing boundary on the top edge (0-x)
// 			if (BC1[0]==2 and i==0) {
// 				i2=INDEX(1,j,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}

// 			// // Absorbing boundary on the bottom edge (x-1)
// 			if (BC2[0]==2 and i==par.NX-1) {
// 				i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}


// 			// // Absorbing boundary on the left edge (0-y)
// 			if (BC1[1]==2 and j==0) {
// 				i2=INDEX(i,1,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}

// 			// // Absorbing boundary on the right edge (y-1)
// 			if (BC2[1]==2 and j==par.NY-1) {
// 				i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}


// 			// below are free surface BC
// 			if (BC1[0]==1 and i==0) {
// 				i2=INDEX(1,j,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC2[0]==1 and i==par.NX-1) {
// 				i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC1[1]==1 and j==0) {
// 				i2=INDEX(i,1,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC2[1]==1 and j==par.NY-1) {
// 				i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}
//     	});
// 		Kokkos::deep_copy(u0d,u1d);
//     	Kokkos::deep_copy(u1d,u2d);
// 		// timer.elapsed();
// 		//Kokkos::deep_copy(uh, u1d);
// 		// sprintf(sprintfBuffer, "out/Kokkos3_%03u.csv", fileIndex);
//       	// FILE* file = fopen(sprintfBuffer, "w");
//       	// if (!file){
//     	// 	std::cout << "Couldn't write the file!\nPlease 'mkdir out/' and rerun" << std::endl;
// 		// 	return 0;
// 		// }
// 		// int ind=0;
//       	// for (unsigned int i = 0; i < par.NX; i += 1) {
//       	// 	ind = INDEX(i,0,par.NX,par.NY);
// 		// 	fprintf(file, "%4.2f", datakokkos.u_h(ind));
//         // 	for (unsigned int j = 1; j < par.NY; j += 1) {
//         // 		ind = INDEX(i,j,par.NX,par.NY);
//         // 		fprintf(file, ", %4.2f", datakokkos.u_h(ind));
//         // 	}
//         // 	fprintf(file, "\n");
//       	// }

//       	// fclose(file);
//       	// ++fileIndex;
// 	}
// 	timer.stop();
// 	#ifdef DEBUG
// 		std::cout << "debug begin" << std::endl; 
//     	for (int index=0; index<par.NX*par.NY; ++index){
//     		if (abs(uh(index)-dataserial.u1[index])>1e-5) 
// 					printf("U1host= %f serial= %f\n",datakokkos.u_h(index),dataserial.u1[index]);
// 		}
// 		std::cout << "debug end" << std::endl; 
//     #endif
// 	double time_seconds = timer.elapsed();
// 	return time_seconds;
// }

// real_t kokkos_4(Params par, DataContext dataserial)
// {
// 	char sprintfBuffer[500];
// 	int fileIndex = 0; 
// 	Timer timer;
// 	int BC1[2]={1,2}; //BC: 0 (rigid); 1 (free surface); 2 (absorbing);
//     int BC2[2]={2,2};
//     vel_module vm(par);
//     vm.readbin_vel((char*)"in/randomvel.bin");
// 	src_module src(par,1);
// 	DataContextKokkos datakokkos(par);
// 	DataArray1d c2h("c2h", par.NX*par.NY);
// 	DataArray1d c2d("c2d", par.NX*par.NY);

// 	DataArray1d ch("ch", par.NX*par.NY);
// 	DataArray1d cd("cd", par.NX*par.NY);
// 	DataArray1d srch("srch", par.NS);
// 	DataArray1d srcd("srcd", par.NS);

// 	DataArray1d u0d("u0d", par.NX *par.NY);
// 	DataArray1d u1d("u1d", par.NX*par.NY);
// 	DataArrayHost uh=Kokkos::create_mirror_view(u1d);

// 	DataArray1d u2d("u2d", par.NX*par.NY);

// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		c2h[i] = vm.courant2[i];
//     }
//     Kokkos::deep_copy(c2d, c2h);

// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		ch[i]= vm.courant[i];
//     }
// 	Kokkos::deep_copy(cd,ch);

// 	Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) {
// 		u0d[index] = 0.;
// 		u1d[index] = 0.;
// 		u2d[index] = 0.; 
//     });
	
//     Kokkos::fence();


//     for (int i=0; i<par.NX*par.NY; ++i){
// 		uh[i] = 0; 
// 	}
// 	fileIndex = 0;
// 	using Range2D = typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<2>>;
// 	Range2D range2d( {{0,0}}, {{par.NX,par.NY}} );

// 	timer.start();
//     for (int it=0; it < par.NT; ++it){
//     	// computing time
//     	par.Time = it * par.dt;
// 		Kokkos::parallel_for( range2d, KOKKOS_LAMBDA(const int& i, const int& j) 
// 		{  
// 			const int ij = INDEX(i,j,par.NX,par.NY);
// 			const int ip1j = INDEX(i+1,j,par.NX,par.NY);
// 			const int im1j = INDEX(i-1,j,par.NX,par.NY);
// 			const int ijp1 = INDEX(i,j+1,par.NX,par.NY);
// 			const int ijm1 = INDEX(i,j-1,par.NX,par.NY);

// 			// if (index == 100) printf("%03d %03d %010d vs ",i,j,ij);
// 			if (i>0 and i<par.NX-1 and j>0 and j<par.NY-1 ){
// 				u2d(ij) = (2-4*c2d(ij))*u1d(ij)+
// 				c2d(ij)*
// 				(u1d(ip1j)+u1d(im1j)+
// 				 u1d(ijp1)+u1d(ijm1))-u0d(ij);
// 			}

//     	});
// 		src.add_src(par.Time);
//     	//deep copy n source from Host to Device
//     	for (int i=0; i<par.NS; ++i){
//     		srch(i) = src.stf[i];
//     		// if (it<3) printf("%d %f\n",it,src.stf[i]);
//     	}

//     	Kokkos::deep_copy(srcd, srch);
//     	Kokkos::fence();

//     	for (int i=0; i<par.NS; ++i){
//     		if (par.Time < src.tlen[i]) {
//     			int ind = INDEX(src.isx[i],src.isy[i],par.NX,par.NY);

//     			Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     				if(index ==  ind) u2d(index)=srcd(i);
//     			});
//     			// printf("%d %d %d %d\n", src.isx[i],src.isy[i],ind,ind2);
//     		}	
//     	}
// 		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     		int i,j,i2;
// 			index2coord(index,i,j,par.NX,par.NY);

// 			// // Absorbing boundary on the top edge (0-x)
// 			if (BC1[0]==2 and i==0) {
// 				i2=INDEX(1,j,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}

// 			// // Absorbing boundary on the bottom edge (x-1)
// 			if (BC2[0]==2 and i==par.NX-1) {
// 				i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}


// 			// // Absorbing boundary on the left edge (0-y)
// 			if (BC1[1]==2 and j==0) {
// 				i2=INDEX(i,1,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}

// 			// // Absorbing boundary on the right edge (y-1)
// 			if (BC2[1]==2 and j==par.NY-1) {
// 				i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}


// 			// below are free surface BC
// 			if (BC1[0]==1 and i==0) {
// 				i2=INDEX(1,j,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC2[0]==1 and i==par.NX-1) {
// 				i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC1[1]==1 and j==0) {
// 				i2=INDEX(i,1,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC2[1]==1 and j==par.NY-1) {
// 				i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}
//     	});
// 		Kokkos::deep_copy(u0d,u1d);
//     	Kokkos::deep_copy(u1d,u2d);
// 		// timer.elapsed();
// 		//Kokkos::deep_copy(uh, u1d);
// 		// sprintf(sprintfBuffer, "out/Kokkos4_%03u.csv", fileIndex);
//       	// FILE* file = fopen(sprintfBuffer, "w");
//       	// if (!file){
//     	// 	std::cout << "Couldn't write the file!\nPlease 'mkdir out/' and rerun" << std::endl;
// 		// 	return 0;
// 		// }
// 		// int ind=0;
//       	// for (unsigned int i = 0; i < par.NX; i += 1) {
//       	// 	ind = INDEX(i,0,par.NX,par.NY);
// 		// 	fprintf(file, "%4.2f", uh(ind));
//         // 	for (unsigned int j = 1; j < par.NY; j += 1) {
//         // 		ind = INDEX(i,j,par.NX,par.NY);
//         // 		fprintf(file, ", %4.2f", uh(ind));
//         // 	}
//         // 	fprintf(file, "\n");
//       	// }

//       	// fclose(file);
//       	// ++fileIndex;
// 	}
// 	timer.stop();
// 	#ifdef DEBUG
// 		std::cout << "debug begin" << std::endl; 
//     	for (int index=0; index<par.NX*par.NY; ++index){
//     		if (abs(uh(index)-dataserial.u1[index])>1e-5) 
// 					printf("U1host= %f serial= %f\n",uh(index),dataserial.u1[index]);
// 		}
// 		std::cout << "debug end" << std::endl; 
//     #endif
// 	double time_seconds = timer.elapsed();
// 	return time_seconds;
// }

// real_t kokkos_5(Params par, DataContext dataserial)
// {
// 	char sprintfBuffer[500];
// 	int fileIndex = 0; 
// 	Timer timer;
// 	int BC1[2]={1,2}; //BC: 0 (rigid); 1 (free surface); 2 (absorbing);
//     int BC2[2]={2,2};
//     vel_module vm(par);
//     vm.readbin_vel((char*)"in/randomvel.bin");
// 	src_module src(par,1);
// 	DataArray1d c2h("c2h", par.NX*par.NY);
// 	DataArray1d c2d("c2d", par.NX*par.NY);

// 	DataArray1d ch("ch", par.NX*par.NY);
// 	DataArray1d cd("cd", par.NX*par.NY);
// 	DataArray1d srch("srch", par.NX*par.NY);
// 	DataArray1d srcd("srcd", par.NX*par.NY);

// 	DataArray1d uh("uh", par.NX*par.NY);

// 	DataArray2d u0d("u0d", par.NX, par.NY); 
// 	DataArray1d u1d("u1d", par.NX,par.NY);
// 	DataArray1d u2d("u2d", par.NX,par.NY);


// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		c2h[i] = vm.courant2[i];
//     }
//     Kokkos::deep_copy(c2d, c2h);

// 	for (int i=0; i<par.NX*par.NY; ++i){
// 		ch[i]= vm.courant[i];
//     }
// 	Kokkos::deep_copy(cd,ch);

// 	Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) {
// 		int i, j; 
// 		index2coord(index, i, j, par.NX, par.NY);
// 		u0d(i,j) = 0.;
// 		u1d(i,j) = 0.;
// 		u2d(i,j) = 0.;
//     });

//     Kokkos::fence();


//     for (int i=0; i<par.NX*par.NY; ++i){
// 		uh[i] = 0; 
// 	}
// 	fileIndex = 0;

//     for (int it=0; it < par.NT; ++it){
// 		timer.start();
//     	// computing time
//     	par.Time = it * par.dt;
// 		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) 
// 		{  
//     		int i,j;
// 			index2coord(index,i,j,par.NX,par.NY);

// 			const int ij = INDEX(i,j,par.NX,par.NY);
// 			const int ip1j = INDEX(i+1,j,par.NX,par.NY);
// 			const int im1j = INDEX(i-1,j,par.NX,par.NY);
// 			const int ijp1 = INDEX(i,j+1,par.NX,par.NY);
// 			const int ijm1 = INDEX(i,j-1,par.NX,par.NY);

// 			// if (index == 100) printf("%03d %03d %010d vs ",i,j,ij);
// 			if (i>0 and i<par.NX-1 and j>0 and j<par.NY-1 ){
// 				u2d(ij) = (2-4*c2d(ij))*u1d(ij)+
// 				c2d(ij)*
// 				(u1d(ip1j)+u1d(im1j)+
// 				 u1d(ijp1)+u1d(ijm1))-u0d(ij);
// 			}

//     	});
// 		src.add_src(par.Time);
//     	//deep copy n source from Host to Device
//     	for (int i=0; i<par.NS; ++i){
//     		srch(i) = src.stf[i];
//     		// if (it<3) printf("%d %f\n",it,src.stf[i]);
//     	}

//     	Kokkos::deep_copy(srcd, srch);
//     	Kokkos::fence();

//     	for (int i=0; i<par.NS; ++i){
//     		if (par.Time < src.tlen[i]) {
//     			int ind = INDEX(src.isx[i],src.isy[i],par.NX,par.NY);

//     			Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     				if(index ==  ind) u2d(index)=srcd(i);
//     			});
//     			// printf("%d %d %d %d\n", src.isx[i],src.isy[i],ind,ind2);
//     		}	
//     	}
// 		Kokkos::parallel_for( par.NX*par.NY, KOKKOS_LAMBDA(const int index) { 
//     		int i,j,i2;
// 			index2coord(index,i,j,par.NX,par.NY);

// 			// // Absorbing boundary on the top edge (0-x)
// 			if (BC1[0]==2 and i==0) {
// 				i2=INDEX(1,j,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}

// 			// // Absorbing boundary on the bottom edge (x-1)
// 			if (BC2[0]==2 and i==par.NX-1) {
// 				i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}


// 			// // Absorbing boundary on the left edge (0-y)
// 			if (BC1[1]==2 and j==0) {
// 				i2=INDEX(i,1,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}

// 			// // Absorbing boundary on the right edge (y-1)
// 			if (BC2[1]==2 and j==par.NY-1) {
// 				i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 				const real_t c = cd(i2);
// 				u2d(index) = c*u1d(i2)+
// 										(1-c)*u1d(index);
// 			}


// 			// below are free surface BC
// 			if (BC1[0]==1 and i==0) {
// 				i2=INDEX(1,j,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC2[0]==1 and i==par.NX-1) {
// 				i2=INDEX(par.NX-2,j,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC1[1]==1 and j==0) {
// 				i2=INDEX(i,1,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}


// 			if (BC2[1]==1 and j==par.NY-1) {
// 				i2=INDEX(i,par.NY-2,par.NX,par.NY);
// 				u2d(index) = u2d(i2);
// 			}
//     	});
// 		timer.elapsed();
// 		Kokkos::deep_copy(u0d,u1d);
//     	Kokkos::deep_copy(u1d,u2d);
// 		Kokkos::deep_copy(uh, u1d);
// 		sprintf(sprintfBuffer, "out/Kokkos_%03u.csv", fileIndex);
//       	FILE* file = fopen(sprintfBuffer, "w");
//       	if (!file){
//     		std::cout << "Couldn't write the file!\nPlease 'mkdir out/' and rerun" << std::endl;
// 			return 0;
// 		}
// 		int ind=0;
//       	for (unsigned int i = 0; i < par.NX; i += 1) {
//       		ind = INDEX(i,0,par.NX,par.NY);
// 			fprintf(file, "%4.2f", uh(ind));
//         	for (unsigned int j = 1; j < par.NY; j += 1) {
//         		ind = INDEX(i,j,par.NX,par.NY);
//         		fprintf(file, ", %4.2f", uh(ind));
//         	}
//         	fprintf(file, "\n");
//       	}

//       	fclose(file);
//       	++fileIndex;
// 	}
// 	timer.stop();
// 	#ifdef DEBUG
// 		std::cout << "debug begin" << std::endl; 
//     	for (int index=0; index<par.NX*par.NY; ++index){
//     		if (abs(uh(index)-dataserial.u1[index])>1e-5) 
// 					printf("U1host= %f serial= %f\n",uh(index),dataserial.u1[index]);
// 		}
// 		std::cout << "debug end" << std::endl; 
//     #endif
// 	double time_seconds = timer.elapsed();
// 	return time_seconds;
// }
