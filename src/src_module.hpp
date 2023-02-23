#ifndef SRC_MODULE_H_
#define SRC_MODULE_H_

#include <cmath>
#include "real_type.h"
#include "params.h"

class src_module{
public:

	int ns; 	//total number of sources
	int *isx;	//source location ix array
	int *isy;	//source location iy array
	real_t *stf; /* array of source pressure values */
	real_t *alpha; //source bandwidth
	real_t *tlen; // source duration in sec
	real_t *ts; //time shift in sec

	src_module(Params &par,int n_):
	ns(n_){
		isx = (int *) malloc(sizeof(int) * ns);
		isy = (int *) malloc(sizeof(int) * ns);
		stf = (real_t *) malloc(sizeof(real_t) * ns);
		alpha = (real_t *) malloc(sizeof(real_t) * ns);
		tlen = (real_t *) malloc(sizeof(real_t) * ns);
		ts = (real_t *) malloc(sizeof(real_t) * ns);


		for (int isrc=0; isrc<ns; ++isrc){
			isx[isrc] = int(par.NX/2.);
			isy[isrc] = int(par.NY/2.);
			stf[isrc] = 0;
			alpha[isrc] = 800;
			tlen[isrc] = 0.45;
			ts[isrc] = 0.2;
		}
		#ifdef DEBUG
		printf("[src_module] #src= %d   alpha=%f tlen=%5.2f(s) ts=%5.2f(s)\n",
			ns,alpha[0],tlen[0],ts[0]);
		#endif
	}

	virtual ~src_module(){
		delete[] isx;
		delete[] isy;
		delete[] stf;
		delete[] alpha;
		delete[] tlen;
		delete[] ts;
	}

	real_t src_func (real_t ,real_t ,real_t );
	void add_src(real_t);
};


real_t src_module::src_func (real_t t_,real_t shift_,real_t alpha_){
	return 3*exp(-alpha_*(t_-shift_)*(t_-shift_));
}


void src_module::add_src(real_t t){
	for (int isrc=0; isrc<ns; ++isrc){
		if (t < tlen[isrc]){
			stf[isrc] = src_func(t,ts[isrc],alpha[isrc]);
		}else{
			stf[isrc]=0.;
		}
	}
}


#endif