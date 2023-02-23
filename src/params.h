#ifndef PARAMS_H_
#define PARAMS_H_

struct Params {

  int NX;
  int NY;
  int NT;
  int NS;
  real_t dt;
  real_t dx;
  real_t Time=0.;
  
  Params() :
    NX(1024), NY(1024),dt(0.1),dx(100),NS(1)
  {}
  
 Params(int NX_, int NY_, int NS_,real_t dt_, real_t dx_) :
    NX(NX_), NY(NY_), NS(NS_), dt(dt_), dx(dx_)
  {}

  Params(int NX_, int NY_, int NS_,real_t dt_, real_t dx_, int nt_) :
    NX(NX_), NY(NY_), NS(NS_), dt(dt_), dx(dx_), NT(nt_)
  {}
  
}; // struct Params

#endif