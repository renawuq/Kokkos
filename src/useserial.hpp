#ifndef USESERIAL_H_
#define USESERIAL_H_

#include <iostream>
#include "real_type.h"

struct DataContext {

  real_t * u0;
  real_t * u1;
  real_t * u2;

  DataContext(Params& params)
  { u0 = new real_t[params.NX*params.NY];
    u1 = new real_t[params.NX*params.NY];
    u2 = new real_t[params.NX*params.NY];
  }

  ~DataContext() {}
  
};


#endif