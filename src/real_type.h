#ifndef REAL_TYPE_H_
#define REAL_TYPE_H_


#ifdef USE_DOUBLE
using real_t = double;
#else
using real_t = float;
#endif // USE_DOUBLE


# define PI 3.14159265358979323846  /* pi */


#endif // REAL_TYPE_H_