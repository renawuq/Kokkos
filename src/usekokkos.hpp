#ifndef USEKOKKOS_H_
#define USEKOKKOS_H_

#include <iostream>
// Include Kokkos Headers
#include<Kokkos_Core.hpp>
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

KOKKOS_INLINE_FUNCTION
int INDEX(int i,  int j,
          int Nx, int Ny)
{
  UNUSED(Nx);
#ifdef KOKKOS_ENABLE_CUDA
  return i + Nx*j; // left layout
#else
  return j + Ny*i; // right layout
#endif
} // INDEX



KOKKOS_INLINE_FUNCTION
void index2coord(int index,
                 int &i, int &j,
                 int Nx, int Ny)
{
  UNUSED(Nx);
#ifdef KOKKOS_ENABLE_CUDA
  j = index / Nx;
  i = index - j*Nx;
#else
  i = index / Ny;
  j = index - i*Ny;
#endif
} // index2coord - 2d


// define Kokkos execution space
using Device = Kokkos::DefaultExecutionSpace;
// Data array for laplace computation
typedef Kokkos::View<real_t*, Device> DataArray;
typedef Kokkos::View<int*, Device> DataInt;

// host mirror
typedef DataArray::HostMirror  DataArrayHost;
typedef DataInt::HostMirror  DataIntHost;


void kokkos_config(){
  #ifdef DEBUG
  std::cout << "##########################\n";
  std::cout << "KOKKOS CONFIG             \n";
  std::cout << "##########################\n";
  
  std::ostringstream msg;
  std::cout << "Kokkos configuration" << std::endl;
  if ( Kokkos::hwloc::available() ) {
    msg << "hwloc( NUMA[" << Kokkos::hwloc::get_available_numa_count()
  << "] x CORE["    << Kokkos::hwloc::get_available_cores_per_numa()
  << "] x HT["      << Kokkos::hwloc::get_available_threads_per_core()
  << "] )"
  << std::endl ;
  }
  Kokkos::print_configuration( msg );
  std::cout << msg.str();
  std::cout << "##########################\n";
  #endif
}


struct DataContextKokkos {
  DataArray c2_d;
  DataArray c_d;
  DataArray u0_d;
  DataArray u1_d;
  DataArray u2_d;
  DataArray src_d;

  DataArrayHost c_h;
  DataArrayHost c2_h;
  DataArrayHost u_h;
  DataArrayHost src_h;


  DataContextKokkos(Params& params) :
  	// DataContextKokkos datakokkos(params);
    c2_d      ("c2_d", params.NX*params.NY),
    c_d       ("c_d", params.NX*params.NY),
    u0_d      ("u0_d", params.NX*params.NY),
    u1_d      ("u1_d", params.NX*params.NY),
    u2_d      ("u2_d", params.NX*params.NY),
    src_d     ("src_d",params.NS)

  { c2_h=Kokkos::create_mirror_view(c2_d);
    c_h=Kokkos::create_mirror_view(c_d);
    u_h=Kokkos::create_mirror_view(u1_d);
    src_h=Kokkos::create_mirror_view(src_d);
  }

  ~DataContextKokkos() {}
  
};

#endif