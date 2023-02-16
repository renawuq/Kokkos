# Acoustic_kokkos

## Download this app

```bash
git clone --recursive https://github.com/wangyf/acoustic_kokkos
```

## How to build ?

### Requirement

- [cmake](https://cmake.org/) version 3.16

- **note**: if you are on a fairly recent OS (ex: Ubuntu 21.10, or any OS using glibc >= 2.34), you may need to turn off linking with libdl when using kokkos/cuda backend. See [this issue](https://github.com/kokkos/kokkos/issues/4824), as nvcc (even version 11.6) apparently doesn't seem to handle empty file `/usr/lib/x86_64-linux-gnu/libdl.a` (stub, libdl is integrated into glibc). Hopefully this will be solved in an upcoming cuda release.

```shell
# run this to know your glibc version
ldd --version
```

### Data Preparation
> 1. In build_xxx, create directories: in/ and out/
> 2. Copy the input velocity structure (heterogeneous) into the in/
> 3. All output .csv will be stored in out/

### Build with target device OpenMP

```bash
mkdir build_openmp
cd build_openmp
CXX=YOUR_COMPILER_HERE cmake -DKokkos_ENABLE_OPENMP=ON ..
make
# then you can run the application
./src/acoustic_kokkos.openmp
```

Optionnally you can enable HWLOC by passing -DKokkos_ENABLE_HWLOC=ON on cmake's command line (or in ccmake curse gui).

### Build with target device CUDA

CMake and Kokkos will set the compiler to `nvcc_wrapper` (located in kokkos sources, cloned as git submodule).

```bash
mkdir build_cuda
cd build_cuda
cmake -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -DKokkos_ARCH_MAXWELL50=ON ..
make
# then you can run the application as before
./src/acoustic_kokkos.cuda
```

Of course, you will need to adapt variable **Kokkos_ARCH** to your actual GPU architecture (use cuda sample device_query to probe the architecture).

Depending on your OS, you may need to set variable **Kokkos_CUDA_DIR** to point to your CUDA SDK (if cmake is not able to figure out by itself); e.g. /usr/local/cuda-9.0

### Build with target device HIP (AMD GPU)

CMake and Kokkos will set the compiler to `hipcc` (located in kokkos sources, cloned as git submodule).

Example:
```bash
mkdir build_hip
cd build_hip
cmake -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_VEGA908=ON ..
make
# then you can run the application as before
./src/acoustic_kokkos.hip
```


### Example complish options:
Please read compile.sh


### Animation example
![Example](/ani.gif)

# Kokkos
