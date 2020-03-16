/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gnufft.cpp]

    Synopsis    [Object that represents a non-uniform discrete Fourier
                    transform implemented via a a CPU and GPU accelerated
                    gridding implementation of the non-uniform Fast Fourier
                    Tranform (NUFFT).]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#include "Gnufft.h"
#include <chrono>
template <typename T1>
Gnufft<T1>::Gnufft(
    uword dataLength, T1 gridos, uword nx, uword ny, uword nz,
    const Col<T1> &k1, const Col<T1> &k2, const Col<T1> &k3, const Col<T1> &i1,
    const Col<T1> &i2,
    const Col<T1> &i3) // Change these arguments as you need to setup the object
{
  
  // cout << "Entering the class constructor for Ggrid" << endl;
  n1 = nx * ny * nz;
  n2 = dataLength;
  Nx = nx;
  Ny = ny;
  Nz = nz;
  ix = i1;
  iy = i2;
  iz = i3;
  //kx = k1;
  //ky = k2;
  //kz = k3;
  gridOS = gridos;

  imageNumElems = Nx * Ny;
  gridNumElems = gridOS * Nx * gridOS * Ny;

  if (Nz > 1) {
    imageNumElems = imageNumElems * Nz;
    gridNumElems = gridNumElems * gridOS * Nz;
  }

  kx = new T1[n2];
  ky = new T1[n2];
  kz = new T1[n2];

  for (int i = 0; i < n2; i++) {
    if (abs(k1(i)) > ((Nx / (T1)2.0) * 1.05) || abs(k2(i)) > ((Ny / (T1)2.0) * 1.05) ||
        abs(k3(i)) > ((Nz / (T1)2.0) * 1.05)) {
      
      cout << "Warning :k-space trajectory out of range [-N/2,N/2]:\n      "
            << "gridding requires that k-space should be contained within the "
            << "window -N/2 to N/2" << endl;
      cout << "kx = " << k1(i) << " ky = " << k2(i)<< " kz = " << k3(i)
           << " i = " << i << endl;
      }
      
      kx[i] = k1(i);
      ky[i] = k2(i);
      kz[i] = k3(i);
  }

  // Set Beta
  kernelWidth = 4.0;
  beta = MRI_PI *
         std::sqrt((gridOS - 0.5) * (gridOS - 0.5) *
                       (kernelWidth * kernelWidth * 4.0) / (gridOS * gridOS) -
                   0.8);
  #ifdef OPENACC_GPU
  stream = acc_get_cuda_stream(acc_async_sync);
  cufftCreate(&plan);
  if(Nz ==1) {
    if (cufftPlan2d(&plan, gridOS*Nx, gridOS*Ny, CUFFT_C2C) != CUFFT_SUCCESS) {
            cout <<  "CUFFT error: Plan creation failed" << endl;
    }
  } else {
    if (cufftPlan3d(&plan, gridOS*Nz, gridOS*Ny, gridOS*Nx, CUFFT_C2C) != CUFFT_SUCCESS) {
            cout << "CUFFT error: Plan creation failed" << endl;
          }
  }
  cufftSetStream(plan, (cudaStream_t)stream);
  #endif

  gridData = new complex<T1>[imageNumElems];
	gridData_d = new complex<T1>[imageNumElems];
  gridData_os = new complex<T1>[gridNumElems];
  gridData_os_d = new complex<T1>[gridNumElems];
  samples = new complex<T1>[n2];
  pGridData = reinterpret_cast<T1 *>(gridData);
  pGridData_d = reinterpret_cast<T1 *>(gridData_d);
  pGridData_os = reinterpret_cast<T1 *>(gridData_os);
  pGridData_os_d = reinterpret_cast<T1 *>(gridData_os_d);
  pSamples = reinterpret_cast<T1 *>(samples);
  // Deal with the LUT
  // Generating Look-Up Table
  // cout << "Calculating look up table" << endl;
  calculateLUT(beta, kernelWidth, LUT, sizeLUT);

#pragma acc enter data copyin(LUT[0 : sizeLUT], kx[0:n2], ky[0:n2], \
  kz[0:n2]) create(pGridData[0:2*imageNumElems], pGridData_d[0:2*imageNumElems], \
  pGridData_os[0:2*gridNumElems], pGridData_os_d[0:2*gridNumElems], pSamples[0:2*n2])
}

// Class destructor to free LUT
template <typename T1> Gnufft<T1>::~Gnufft() {
  RANGE()
  #ifdef OPENACC_GPU
    cufftDestroy(plan);
  #endif
  #pragma acc exit data delete(pGridData[0:2*imageNumElems], \
   pGridData_d[0:2*imageNumElems], pGridData_os[0:2*gridNumElems], \
   pGridData_os_d[0:2*gridNumElems], kx[0:n2], ky[0:n2], kz[0:n2], \
   pSamples[0:2*n2])

  delete[] samples;
  delete[] gridData;
  delete[] gridData_d;
  delete[] gridData_os;
  delete[] gridData_os_d;
  delete[] kx;
  delete[] ky;
  delete[] kz;

  if (LUT) {
#pragma acc exit data delete (LUT)
    free(LUT);
  }
}

// Overloaded methods for forward and adjoint transform
// Forward transform operation using gridding
template <typename T1>
inline Col<complex<T1>> Gnufft<T1>::
operator*(const Col<complex<T1>> &d) const // Don't change these arguments
{
RANGE()

  const T1 *dataPtr = reinterpret_cast<const T1 *>(d.memptr());

  // T2 gridOS = 2.0;
  // cout << "About to call the forward gridding routine." << endl;
  #ifdef OPENACC_GPU
    cufftHandle *nPlan = const_cast<cufftHandle *>(&plan);
  #else
    void* nPlan = NULL;
  #endif
  computeFd_CPU_Grid<T1>(n2, kx, ky, kz, dataPtr,
                         Nx, Ny, Nz, gridOS,
                         kernelWidth, beta, LUT, sizeLUT,
                         stream, nPlan, pGridData, pGridData_d, pGridData_os,
                         pGridData_os_d, pSamples);

  Col<CxT1> temp(reinterpret_cast<CxT1 *>(pSamples), n2, false, true);
  return temp; // Return a vector of type T1
}

// Adjoint transform operation
template <typename T1>
inline Col<complex<T1>> Gnufft<T1>::operator/(const Col<complex<T1>> &d) const {
  // uword dataLength = n2;
  // Let's trim the operations to avoid data overhead and transfers
  // Basically if we know that the data points are zero, they have no impact
  // on the transform

  uword dataLength = this->n2;

  //Col<T1> realData = real(d).eval();
  //Col<T1> imagData = imag(d).eval();

  const T1 *dataPtr = reinterpret_cast<const T1 *>(d.memptr());

  // Process data here, like calling a brute force transform, dft...
  // I assume you create the pointers to the arrays where the transformed data
  // will be stored
  // realXformedDataPtr and imagXformedDataPtr and they are of type float*

  // T2 gridOS = 2.0;
  #ifdef OPENACC_GPU
  cufftHandle *nPlan = const_cast<cufftHandle *>(&plan);
  #else
  void *nPlan = NULL;
  #endif
  computeFH_CPU_Grid<T1>(dataLength, kx, ky, kz,
                         dataPtr, Nx, Ny, Nz, gridOS,
                         kernelWidth,
                         beta, LUT, sizeLUT, stream, nPlan, pGridData,
                         pGridData_d, pGridData_os, pGridData_os_d);
  
  Col<CxT1> temp(reinterpret_cast<CxT1 *>(pGridData), n1, false, true);
  return temp; // Return a vector of type T1
}

template <typename T1>
inline Col<complex<T1>> Gnufft<T1>::forwardSpatialInterp(const Col<complex<T1>> &d) const {
  uword dataLength = this->n2;

  const T1 *dataPtr = reinterpret_cast<const T1 *>(d.memptr());

  parameters<T1> params;
  params.sync = 0;
  params.binsize = 128;
  params.useLUT = 1;
  params.kernelWidth = kernelWidth;
  params.gridOS = gridOS;
  params.imageSize[0] = Nx; // gridSize is gridOS times larger than imageSize.
  params.imageSize[1] = Ny;
  params.imageSize[2] = Nz;
  params.gridSize[0] = std::ceil(gridOS * (T1)Nx);
  params.gridSize[1] = std::ceil(gridOS * (T1)Ny);
  if (params.gridSize[0] % 2) // 3D case, gridOS is adjusted on the z dimension:
    params.gridSize[0] += 1; // That why we need to make sure here that the xy
  if (params.gridSize[1] % 2) // dimensions have even sizes.
    params.gridSize[1] += 1;
  params.gridSize[2] = (Nz == 1) ? Nz : (std::ceil(gridOS * (T1)Nz)); // 2D or 3D
  params.numSamples = dataLength;

  unsigned int n = params.numSamples;
  
  memcpy(pGridData_os, dataPtr, sizeof(T1) * 2 * gridNumElems);
  #pragma acc update device(pGridData_os[0:2 * gridNumElems])

  //#pragma acc parallel loop present(pSamples [0:2 * n])
    for (int ii = 0; ii < 2 * n; ii++) {
        pSamples[ii] = (T1)0.0;
    }
    #pragma acc update device(pSamples [0:2 * n])

    if (Nz == 1) {
      gridding_forward_2D(n, params, kx,
                         ky, beta, pSamples, LUT,
                          sizeLUT, pGridData_os);
    } else {
      gridding_forward_3D(n, params, kx,
                      ky, kz, beta, pSamples, LUT,
                      sizeLUT, pGridData_os);
    }

  #pragma acc update host(pSamples [0:2 * n])

  Col<CxT1> temp(reinterpret_cast<CxT1 *>(pSamples), n, false, true);
  return temp; // Return a vector of type T1

}

template <typename T1>
inline Col<complex<T1>> Gnufft<T1>::adjointSpatialInterp(const Col<complex<T1>> &d) const {

  uword dataLength = this->n2;

  const T1 *dataPtr = reinterpret_cast<const T1 *>(d.memptr());

  parameters<T1> params;
  params.sync = 0;
  params.binsize = 128;
  params.useLUT = 1;
  params.kernelWidth = kernelWidth;
  params.gridOS = gridOS;
  params.imageSize[0] = Nx; // gridSize is gridOS times larger than imageSize.
  params.imageSize[1] = Ny;
  params.imageSize[2] = Nz;
  params.gridSize[0] = std::ceil(gridOS * (T1)Nx);
  params.gridSize[1] = std::ceil(gridOS * (T1)Ny);
  if (params.gridSize[0] % 2) // 3D case, gridOS is adjusted on the z dimension:
    params.gridSize[0] += 1; // That why we need to make sure here that the xy
  if (params.gridSize[1] % 2) // dimensions have even sizes.
    params.gridSize[1] += 1;
  params.gridSize[2] = (Nz == 1) ? Nz : (std::ceil(gridOS * (T1)Nz)); // 2D or 3D
  params.numSamples = dataLength;

  unsigned int n = params.numSamples;

  ReconstructionSample<T1>* samples; // Input Data
  // allocate samples
  samples = (ReconstructionSample<T1>*)malloc(
  params.numSamples * sizeof(ReconstructionSample<T1>));

  if (samples == NULL) {
    printf("ERROR: Unable to allocate memory for input data\n");
    exit(1);
  }

  //
  for (int i = 0; i < params.numSamples; i++) {

    samples[i].kX = kx[i];
    samples[i].kY = ky[i];
    samples[i].kZ = kz[i];

    samples[i].real = dataPtr[2 * i];
    samples[i].imag = dataPtr[2 * i + 1];

    samples[i].sdc = (T1)1.0;
    // samples[i].t = t[i];
  }

  #pragma acc enter data copyin(samples [0:n])

  //#pragma acc parallel loop 
    for (int i = 0; i < gridNumElems; i++) {
        pGridData_os[2 * i] = (T1)0.0;
        pGridData_os[2 * i + 1] = (T1)0.0;
    }
  
  #pragma acc update device(pGridData_os[0:2*gridNumElems])

  if (Nz == 1) {
    gridding_adjoint_2D(n, params, beta, samples,
                        LUT, sizeLUT, pGridData_os);
  } else {
    gridding_adjoint_3D(n, params, beta, samples,
                        LUT, sizeLUT, pGridData_os);
  }

  #pragma acc update host(pGridData_os[0:2*gridNumElems])
  Col<CxT1> temp(reinterpret_cast<CxT1 *>(pGridData_os), gridNumElems, false, true);
  return temp; // Return a vector of type T1
}

// Explicit Instantiation
template class Gnufft<float>;
template class Gnufft<double>;
