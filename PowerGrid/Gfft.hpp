/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gfft.hpp]

    Synopsis    [Object that represents a uniform discrete Fourier
                    transform implemented via a fast Fourier transform (FFT).]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef __PowerGrid__Gfft__hpp
#define __PowerGrid__Gfft__hpp

#ifdef _OPENACC // GPU Version

#include "cufft.h"
#include "fftGPU.hpp"
#include "openacc.h"

#else // CPU Version

#include "fftCPU.hpp"

#endif

// We want to operate on many types of variables (assume of type Col<type>)
template <typename T1> class Gfft {
  typedef complex<T1> CxT1;

public:
  // Default Class Constructor and Destructor
  Gfft();
  // Class Constructor
  Gfft(uword ix, uword iy, uword iz) {
    Nx = ix;
    Ny = iy;
    Nz = iz;
  }

  // Class variables go here.
  uword Nx = 0; // Size in x dimension
  uword Ny = 0; // Size in y dimension
  uword Nz = 0; // Size in z dimension

  // Overloaded methods for forward and adjoint transform
  // Forward transform operation
  Col<CxT1> operator*(const Col<CxT1> &d) const {
    uword stx = (this->n2 - this->n1) / 2;

    Col<T1> realData = real(d);
    Col<T1> imagData = imag(d);

    T1 *realDataPtr = realData.memptr();
    T1 *imagDataPtr = imagData.memptr();

    int imageNumElems = Nx * Ny * Nz;

    Col<T1> realXformedData;
    Col<T1> imagXformedData;
    realXformedData.zeros(imageNumElems);
    imagXformedData.zeros(imageNumElems);

    T1 *realXformedDataPtr = realXformedData.memptr();
    T1 *imagXformedDataPtr = imagXformedData.memptr();
    // allocate gridData
    CxT1 *gridData = new CxT1[imageNumElems];
    CxT1 *gridData_d = new CxT1[imageNumElems];

    // Have to set 'gridData' to zero.
    // Because they will be involved in accumulative operations
    // inside gridding functions.
    for (int i = 0; i < imageNumElems; i++) {
      gridData[i].real(realDataPtr[i]);
      gridData[i].imag(imagDataPtr[i]);
    }

    T1 *pGridData = reinterpret_cast<T1 *>(gridData);
    T1 *pGridData_d = reinterpret_cast<T1 *>(gridData_d);

// fftn(gridData)
#ifdef _OPENACC // We're on GPU
                // Inside this region the device data pointer will be used
// cout << "about to reach openacc region in forward transform" << endl;

#pragma acc data copyin(pGridData[0 : 2 * imageNumElems]) create(              \
    pGridData_d[0 : 2 * imageNumElems])                                        \
        copyout(realXformedData[0 : imageNumElems],                            \
                                imagXformedData[0 : imageNumElems])
    {

#pragma acc host_data use_device(pGridData_d, pGridData)
      {
        // Query OpenACC for CUDA stream
        void *stream = acc_get_cuda_stream(acc_async_sync);

        // Launch FFT on the GPU
        if (Nz == 1) {
          fftshift2<T1>(pGridData_d, pGridData, Nx, Ny);
          fft2dGPU(pGridData_d, Nx, Ny, stream);
          fftshift2<T1>(pGridData, pGridData_d, Nx, Ny);
          deinterleave_data2d<T1>(pGridData, realXformedData, imagXformedData,
                                  Nx, Ny);
        } else {
          fftshift3<T1>(pGridData_d, pGridData, Nx, Ny, Nz);
          fft3dGPU(pGridData_d, Nx, Ny, Nz, stream);
          fftshift3<T1>(pGridData, pGridData_d, Nx, Ny, Nz);
          deinterleave_data3d<T1>(pGridData, realXformedData, imagXformedData,
                                  Nx, Ny, Nz);
        }
      }
    }

#else // We're on CPU
    if (Nz == 1) {
      fftshift2<T1>(pGridData_d, pGridData, Nx, Ny);
      fft2dCPU(pGridData_d, Nx, Ny);
      fftshift2<T1>(pGridData, pGridData_d, Nx, Ny);
      deinterleave_data2d<T1>(pGridData, realXformedData, imagXformedData, Nx,
                              Ny);
    } else {
      fftshift3<T1>(pGridData_d, pGridData, Nx, Ny, Nz);
      fft3dCPU(pGridData_d, Nx, Ny, Nz);
      fftshift3<T1>(pGridData, pGridData_d, Nx, Ny, Nz);
      deinterleave_data3d<T1>(pGridData, realXformedData, imagXformedData, Nx,
                              Ny, Nz);
    }
#endif

    Col<CxT1> XformedData(Nx * Ny * Nz);
    XformedData.set_real(realXformedData);
    XformedData.set_imag(imagXformedData);

    return conv_to<Col<CxT1>>::from(XformedData); // Return a vector of type T1
  }

  // Adjoint transform operation
  Col<CxT1> operator/(const Col<CxT1> &d) const {
    Col<T1> realData = real(d);
    Col<T1> imagData = imag(d);

    T1 *realDataPtr = realData.memptr();
    T1 *imagDataPtr = imagData.memptr();

    int imageNumElems = Nx * Ny * Nz;

    Col<T1> realXformedData;
    Col<T1> imagXformedData;
    realXformedData.zeros(imageNumElems);
    imagXformedData.zeros(imageNumElems);

    T1 *realXformedDataPtr = realXformedData.memptr();
    T1 *imagXformedDataPtr = imagXformedData.memptr();

    // allocate gridData
    CxT1 *gridData = new CxT1[imageNumElems];
    CxT1 *gridData_d = new CxT1[imageNumElems];

    // Have to set 'gridData' to zero.
    // Because they will be involved in accumulative operations
    // inside gridding functions.
    for (int i = 0; i < imageNumElems; i++) {
      gridData[i].real(realDataPtr[i]);
      gridData[i].imag(imagDataPtr[i]);
    }

    T1 *pGridData = reinterpret_cast<T1 *>(gridData);
    T1 *pGridData_d = reinterpret_cast<T1 *>(gridData_d);
// fftn(gridData)

#ifdef _OPENACC // We're on GPU
                // Inside this region the device data pointer will be used
// cout << "about to reach openacc region in forward transform" << endl;

#pragma acc data copyin(pGridData[0 : 2 * imageNumElems]) create(              \
    pGridData_d[0 : 2 * imageNumElems])                                        \
        copyout(realXformedData[0 : imageNumElems],                            \
                                imagXformedData[0 : imageNumElems])
    {

#pragma acc host_data use_device(pGridData_d, pGridData)
      {
        // Query OpenACC for CUDA stream
        void *stream = acc_get_cuda_stream(acc_async_sync);

        // Launch FFT on the GPU
        if (Nz == 1) {
          ifftshift2<T1>(pGridData_d, pGridData, Nx, Ny);
          ifft2dGPU(pGridData_d, Nx, Ny, stream);
          ifftshift2<T1>(pGridData, pGridData_d, Nx, Ny);
          deinterleave_data2d<T1>(pGridData, realXformedData, imagXformedData,
                                  Nx, Ny);
        } else {
          ifftshift3<T1>(pGridData_d, pGridData, Nx, Ny, Nz);
          ifft3dGPU(pGridData_d, Nx, Ny, Nz, stream);
          ifftshift3<T1>(pGridData, pGridData_d, Nx, Ny, Nz);
          deinterleave_data3d<T1>(pGridData, realXformedData, imagXformedData,
                                  Nx, Ny, Nz);
        }
      }
    }

#else // We're on CPU
    if (Nz == 1) {
      ifftshift2<T1>(pGridData_d, pGridData, Nx, Ny);
      ifft2dCPU(pGridData_d, Nx, Ny);
      ifftshift2<T1>(pGridData, pGridData_d, Nx, Ny);
      deinterleave_data2d<T1>(pGridData, realXformedData, imagXformedData, Nx,
                              Ny);
    } else {
      ifftshift3<T1>(pGridData_d, pGridData, Nx, Ny, Nz);
      ifft3dCPU(pGridData_d, Nx, Ny, Nz);
      ifftshift3<T1>(pGridData, pGridData_d, Nx, Ny, Nz);
      deinterleave_data3d<T1>(pGridData, realXformedData, imagXformedData, Nx,
                              Ny, Nz);
    }
#endif
  }
};
#endif
