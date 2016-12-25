/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [fftCPU.h]

    Synopsis    [Wrappers to the FFTW library supporting single and double
                                precision]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]
                [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/13/2016]

 *****************************************************************************/

#ifndef PowerGrid_fftCPU_hpp
#define PowerGrid_fftCPU_hpp

#include "fftw3.h"
#include <complex>
#include <type_traits>

// Like Armadillo, we're using SFINAE here to choose between float and double.
// (Maybe FP16 some day in the future)
// We need enable_if to choose which version to run based on the type of the
// template parameter.
template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void ifft2dCPU(T1 *d_data, int nx, int ny);

template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void fft2dCPU(T1 *d_data, int nx, int ny);

template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void ifft3dCPU(T1 *d_data, int nx, int ny, int nz);

template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void fft3dCPU(T1 *d_data, int nx, int ny, int nz);

template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void ifft2dCPU(T1 *d_data, int nx, int ny);

template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void fft2dCPU(T1 *d_data, int nx, int ny);

template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void ifft3dCPU(T1 *d_data, int nx, int ny, int nz);

template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void fft3dCPU(T1 *d_data, int nx, int ny, int nz);

// Explicit Instantiations
template void ifft2dCPU<float>(float *, int, int);
template void ifft2dCPU<double>(double *, int, int);

template void fft2dCPU<float>(float *, int, int);
template void fft2dCPU<double>(double *, int, int);

template void ifft3dCPU<float>(float *, int, int, int);
template void ifft3dCPU<double>(double *, int, int, int);

template void fft3dCPU<float>(float *, int, int, int);
template void fft3dCPU<double>(double *, int, int, int);

#endif // PowerGrid_fftCPU_hpp
