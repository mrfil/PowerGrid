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

#include "PGIncludes.h"

#include "fftw3.h"

#include <complex>
#include <iostream>
#include <type_traits>

using namespace std;
using namespace arma;
// Like Armadillo, we're using SFINAE here to choose between float and double.
// (Maybe FP16 some day in the future)
// We need enable_if to choose which version to run based on the type of the
// template parameter.

template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               uword>::type = 0>
void fft1dCPU(T1 *d_data, uword nx);

template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               uword>::type = 0>
void ifft2dCPU(T1 *d_data, uword nx, uword ny);

template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               uword>::type = 0>
void fft2dCPU(T1 *d_data, uword nx, uword ny);

template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               uword>::type = 0>
void ifft3dCPU(T1 *d_data, uword nx, uword ny, uword nz);

template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               uword>::type = 0>
void fft3dCPU(T1 *d_data, uword nx, uword ny, uword nz);

template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               uword>::type = 0>
void ifft2dCPU(T1 *d_data, uword nx, uword ny);

template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               uword>::type = 0>
void fft1dCPU(T1 *d_data, uword nx);


template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               uword>::type = 0>
void fft2dCPU(T1 *d_data, uword nx, uword ny);

template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               uword>::type = 0>
void ifft3dCPU(T1 *d_data, uword nx, uword ny, uword nz);

template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               uword>::type = 0>
void fft3dCPU(T1 *d_data, uword nx, uword ny, uword nz);

// Explicit Instantiations
extern template void fft1dCPU<float>(float *, uword);
extern template void fft1dCPU<double>(double *, uword);

extern template void ifft2dCPU<float>(float *, uword, uword);
extern template void ifft2dCPU<double>(double *, uword, uword);

extern template void fft2dCPU<float>(float *, uword, uword);
extern template void fft2dCPU<double>(double *, uword, uword);

extern template void ifft3dCPU<float>(float *, uword, uword, uword);
extern template void ifft3dCPU<double>(double *, uword, uword, uword);

extern template void fft3dCPU<float>(float *, uword, uword, uword);
extern template void fft3dCPU<double>(double *, uword, uword, uword);

#endif // PowerGrid_fftCPU_hpp
