/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gfft.h]

    Synopsis    [Object that represents a uniform discrete Fourier
                    transform implemented via a fast Fourier transform (FFT).]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef __PowerGrid__Gfft__h
#define __PowerGrid__Gfft__h

#ifdef _OPENACC // GPU Version

#include "cufft.h"
#include "fftGPU.h"
#include "gridding.h"
#include "openacc.h"
#else // CPU Version

#include "fftCPU.h"
#include "gridding.h"
#endif

// We want to operate on many types of variables (assume of type Col<type>)
template <typename T1> class Gfft {
  typedef complex<T1> CxT1;

public:
  // Default Class Constructor and Destructor
  Gfft();
  // Class Constructor
  Gfft(uword ix, uword iy, uword iz);

  // Class variables go here.
  uword Nx = 0; // Size in x dimension
  uword Ny = 0; // Size in y dimension
  uword Nz = 0; // Size in z dimension

  // Overloaded methods for forward and adjoint transform
  // Forward transform operation
  Col<CxT1> operator*(const Col<CxT1> &d) const;

  // Adjoint transform operation
  Col<CxT1> operator/(const Col<CxT1> &d) const;
};

// Explicit Instantiations
extern template class Gfft<float>;
extern template class Gfft<double>;

#endif // __PowerGrid__Gfft__h
