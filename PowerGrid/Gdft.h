/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gdft.h]

    Synopsis    [Object that represents a non-uniform field corrected discrete
                    Fourier tranform.]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/2/2016]

 *****************************************************************************/

#ifndef PowerGrid_Gdft_h
#define PowerGrid_Gdft_h

#include "PGIncludes.h"
#include "ftCpu.h"

using namespace arma;
using namespace std;

template <typename T1> // This is of type complex<double> or complex<float>, or
// any other type like float or single
class Gdft {
  typedef complex<T1> CxT1;

public:
  // Default Class Constructor and Destructor
  Gdft();
  // Class Constructor
  Gdft(uword a, uword b, const Col<T1> &k1, const Col<T1> &k2,
       const Col<T1> &k3, const Col<T1> &i1, const Col<T1> &i2,
       const Col<T1> &i3, const Col<T1> &f1, const Col<T1> &t1);

  // Class variables go here. Change as necessary
  uword n1 = 0;
  uword n2 = 0;

  Col<T1> kx; // k-space coordinates
  Col<T1> ky;
  Col<T1> kz;
  Col<T1> ix; // image space coordinates
  Col<T1> iy;
  Col<T1> iz;
  Col<T1> FM;
  Col<T1> t;

  // Overloaded methods for forward and adjoint transform
  // Forward transform operation
  Col<CxT1> operator*(const Col<CxT1> &d) const;
  // Adjoint transform operation
  Col<CxT1> operator/(const Col<CxT1> &d) const;
};

extern template class Gdft<float>;
extern template class Gdft<double>;

#endif // PowerGrid_Gdft_h
