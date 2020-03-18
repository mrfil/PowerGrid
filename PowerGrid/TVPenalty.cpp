/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [TVPenalty.cpp]

    Synopsis    [Implementation of an approximate total variational penalty.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#include "TVPenalty.h"

// It was declared as type Mat<uword> and the 3D type was a cube. We need to
// vectorize it before it is passed to QuadPenalty.
// Custom Class Constructor
template <typename T1>
TVPenalty<T1>::TVPenalty(uword nx, uword ny, uword nz, T1 beta, T1 delta, uword dims2penalize) {
  // Set Class Members
  this->Nx = nx;
  this->Ny = ny;
  this->Nz = nz;
  this->Beta = beta;
  this->Delta = delta;
  this->Dims2Penalize = dims2penalize;
}

// Class Methods
template <typename T1>
Col<complex<T1>> TVPenalty<T1>::wpot(const Col<complex<T1>> &d) const {
  RANGE()
  Col<T1> temp = abs(d / this->Delta);
  Col<complex<T1>> out =
      1.0 / sqrt(1.0 + conv_to<Col<complex<T1>>>::from(temp % temp));
  return out;
}

template <typename T1>
Col<complex<T1>> TVPenalty<T1>::dpot(const Col<complex<T1>> &d) const {
  RANGE()
  Col<T1> temp = abs(d / this->Delta);
  Col<complex<T1>> out =
      d / sqrt(1.0 + conv_to<Col<complex<T1>>>::from(temp % temp));
  return out;
}

template <typename T1>
Col<complex<T1>> TVPenalty<T1>::pot(const Col<complex<T1>> &d) const {
  RANGE()
  Col<complex<T1>> out = conv_to<Col<complex<T1>>>::from(
      this->Delta * this->Delta *
      (sqrt(1.0 + (abs(d / this->Delta) % abs(d / this->Delta))) - 1.0));
  return out;
}

// Explicit Instantiation
template class TVPenalty<double>;
template class TVPenalty<float>;
