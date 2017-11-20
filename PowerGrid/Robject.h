/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Robject.h]

    Synopsis    [Base class of implementing quadratic and non-quadratic
                    regularization.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef PowerGrid_Robject_h
#define PowerGrid_Robject_h

#include "PGIncludes.h"

using namespace arma;
using namespace std;

template <typename T1> class Robject {
  typedef complex<T1> CxT1;

public:
  // Robject();
  Robject(){};
  // Class members
  uword Nx;
  uword Ny;
  uword Nz;
  T1 DeltaX;
  T1 DeltaY;
  T1 DeltaZ;
  T1 Beta;

  // It was declared as type Mat<uword> and the 3D type was a cube. We need to
  // vectorize it before it is passed to QuadPenalty.
  // Custom Class Constructor
  Robject(uword nx, uword ny, uword nz, T1 beta);

  // Class Methods - Declared virtual so they can be implemented in the base
  // classes. Also they are virtual so that if you try to call Robject, things
  // crash rather than give un results.
  virtual Col<CxT1> wpot(const Col<CxT1> &d) const {
    return ones<Col<CxT1>>(d.n_rows);
  }

  virtual Col<CxT1> dpot(const Col<CxT1> &d) const { return d; }

  virtual Col<CxT1> pot(const Col<CxT1> &d) const {
    Col<T1> temp = abs(d) % abs(d) / 2.0;
    return conv_to<Col<CxT1>>::from(temp);
  }

  Col<CxT1> Cd(const Col<CxT1> &d, uword dim) const;

  Col<CxT1> Ctd(const Col<CxT1> &d, uword dim) const;

  T1 Penalty(const Col<CxT1> &x) const;

  Col<CxT1> Gradient(const Col<CxT1> &x) const;

  CxT1 Denom(const Col<CxT1> &ddir, const Col<CxT1> &x) const;
};

// Explicit Instantiation
extern template class Robject<float>;
extern template class Robject<double>;

#endif // PowerGrid_Robject_h
