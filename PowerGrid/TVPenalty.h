/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [TVPenalty.h]

    Synopsis    [Implementation of an approximate total variational penalty.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef PowerGrid_TVPenalty_h
#define PowerGrid_TVPenalty_h

#include "Robject.h"

using namespace arma;

template <typename T1> class TVPenalty : public Robject<T1> {
  typedef complex<T1> CxT1;

public:
  TVPenalty();

  // It was declared as type Mat<uword> and the 3D type was a cube. We need to
  // vectorize it before it is passed to QuadPenalty.
  // Custom Class Constructor
  TVPenalty(uword nx, uword ny, uword nz, T1 beta, T1 delta, uword dims2penalize = 3);

  // Class Methods

  Col<CxT1> wpot(const Col<CxT1> &d) const;

  virtual Col<CxT1> dpot(const Col<CxT1> &d) const;

  Col<CxT1> pot(const Col<CxT1> &d) const;

private:
  T1 Delta;
};

// Explicit Instantiation
extern template class TVPenalty<double>;
extern template class TVPenalty<float>;

#endif // PowerGrid_TVPenalty_h
