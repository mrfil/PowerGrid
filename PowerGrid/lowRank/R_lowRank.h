/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [R_lowRank.h]

    Synopsis    [Wrapper class to extend penalization in the low rank recon case.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [03/16/2020]

 *****************************************************************************/


#ifndef PowerGrid_R_lowRank_h
#define PowerGrid_R_lowRank_h

#include "PGIncludes.h"
#include "QuadPenalty.h"
#include "TVPenalty.h"

using namespace arma;
using namespace std;

template <typename T1, typename Robj> class R_lowRank {
  typedef complex<T1> CxT1;
  
public:
  // R_lowRank();
  R_lowRank(){};
  // Class members
  uword Nrank; // Rank of the problem
  const Robj* R;
  uword Nimg; //
  Mat<CxT1> v; //object storing temporal basis

  R_lowRank(Robj const &R1, uword rank, uword nImages, Mat<CxT1> vBasis);

  T1 Penalty(const Col<CxT1> &x) const;

  Col<CxT1> Gradient(const Col<CxT1> &x) const;

  CxT1 Denom(const Col<CxT1> &ddir, const Col<CxT1> &x) const;

};

// Explicit Instantiation
extern template class R_lowRank<float, QuadPenalty<float>>;
extern template class R_lowRank<float, TVPenalty<float>>;
extern template class R_lowRank<double, QuadPenalty<double>>;
extern template class R_lowRank<double, TVPenalty<double>>;

#endif //PowerGrid_R_lowRank_h
