/*
(C) Copyright 2015-2020 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [LRobj.h]

    Synopsis    [Class to implement partial separability model via low rank
                 representations. ]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [03/16/2020]

 *****************************************************************************/

#ifndef PowerGrid_LRobj_h
#define PowerGrid_LRobj_h

#include "PGIncludes.h"
#include "pcSENSE.h"
#include "pcSenseTimeSeg.h"

using namespace arma;
using namespace std;


template <typename T1, typename Gobj> class LRobj {
  typedef complex<T1> CxT1;

public:
    Mat<CxT1> v; // Temporal basis functions
    Gobj **A_all = NULL; // List of objects
    uword Nrank;    // Number of ranks
    uword Nimg; // Number of images
    uword Ndata; // Shot length
    uword N; // Image size

    LRobj(uword NimgSize, uword Ndata, uword n_rank, uword n_img, Mat<CxT1> vBasis, Gobj** All);

    Col<CxT1> operator*(const Col<CxT1> &d) const;

    Col<CxT1> operator/(const Col<CxT1> &d) const;

};


#endif //PowerGrid_LRobj_h

extern template class LRobj<float, pcSENSE<float>>;
extern template class LRobj<float, pcSenseTimeSeg<float>>;
extern template class LRobj<double, pcSENSE<double>>;
extern template class LRobj<double, pcSenseTimeSeg<double>>;
