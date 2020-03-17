/*
   (C) Copyright 2015-2020 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [LRobj.cpp]

    Synopsis    [Class to implement partial separability model via low rank
                 representations. ]

    Description []

    Revision    [0.1.0; Joseph Holtrop, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

#include "LRobj.h"

template <typename T1, typename Gobj>
LRobj<T1,Gobj>::LRobj(uword NimgSize, uword N_data, uword n_rank, uword n_img, Mat<CxT1> vBasis, Gobj** All) {
    this->Nimg = n_img;
    this->N = NimgSize;
    this->Nrank = n_rank;
    this->Ndata = N_data;
    this->v = vBasis;
    this->A_all = All;
}

template <typename T1, typename Gobj>
Col<complex<T1> > LRobj<T1,Gobj>::operator*(const Col<complex<T1> > &d) const {
    RANGE("pcSENSE::operator*")
        
    Mat<complex<T1>> tempData = reshape(d,this->N, this->Nrank);
    Mat<complex<T1>> tempOut(this->Ndata, this->Nimg);
    Mat<complex<T1>> tempFullRank = tempData * this->v;

    for(uword ii = 0; ii < this->Nimg; ii++) {
        tempOut.col(ii) = (*A_all[ii]) * tempFullRank.col(ii);
    }

    return vectorise(tempOut);
}

template <typename T1, typename Gobj>
Col<complex<T1> > LRobj<T1,Gobj>::operator/(const Col<complex<T1> > &d) const {

    Mat<complex<T1>> tempData = reshape(d, this->Ndata, this->Nimg);
    Mat<complex<T1>> tempOut = zeros<Mat<complex<T1>>>(this->N, this->Nimg);

    for(uword ii = 0; ii < this->Nimg; ii++) {
        tempOut.col(ii) = (*A_all[ii]) / tempData.col(ii);
    }
    Mat<complex<T1>> compressedOut = tempOut * trans(this->v);

    return vectorise(compressedOut);
}

template class LRobj<float, pcSENSE<float>>;
template class LRobj<float, pcSenseTimeSeg<float>>;
template class LRobj<double, pcSENSE<double>>;
template class LRobj<double, pcSenseTimeSeg<double>>;
