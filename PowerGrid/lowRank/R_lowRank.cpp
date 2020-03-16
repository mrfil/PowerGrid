/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [R_lowRank.cpp]

    Synopsis    [Wrapper class to extend penalization in the low rank recon case.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [06/16/2020]

*****************************************************************************/

#include "R_lowRank.h"
#include <chrono>  // for high_resolution_clock


template <typename T1, typename Robj>
R_lowRank<T1, Robj>::R_lowRank(Robj const &R1, uword rank, uword nImages, Mat<CxT1> vBasis) {
        // Set Class Memebers
        this->R = &R1;
        this->Nrank = rank;
        this->Nimg = nImages;
        this->v = vBasis;
}

template <typename T1, typename Robj>
T1 R_lowRank<T1, Robj>::Penalty(const Col<complex<T1> > &x) const {
        RANGE()
        
        Col<CxT1> imgOut(Nimg);

        uword N = x.n_rows/this->Nrank;
        Mat<CxT1> spatialBasis = reshape(x, N, this->Nrank);
        Mat<CxT1> imgStack = reshape(spatialBasis * trans(this->v.cols(0,Nrank-1)), N, Nimg);

        for(uword ii = 0; ii < Nimg; ii++) {
            imgOut(ii) = this->R->Penalty(imgStack.col(ii));
        }

        return abs(sum(imgOut));
}
template <typename T1, typename Robj>
Col<complex<T1>> R_lowRank<T1, Robj>::Gradient(const Col<complex<T1> > &x) const {
        RANGE()
        
        uword N = x.n_rows/this->Nrank;
        Mat<CxT1> spatialBasis = reshape(x, N, this->Nrank);
        Mat<CxT1> imgStack = reshape(spatialBasis * trans(this->v.cols(0,Nrank-1)), N, Nimg);
        Mat<CxT1> yOut(N,this->Nimg);
        
        for(uword ii = 0; ii < Nimg; ii++) {
            yOut.col(ii) = this->R->Gradient(imgStack.col(ii));
        }

        yOut = strans(this->v.cols(0,(this->Nrank-1))) * strans(yOut);
        
        return vectorise(yOut);
}

template <typename T1, typename Robj>
complex<T1> R_lowRank<T1, Robj>::Denom(const Col<complex<T1> > &ddir,
                               const Col<complex<T1> > &x) const {
        RANGE()
        

        uword N = x.n_rows/this->Nrank;
        Mat<CxT1> spatialBasis = reshape(x, N, this->Nrank);
        Mat<CxT1> imgStack = reshape(spatialBasis * trans(this->v.cols(0,Nrank-1)), N, Nimg);
        
        Mat<CxT1> ddirReshape = reshape(ddir,N,this->Nrank);
        Mat<CxT1> spatialDdir = reshape(ddirReshape * trans(this->v.cols(0,this->Nrank-1)), N, Nimg);

        Col<CxT1> temp = zeros<Col<CxT1>>(this->Nimg);

        for(uword ii = 0; ii < this->Nimg; ii++) {
            temp(ii) = this->R->Denom(spatialDdir.col(ii), imgStack.col(ii));
        }

        return this->R->Beta * sum(temp);
}

// Explicit Instantiation
template class R_lowRank<float, QuadPenalty<float>>;
template class R_lowRank<float, TVPenalty<float>>;
template class R_lowRank<double, QuadPenalty<double>>;
template class R_lowRank<double, TVPenalty<double>>;



