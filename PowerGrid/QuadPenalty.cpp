/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [QuadPenalty.cpp]

    Synopsis    [Quadratic Penalty using Robject interface to maintain
                    commonality with the Image Reconstruction Toolkit (IRT).]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

#include "QuadPenalty.h"

// It was declared as type Mat<uword> and the 3D type was a cube. We need to
// vectorize it before it is passed to QuadPenalty.
// Custom Class Constructor
template <typename T1>
QuadPenalty<T1>::QuadPenalty(uword nx, uword ny, uword nz, T1 beta, uword dims2penalize) {

        // Set Class Memebers
        this->Nx = nx;
        this->Ny = ny;
        this->Nz = nz;
        this->Beta = beta;
        this->Dims2Penalize = dims2penalize;
}

// Class Methods

template <typename T1>
inline
Col<complex<T1> > QuadPenalty<T1>::wpot(const Col<complex<T1> > &d) const {
        RANGE()
        return ones<Col<complex<T1> > >(d.n_rows);
}

template <typename T1>
inline
Col<complex<T1> > QuadPenalty<T1>::dpot(const Col<complex<T1> > &d) const {
        RANGE()
        return d;
}
template <typename T1>
inline
Col<complex<T1> > QuadPenalty<T1>::pot(const Col<complex<T1> > &d) const {
        RANGE()
        Col<T1> temp = abs(d) % abs(d) / 2.0;
        return conv_to<Col<complex<T1> > >::from(temp);
}

// Explicit Instantiation
template class QuadPenalty<float>;
template class QuadPenalty<double>;
