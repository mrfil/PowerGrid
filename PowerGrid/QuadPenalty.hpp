/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [QuadPenalty.hpp]

    Synopsis    [Quadratic Penalty using Robject interface to maintain
                    commonality with the Image Reconstruction Toolkit (IRT).]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef PowerGrid_QuadPenalty_hpp
#define PowerGrid_QuadPenalty_hpp

//#include <armadillo>

using namespace arma;

template <typename T1>
class QuadPenalty: public Robject<T1>
{
typedef complex<T1> CxT1;
public:
    QuadPenalty();

    // It was declared as type Mat<uword> and the 3D type was a cube. We need to vectorize it before it is passed to QuadPenalty.
    //Custom Class Constructor
    QuadPenalty(uword nx,uword ny,uword nz, double beta)
    {
        //Set Class Memebers
        this->Nx = nx;
        this->Ny = ny;
        this->Nz = nz;
        this->Beta = beta;

    }

    //Class Methods

    Col<CxT1> wpot(const Col<CxT1>& d) const
    {
        return ones < Col < CxT1 >> (d.n_rows);
    }

    Col <CxT1> dpot(const Col <CxT1>& d) const
    {
        return d;
    }
    Col<CxT1> pot(const Col< CxT1 >& d) const
    {
        Col <T1> temp = abs(d)%abs(d)/2.0;
        return conv_to<Col<CxT1 >> ::from(temp);
    }


};

#endif
