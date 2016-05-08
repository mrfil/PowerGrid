/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [TVPenalty.hpp]

    Synopsis    [Implementation of an approximate total variational penalty.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/


#ifndef PowerGrid_TVPenalty_hpp
#define PowerGrid_TVPenalty_hpp

//#include <armadillo>

using namespace arma;

template<typename T1>
class TVPenalty : public Robject<T1> {
    typedef complex <T1> CxT1;
public:
    TVPenalty();

    // It was declared as type Mat<uword> and the 3D type was a cube. We need to vectorize it before it is passed to QuadPenalty.
    //Custom Class Constructor
    TVPenalty(uword nx, uword ny, uword nz, T1 beta, T1 delta)
    {
        //Set Class Members
        this->Nx = nx;
        this->Ny = ny;
        this->Nz = nz;
        this->Beta = beta;
        this->Delta = delta;

    }

    //Class Methods

    Col <CxT1> wpot(const Col <CxT1>& d) const
    {
        Col <T1> temp = abs(d/this->Delta);
        Col <CxT1> out = 1.0/sqrt(1.0+conv_to<Col<CxT1 >> ::from(temp%temp));
        return out;
    }

    virtual Col <CxT1> dpot(const Col <CxT1>& d) const
    {
        Col <T1> temp = abs(d/this->Delta);
        Col <CxT1> out = d/sqrt(1.0+conv_to<Col<CxT1 >> ::from(temp%temp));
        return out;
    }

    Col <CxT1> pot(const Col <CxT1>& d) const
    {
        Col <CxT1> out = conv_to<Col
                <CxT1 >> ::from(this->Delta*this->Delta*(sqrt(1.0+(abs(d/this->Delta)%abs(d/this->Delta)))-1.0));
        return out;
    }

private:
    T1 Delta;
};

#endif //PowerGrid_TVPenalty_hpp
