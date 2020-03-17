/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [pcSENSE.h]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Joseph Holtrop, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

#ifndef PowerGrid_pcSENSE_hpp
#define PowerGrid_pcSENSE_hpp

#include "PGIncludes.h"
#include "Gdft.h"
#include "Gnufft.h"
#include "TimeSegmentation.h"

using namespace arma;
//using namespace PowerGrid;

template <typename T1> class pcSENSE {
typedef std::complex<T1> CxT1;

public:
pcSENSE();

~pcSENSE();

// Class variables go here
uword Nd = 0;   // Data size  (the size of one gdft or Gnufft object, length of
                // a single shot)
uword Ni = 0;   // Image size
uword Nc = 0;   // number of coils
uword Ns = 0;   // number of shots
Mat<CxT1> SMap;   // coil sensitivity, dimensions Image size(n1) by number of
                  // coils (nc)
Mat<T1> PMap;   // shot phase, dimensions Image size(n1) by number of shots. in
                // radians.
Col<T1> FMap;   // Fieldmap
Mat<T1> Kx;     // kspace coordinates in x direction
Mat<T1> Ky;     // kspace coordinates in y direction
Mat<T1> Kz;     // kspace coordinates in z direction
Mat<T1> Tvec;   // timing vector for a single shot (all shots assumed to have
                // same timing vector)
uword Nx;
uword Ny;
uword Nz;
Col<T1> Ix;
Col<T1> Iy;
Col<T1> Iz;
CxT1 i = CxT1(0., 1.);
uword type = 1;   // 2 for min max time seg and 1 for Hanning
uword L = 20;
Gdft<T1> **AObj = NULL;
	//TimeSegmentation <T1, Gnufft<T1>> **AObj = NULL;

// Class constructor
pcSENSE(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx, uword ny, uword nz,
        uword nc, Col<T1> t, Col<CxT1> SENSEmap, Col<T1> FieldMap,
        Col<T1> ShotPhaseMap);

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
Col<CxT1> operator*(const Col<CxT1> &d) const;
// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
Col<CxT1> operator/(const Col<CxT1> &d) const;
};

// Explicit Instantiation
extern template class pcSENSE<float>;
extern template class pcSENSE<double>;
#endif
