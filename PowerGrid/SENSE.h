/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [SENSE.h]

    Synopsis    [Object implementing sensitivity encoding reconstructions. ]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef PowerGrid_SENSE_hpp
#define PowerGrid_SENSE_hpp
using namespace arma;

// We are using two template types at the moment. One for the type of data to be
// processed (ie Col<cx_double>) and one for the type of G object (ie
// Gfft<Col<cx_double>>
template <typename T1, typename Tobj> class SENSE {
  typedef complex<T1> CxT1;

public:
  SENSE();

  // Class variables go here
  uword n1 = 0; // Data size
  uword n2 = 0; // Image size
  uword nc = 0; // number of coils
  Tobj *G_obj;
  Mat<CxT1> SMap; // dimensions Image size b (n1 by number of coils (nc)

  // Class constructor
  SENSE(Tobj &G, Col<CxT1> SENSEmap, uword a, uword b, uword c);

  // Overloaded operators go here

  // Forward transformation is *
  // d is the vector of data of type T1, note it is const, so we don't modify it
  // directly rather return another vector of type T1
  Col<CxT1> operator*(const Col<CxT1> &d) const;

  // For the adjoint operation, we have to weight the adjoint transform of the
  // coil data by the SENSE map.
  Col<CxT1> operator/(const Col<CxT1> &d) const;
};

template class SENSE<float, Gnufft<float>>;
template class SENSE<double, Gnufft<double>>;
template class SENSE<float, Gdft<float>>;
template class SENSE<double, Gdft<double>>;

#endif
