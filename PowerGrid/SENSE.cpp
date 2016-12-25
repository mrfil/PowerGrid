/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [SENSE.cpp]

    Synopsis    [Object implementing sensitivity encoding reconstructions. ]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

using namespace arma;

// We are using two template types at the moment. One for the type of data to be
// processed (ie Col<cx_double>) and one for the type of G object (ie
// Gfft<Col<cx_double>>
template <typename T1, typename Tobj>
SENSE<T1, Tobj>::SENSE(Tobj &G, Col<CxT1> SENSEmap, uword a, uword b, uword c) {
  n1 = a;
  n2 = b;
  nc = c;
  G_obj = &G;
  SMap = reshape(SENSEmap, n2, nc);
}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1, typename Tobj>
Col<CxT1> SENSE<T1, Tobj>::operator*(const Col<CxT1> &d) const {

  Mat<CxT1> outData = zeros<Mat<CxT1>>(this->n1, this->nc);
  // Col<CxT1> temp;
  // In SENSE we store coil data using the columns of the data matrix, and we
  // weight the data by the coil sensitivies from the SENSE map

  for (unsigned int ii = 0; ii < this->nc; ii++) {

    outData.col(ii) = (*this->G_obj) * (d % (this->SMap.col(ii)));
  }
  Col<CxT1> out = vectorise(outData);
  // equivalent to returning col(output) in MATLAB with IRT
  return out;
}

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
template <typename T1, typename Tobj>
Col<CxT1> SENSE<T1, Tobj>::operator/(const Col<CxT1> &d) const {

  Mat<CxT1> inData = reshape(d, this->n1, this->nc);

  Col<CxT1> outData = zeros<Col<CxT1>>(this->n2);
  // Mat <CxT1> coilImages(n2,nc);

  for (unsigned int ii = 0; ii < this->nc; ii++) {
    // coilImages.col(ii) = (*this->G_obj)/inData.col(ii);
    outData += conj(this->SMap.col(ii)) % ((*this->G_obj) / inData.col(ii));
  }
  // outData = sum(conj(SMap)%coilImages,2);
  // equivalent to returning col(output) in MATLAB with IRT
  return vectorise(outData);
}
