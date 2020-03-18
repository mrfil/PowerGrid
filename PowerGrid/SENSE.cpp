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
#include "SENSE.h"

using namespace arma;

// We are using two template types at the moment. One for the type of data to be
// processed (ie Col<cx_double>) and one for the type of G object (ie
// Gfft<Col<cx_double>>
template <typename T1, typename Tobj>
SENSE<T1, Tobj>::SENSE(Tobj& G, Col<complex<T1>> SENSEmap, uword a, uword b,
    uword c)
{
    n1 = a;
    n2 = b;
    nc = c;
    G_obj = &G;
    SMap = reshape(SENSEmap, n2, nc);
    conjSMap = conj(SMap);
    outData.set_size(n1, nc);
    outImg.set_size(n2, nc);
}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1, typename Tobj>
Col<complex<T1>> SENSE<T1, Tobj>::operator*(const Col<complex<T1>>& d) const
{
    RANGE("SENSE::operator*")
    //Mat<complex<T1>> outData = zeros<Mat<complex<T1>>>(this->n1, this->nc);
    // Col<complex<T1>> temp;
    // In SENSE we store coil data using the columns of the data matrix, and we
    // weight the data by the coil sensitivies from the SENSE map
    outImg = this->SMap;
#pragma omp parallel for schedule(dynamic) shared(outData, d, SMap)
    for (unsigned int ii = 0; ii < this->nc; ii++) {
        outImg.unsafe_col(ii) %= d;
    }

    for (unsigned int ii = 0; ii < this->nc; ii++) {
        outData.unsafe_col(ii) = (*this->G_obj) * outImg.unsafe_col(ii);
    }

    // equivalent to returning col(output) in MATLAB with IRT
    return vectorise(outData);
}

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
template <typename T1, typename Tobj>
Col<complex<T1>> SENSE<T1, Tobj>::operator/(const Col<complex<T1>>& d) const
{
    RANGE("SENSE::operator/")
    //Mat<complex<T1>> inData = reshape(d, this->n1, this->nc);

    //Col<complex<T1>> outData = zeros<Col<complex<T1>>>(this->n2);
    // Mat <complex<T1>> coilImages(n2,nc);

    for (unsigned int ii = 0; ii < this->nc; ii++) {
        // coilImages.col(ii) = (*this->G_obj)/inData.col(ii);
        outImg.unsafe_col(ii) = (*this->G_obj) / d.subvec((ii)*n1, ((ii + 1) * n1) - 1);
    }

#pragma omp parallel for schedule(dynamic) shared(outImg, conjSMap)
    for (unsigned int ii = 0; ii < this->nc; ii++) {
        outImg.unsafe_col(ii) %= this->conjSMap.unsafe_col(ii);
    }
    // outData = sum(conj(SMap)%coilImages,2);

    // equivalent to returning col(output) in MATLAB with IRT
    return sum(outImg, 1);
}

// Explicit Instantiations
template class SENSE<float, Gnufft<float>>;
template class SENSE<float, TimeSegmentation<float, Gnufft<float>>>;
template class SENSE<double, Gnufft<double>>;
template class SENSE<double, TimeSegmentation<double, Gnufft<double>>>;
template class SENSE<float, Gdft<float>>;
template class SENSE<double, Gdft<double>>;
template class SENSE<float, GdftR2<float>>;
template class SENSE<double, GdftR2<double>>;
