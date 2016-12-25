/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gdft.cpp]

    Synopsis    [Object that represents a non-uniform field corrected discrete
                    Fourier tranform.]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/2/2016]

 *****************************************************************************/

using namespace arma;

template <typename T1>
Gdft::Gdft(
    uword a, uword b, const Col<T1> &k1, const Col<T1> &k2, const Col<T1> &k3,
    const Col<T1> &i1, const Col<T1> &i2, const Col<T1> &i3, const Col<T1> &f1,
    const Col<T1> &t1) // Change these arguments as you need to setup the object
{
  n1 = a;
  n2 = b;
  kx = k1;
  ky = k2;
  kz = k3;
  ix = i1;
  iy = i2;
  iz = i3;
  FM = f1;
  t = t1;
}

// Overloaded methods for forward and adjoint transform
// Forward transform operation
template <typename CxT1> Col<CxT1> Gdft::operator*(const Col<CxT1> &d) const {
  // This is just specifying size assuming things are the same size, change as
  // necessary
  Col<T1> realData = real(d);
  Col<T1> imagData = imag(d);
  // Now we grab the data out of armadillo with the memptr() function
  // This returns a pointer of the type of the elements of the
  // array/vector/matrix/cube (3d matrix)
  // Armadillo uses column major like MATLAB and Fortran, but different from
  // 2D C++ arrays which are row major.
  T1 *realDataPtr = realData.memptr();
  T1 *imagDataPtr = imagData.memptr();

  Col<T1> realXformedData;
  Col<T1> imagXformedData;
  realXformedData.zeros(this->n1);
  imagXformedData.zeros(this->n1);

  T1 *realXformedDataPtr = realXformedData.memptr();
  T1 *imagXformedDataPtr = imagXformedData.memptr();
  // Process data here, like calling a brute force transform, dft...
  // I assume you create the pointers to the arrays where the transformed data
  // will be stored
  // realXformedDataPtr and imagXformedDataPtr and they are of type float*
  ftCpu<T1>(realXformedDataPtr, imagXformedDataPtr, realDataPtr, imagDataPtr,
            kx.memptr(), ky.memptr(), kz.memptr(), ix.memptr(), iy.memptr(),
            iz.memptr(), FM.memptr(), t.memptr(), this->n1, this->n2);

  // To return data, we need to put our data back into Armadillo objects
  // We are telling the object how long it is because it will copy the data
  // back into managed memory
  // realXformedData(realXformedDataPtr, dataLength);
  // imagXformedData(imagXformedDataPtr, dataLength);

  // We can free the realDataXformPtr and imagDataXformPtr at this point and
  // Armadillo will manage armadillo object memory as things change size or go
  // out of scope and need to be destroyed

  Col<CxT1> XformedData(this->n1);
  XformedData.set_real(realXformedData);
  XformedData.set_imag(imagXformedData);

  return conv_to<Col<CxT1>>::from(XformedData); // Return a vector of type T1
}

// Adjoint transform operation
template <typename CxT1> Col<CxT1> Gdft::operator/(const Col<CxT1> &d) const {

  Col<T1> realData = real(d);
  Col<T1> imagData = imag(d);

  T1 *realDataPtr = realData.memptr();
  T1 *imagDataPtr = imagData.memptr();

  Col<T1> realXformedData;
  Col<T1> imagXformedData;
  realXformedData.zeros(this->n2);
  imagXformedData.zeros(this->n2);

  T1 *realXformedDataPtr = realXformedData.memptr();
  T1 *imagXformedDataPtr = imagXformedData.memptr();
  // Process data here, like calling a brute force transform, dft...
  // I assume you create the pointers to the arrays where the transformed data
  // will be stored
  // realXformedDataPtr and imagXformedDataPtr and they are of type float*
  iftCpu<T1>(realXformedDataPtr, imagXformedDataPtr, realDataPtr, imagDataPtr,
             kx.memptr(), ky.memptr(), kz.memptr(), ix.memptr(), iy.memptr(),
             iz.memptr(), FM.memptr(), t.memptr(), this->n1, this->n2);

  // realXformedData(realXformedDataPtr, dataLength);
  // imagXformedData(imagXformedDataPtr, dataLength);

  // We can free the realDataXformPtr and imagDataXformPtr at this point and
  // Armadillo will manage armadillo object memory as things change size or go
  // out of scope and need to be destroyed

  Col<CxT1> XformedData(this->n2);
  XformedData.set_real(realXformedData);
  XformedData.set_imag(imagXformedData);

  return conv_to<Col<CxT1>>::from(XformedData); // Return a vector of type T1
}

template class Gdft<std::complex<float>>;
template class Gdft<std::complex<double>>;
