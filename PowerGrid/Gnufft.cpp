/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gnufft.cpp]

    Synopsis    [Object that represents a non-uniform discrete Fourier
                    transform implemented via a a CPU and GPU accelerated
                    gridding implementation of the non-uniform Fast Fourier
                    Tranform (NUFFT).]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#include "Gnufft.h"

template <typename T1>
Gnufft<T1>::Gnufft(
    uword dataLength, T1 gridos, uword nx, uword ny, uword nz,
    const Col<T1> &k1, const Col<T1> &k2, const Col<T1> &k3, const Col<T1> &i1,
    const Col<T1> &i2,
    const Col<T1> &i3) // Change these arguments as you need to setup the object
{
  // cout << "Entering the class constructor for Ggrid" << endl;
  n1 = nx * ny * nz;
  n2 = dataLength;
  Nx = nx;
  Ny = ny;
  Nz = nz;
  ix = i1;
  iy = i2;
  iz = i3;
  kx = k1;
  ky = k2;
  kz = k3;
  gridOS = gridos;
  // Set Beta
  kernelWidth = 4.0;
  beta = MRI_PI *
         std::sqrt((gridOS - 0.5) * (gridOS - 0.5) *
                       (kernelWidth * kernelWidth * 4.0) / (gridOS * gridOS) -
                   0.8);

  // Deal with the LUT
  // Generating Look-Up Table
  // cout << "Calculating look up table" << endl;
  calculateLUT(beta, kernelWidth, LUT, sizeLUT);
#pragma acc enter data copyin(LUT[0 : sizeLUT])
}

// Class destructor to free LUT
template <typename T1> Gnufft<T1>::~Gnufft() {
  if (LUT) {
#pragma acc exit data delete (LUT)
#pragma acc exit data delete (LUT)
    free(LUT);
  }
}

// Overloaded methods for forward and adjoint transform
// Forward transform operation using gridding
template <typename T1>
Col<complex<T1>> Gnufft<T1>::
operator*(const Col<complex<T1>> &d) const // Don't change these arguments
{
  // cout << "Entering forward operator overload in Ggrid." << endl;
  // This is just specifying size assuming things are the same size, change as
  // necessary
  // uword dataLength = d.n_rows;

  // cout << "Seperating real and imaginary data." << endl;

  Col<T1> realData = real(d);
  Col<T1> imagData = imag(d);
  // Now we grab the data out of armadillo with the memptr() function
  // This returns a pointer of the type of the elements of the
  // array/vector/matrix/cube (3d matrix)
  // Armadillo uses column major like MATLAB and Fortran, but different from
  // 2D C++ arrays which are row major.
  // cout << "Grabbing pointers." << endl;

  T1 *realDataPtr = realData.memptr();
  T1 *imagDataPtr = imagData.memptr();

  // cout << "Allocating memory for transformed data." << endl;

  Col<T1> realXformedData(this->n2);
  Col<T1> imagXformedData(this->n2);
  // realXformedData.zeros();
  // imagXformedData.zeros();
  // cout << "Grabbing pointers for new memory." << endl;

  T1 *realXformedDataPtr = realXformedData.memptr();
  T1 *imagXformedDataPtr = imagXformedData.memptr();

  // Process data here, like calling a brute force transform, dft...
  // I assume you create the pointers to the arrays where the transformed data
  // will be stored
  // realXformedDataPtr and imagXformedDataPtr and they are of type float*
  /*
  ftCpu<T2>(realXformedDataPtr,imagXformedDataPtr,
            realDataPtr, imagDataPtr, kx.memptr(),
            ky.memptr(), kz.memptr(),
            ix.memptr(), iy.memptr(), iz.memptr(),
            FM.memptr(), t.memptr(),
            this->n2, this->n1
  );
  */
  // T2 gridOS = 2.0;
  // cout << "About to call the forward gridding routine." << endl;
  computeFd_CPU_Grid<T1>(n2, kx.memptr(), ky.memptr(), kz.memptr(), realDataPtr,
                         imagDataPtr, Nx, Ny, Nz, gridOS, realXformedDataPtr,
                         imagXformedDataPtr, kernelWidth, beta, LUT, sizeLUT);

  // To return data, we need to put our data back into Armadillo objects
  // We are telling the object how long it is because it will copy the data
  // back into managed memory
  // realXformedData(realXformedDataPtr, dataLength);
  // imagXformedData(imagXformedDataPtr, dataLength);

  // We can free the realDataXformPtr and imagDataXformPtr at this point and
  // Armadillo will manage armadillo object memory as things change size or go
  // out of scope and need to be destroyed

  Col<complex<T1>> XformedData(this->n2);
  XformedData.set_real(realXformedData);
  XformedData.set_imag(imagXformedData);

  return conv_to<Col<complex<T1>>>::from(
      XformedData); // Return a vector of type T1
}

// Adjoint transform operation
template <typename T1>
Col<complex<T1>> Gnufft<T1>::operator/(const Col<complex<T1>> &d) const {

  // uword dataLength = n2;
  // Let's trim the operations to avoid data overhead and transfers
  // Basically if we know that the data points are zero, they have no impact
  // on the transform

  uword dataLength = this->n2;

  Col<T1> realData = real(d);
  Col<T1> imagData = imag(d);

  T1 *realDataPtr = realData.memptr();
  T1 *imagDataPtr = imagData.memptr();

  Col<T1> realXformedData(n1);
  Col<T1> imagXformedData(n1);

  // realXformedData.zeros();
  // imagXformedData.zeros();

  T1 *realXformedDataPtr = realXformedData.memptr();
  T1 *imagXformedDataPtr = imagXformedData.memptr();
  // Process data here, like calling a brute force transform, dft...
  // I assume you create the pointers to the arrays where the transformed data
  // will be stored
  // realXformedDataPtr and imagXformedDataPtr and they are of type float*

  // T2 gridOS = 2.0;

  computeFH_CPU_Grid<T1>(dataLength, kx.memptr(), ky.memptr(), kz.memptr(),
                         realDataPtr, imagDataPtr, Nx, Ny, Nz, gridOS,
                         realXformedDataPtr, imagXformedDataPtr, kernelWidth,
                         beta, LUT, sizeLUT);
  /*
  iftCpu<T2>(realXformedDataPtr,imagXformedDataPtr,
             realDataPtr, imagDataPtr, kx.memptr(),
             ky.memptr(), kz.memptr(),
             ix.memptr(), iy.memptr(), iz.memptr(),
             FM.memptr(), t.memptr(),
             this->n2, this->n1
             );
  */
  // realXformedData(realXformedDataPtr, dataLength);
  // imagXformedData(imagXformedDataPtr, dataLength);

  // We can free the realDataXformPtr and imagDataXformPtr at this point and
  // Armadillo will manage armadillo object memory as things change size or go
  // out of scope and need to be destroyed

  Col<complex<T1>> XformedData(n1);
  XformedData.set_real(realXformedData);
  XformedData.set_imag(imagXformedData);
  // XformedData.elem(dataMaskTrimmed) = XformedDataTrimmed;
  // savemat("/shared/mrfil-data/data/PowerGridTest/64_64_16_4coils/ggrid.mat","img",XformedData);

  return conv_to<Col<complex<T1>>>::from(
      XformedData); // Return a vector of type T1
}

template <typename T1>
Col<complex<T1>> Gnufft<T1>::trimmedForwardOp(
    const Col<complex<T1>> &d,
    const Col<complex<T1>> &tempInterp) const // Don't change these arguments
{
  // Let's trim the operations to avoid data overhead and transfers
  // Basically if we know that the data points are zero, they have no impact
  // on the transform
  uword dataLength = this->n2;
  uvec dataMaskTrimmed = find(abs(tempInterp) > 0);

  Col<T1> kxTrimmed = kx.elem(dataMaskTrimmed);
  Col<T1> kyTrimmed = ky.elem(dataMaskTrimmed);

  Col<T1> kzTrimmed = kz.elem(dataMaskTrimmed);
  uword dataLengthTrimmed = kxTrimmed.n_rows;
  // std::cout << "Length of DataMaskTrimmed = " << dataLengthTrimmed <<
  // std::endl;
  // cout << "Separating real and imaginary data." << endl;

  Col<T1> realData = real(d);
  Col<T1> imagData = imag(d);

  // Now we grab the data out of armadillo with the memptr() function
  // This returns a pointer of the type of the elements of the
  // array/vector/matrix/cube (3d matrix)
  // Armadillo uses column major like MATLAB and Fortran, but different from
  // 2D C++ arrays which are row major.

  T1 *realDataPtr = realData.memptr();
  T1 *imagDataPtr = imagData.memptr();

  // cout << "Allocating memory for transformed data." << endl;
  Col<T1> realXformedDataTrimmed(dataLengthTrimmed);
  Col<T1> imagXformedDataTrimmed(dataLengthTrimmed);

  // realXformedData.zeros();
  // imagXformedData.zeros();
  // cout << "Grabbing pointers for new memory." << endl;

  T1 *realXformedDataPtr = realXformedDataTrimmed.memptr();
  T1 *imagXformedDataPtr = imagXformedDataTrimmed.memptr();

  // Process data here, like calling a brute force transform, dft...
  // I assume you create the pointers to the arrays where the transformed data
  // will be stored
  // realXformedDataPtr and imagXformedDataPtr and they are of type float*

  // cout << "About to call the forward gridding routine." << endl;
  computeFd_CPU_Grid<T1>(dataLengthTrimmed, kxTrimmed.memptr(),
                         kyTrimmed.memptr(), kzTrimmed.memptr(), realDataPtr,
                         imagDataPtr, Nx, Ny, Nz, gridOS, realXformedDataPtr,
                         imagXformedDataPtr, kernelWidth, beta, LUT, sizeLUT);

  // To return data, we need to put our data back into Armadillo objects
  // We are telling the object how long it is because it will copy the data
  // back into managed memory
  // realXformedData(realXformedDataPtr, dataLength);
  // imagXformedData(imagXformedDataPtr, dataLength);

  // We can free the realDataXformPtr and imagDataXformPtr at this point and
  // Armadillo will manage armadillo object memory as things change size or go
  // out of scope and need to be destroyed
  Col<complex<T1>> XformedDataTrimmed(dataLengthTrimmed);
  Col<complex<T1>> XformedData(dataLength);
  XformedDataTrimmed.set_real(realXformedDataTrimmed);
  XformedDataTrimmed.set_imag(imagXformedDataTrimmed);
  XformedData.elem(dataMaskTrimmed) = XformedDataTrimmed;

  return conv_to<Col<complex<T1>>>::from(
      XformedData); // Return a vector of type T1
}

// Adjoint transform operation
template <typename T1>
Col<complex<T1>>
Gnufft<T1>::trimmedAdjointOp(const Col<complex<T1>> &d,
                             const Col<complex<T1>> &tempInterp) const {

  // uword dataLength = n2;
  // Let's trim the operations to avoid data overhead and transfers
  // Basically if we know that the data points are zero, they have no impact
  // on the transform

  uword dataLength = this->n2;
  uvec dataMaskTrimmed = find(abs(tempInterp) > 0);
  uword dataLengthTrimmed = dataMaskTrimmed.n_rows;
  // std::cout << "Length of DataMaskTrimmed = " << dataLengthTrimmed <<
  // std::endl;

  Col<complex<T1>> dTrimmed = d.elem(dataMaskTrimmed);
  Col<T1> kxTrimmed = kx.elem(dataMaskTrimmed);
  Col<T1> kyTrimmed = ky.elem(dataMaskTrimmed);
  Col<T1> kzTrimmed = kz.elem(dataMaskTrimmed);

  Col<T1> realData = real(dTrimmed);
  Col<T1> imagData = imag(dTrimmed);

  T1 *realDataPtr = realData.memptr();
  T1 *imagDataPtr = imagData.memptr();

  Col<T1> realXformedData(n1);
  Col<T1> imagXformedData(n1);

  // realXformedData.zeros();
  // imagXformedData.zeros();

  T1 *realXformedDataPtr = realXformedData.memptr();
  T1 *imagXformedDataPtr = imagXformedData.memptr();
  // Process data here, like calling a brute force transform, dft...
  // I assume you create the pointers to the arrays where the transformed data
  // will be stored
  // realXformedDataPtr and imagXformedDataPtr and they are of type float*

  // T2 gridOS = 2.0;

  computeFH_CPU_Grid<T1>(dataLengthTrimmed, kxTrimmed.memptr(),
                         kyTrimmed.memptr(), kzTrimmed.memptr(), realDataPtr,
                         imagDataPtr, Nx, Ny, Nz, gridOS, realXformedDataPtr,
                         imagXformedDataPtr, kernelWidth, beta, LUT, sizeLUT);
  /*
  iftCpu<T2>(realXformedDataPtr,imagXformedDataPtr,
             realDataPtr, imagDataPtr, kx.memptr(),
             ky.memptr(), kz.memptr(),
             ix.memptr(), iy.memptr(), iz.memptr(),
             FM.memptr(), t.memptr(),
             this->n2, this->n1
             );
  */
  // realXformedData(realXformedDataPtr, dataLength);
  // imagXformedData(imagXformedDataPtr, dataLength);

  // We can free the realDataXformPtr and imagDataXformPtr at this point and
  // Armadillo will manage armadillo object memory as things change size or go
  // out of scope and need to be destroyed

  Col<complex<T1>> XformedData(n1);
  XformedData.set_real(realXformedData);
  XformedData.set_imag(imagXformedData);
  // XformedData.elem(dataMaskTrimmed) = XformedDataTrimmed;
  // savemat("/shared/mrfil-data/data/PowerGridTest/64_64_16_4coils/ggrid.mat","img",XformedData);

  return conv_to<Col<complex<T1>>>::from(
      XformedData); // Return a vector of type T1
}

// Explicit Instantiation
template class Gnufft<float>;
template class Gnufft<double>;
