/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [griddingSupport.cpp]

    Synopsis    [Support functions for the forward and adjoint non-uniform
                    Fast Fourier Transform (NUFFT) on CPU and GPU via
                    OpenACC.]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/2/2016]

*****************************************************************************/

#include "griddingSupport.h"

///*From Numerical Recipes in C, 2nd Edition
// Just a vanilla I(0,x) function approximation
#pragma acc routine seq
template <typename T1> inline T1 bessi0(T1 x) {
  T1 ax, ans;
  T1 y;

  if ((ax = std::abs(x)) < 3.75) {
    y = x / 3.75;
    y = y * y;
    ans =
        1.0 +
        y * (3.5156229 +
             y * (3.0899424 +
                  y * (1.2067492 +
                       y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
  } else {
    y = 3.75 / ax;
    ans = (std::exp(ax) / std::sqrt(ax)) *
          (0.39894228 +
           y * (0.1328592e-1 +
                y * (0.225319e-2 +
                     y * (-0.157565e-2 +
                          y * (0.916281e-2 +
                               y * (-0.2057706e-1 +
                                    y * (0.2635537e-1 +
                                         y * (-0.1647633e-1 +
                                              y * 0.392377e-2))))))));
  }
  return ans;
}

// We use this lookup table as a faster way to calculate the Kasier-Bessel
// window
template <typename T1>
void calculateLUT(T1 beta, T1 width, T1 *&LUT, uword &sizeLUT) {
  T1 v;
  // T1 _width2_4 = (width*width)/4.0;

  if (width > 0) {
    // compute size of LUT based on kernel width
    sizeLUT = (uword)(10000 * width);
    // allocate memory

    LUT = (T1 *)malloc(sizeLUT * sizeof(T1));
    //#pragma acc enter data create(LUT[0:sizeLUT])

    for (int k = 0; k < sizeLUT; ++k) {
      // compute value to evaluate kernel at
      // v in the range 0:(_width/2)^2
      v = static_cast<T1>(k) / static_cast<T1>(sizeLUT);

      // compute kernel value and store
      LUT[k] = bessi0(beta * std::sqrt(1.0 - (v))) / width;
    }
  }
}
/*
#pragma acc routine seq
template <typename T1>
T1 kernel_value_LUT(T1 dist, const T1 *LUT, uword sizeLUT,
                    T1 width) { // v is between [0,width/2.0]
  uword k0;
  // T1 v0;
  T1 _width2_4 =
      4.0 /
      (width * width); // Reciprocal of _width2_4 from calculateLUT function

  k0 = (uword)((dist * dist * _width2_4) * (T1)sizeLUT);

  // cout << "dist = " << dist << " k0 =" << k0 << endl;

  if (k0 >= sizeLUT) {
    return 0;
  }  else {
   //cout << "about to access the look up table" << endl;
  return LUT[k0];
  }
}
*/
template <typename T1>
void deinterleave_data2d( /// NAIVE
    T1 *__restrict pSrc, T1 *__restrict outR_d, T1 *__restrict outI_d,
    int imageX, int imageY) {
  int lIndex;
#pragma acc parallel loop collapse(2) independent present(                     \
    pSrc[0 : 2 * imageX *imageY])                                              \
        present(outR_d[0 : imageX *imageY], outI_d[0 : imageX *imageY])
  for (int X = 0; X < imageX; X++) {
    for (int Y = 0; Y < imageY; Y++) {
      lIndex = Y + X * imageY;
      outR_d[lIndex] = pSrc[2 * lIndex];
      outI_d[lIndex] = pSrc[2 * lIndex + 1];
    }
  }
}

template <typename T1>
void deinterleave_data3d(T1 *__restrict pSrc, T1 *__restrict outR_d,
                         T1 *__restrict outI_d, int imageX, int imageY,
                         int imageZ) {
  int lIndex, X, Y, Z;
#pragma acc parallel loop collapse(                                            \
    3) independent present(pSrc[0 : 2 * imageX *imageY *imageZ]) present(      \
    outR_d[0 : imageX *imageY *imageZ], outI_d[0 : imageX *imageY *imageZ])
  for (Z = 0; Z < imageZ; Z++) {
    for (X = 0; X < imageX; X++) {
      for (Y = 0; Y < imageY; Y++) {
        lIndex = Y + X * imageY + Z * imageY * imageX;
        outR_d[lIndex] = pSrc[2 * lIndex];
        outI_d[lIndex] = pSrc[2 * lIndex + 1];
      }
    }
  }
}

// Deapodizes 2d data by FT of the Kasier-Bessel kernel
template <typename T1>
void deapodization2d(T1 *__restrict pDst, T1 *__restrict pSrc, int imageX,
                     int imageY, T1 kernelWidth, T1 beta, T1 gridOS) {

  // int imageNumElems = imageX * imageY;

  int Y;
  int X;

  T1 gridKernelY;
  T1 gridKernelX;

  T1 common_exprX;
  T1 common_exprY;

  T1 common_exprX1;
  T1 common_exprY1;

  T1 gridKernel;
  int common_index;
  T1 gridOS2;

  int destSize = imageX * imageY;

#pragma acc kernels present(pDst[0 : 2 * imageX *imageY])
  {
#pragma acc loop
    for (int ii = 0; ii < 2 * destSize; ii++) {
      pDst[ii] = (T1)0.0;
    }
  }
#pragma acc parallel loop collapse(2) independent present(                     \
    pSrc[0 : 2 * imageX *imageY]) present(pDst[0 : 2 * imageX *imageY])
  for (X = 0; X < imageX; X++) {
    for (Y = 0; Y < imageY; Y++) {

      gridKernelY = (T1)(Y - (imageY / 2)) / (T1)imageY;
      gridKernelX = (T1)(X - (imageX / 2)) / (T1)imageX;

      common_exprX = (MRI_PI * MRI_PI * kernelWidth * kernelWidth *
                          gridKernelX * gridKernelX -
                      beta * beta);
      common_exprY = (MRI_PI * MRI_PI * kernelWidth * kernelWidth *
                          gridKernelY * gridKernelY -
                      beta * beta);

      if (common_exprX >= 0)
        common_exprX1 =
            (SIN(std::sqrt(common_exprX)) / std::sqrt(common_exprX));
      else
        common_exprX1 = (SINH(std::sqrt((T1)-1.0 * common_exprX)) /
                         std::sqrt((T1)-1.0 * common_exprX));

      if (common_exprY >= 0)
        common_exprY1 =
            (SIN(std::sqrt(common_exprY)) / std::sqrt(common_exprY));
      else
        common_exprY1 = (SINH(std::sqrt((T1)-1.0 * common_exprY)) /
                         std::sqrt((T1)-1.0 * common_exprY));

      gridKernel = common_exprX1 * common_exprY1;

      if (gridKernel == gridKernel) // Check for NaN
      {

        common_index = Y + X * imageY;
        gridOS2 = gridOS * gridOS;

        // dst[common_index] =
        // (src[common_index]) / gridKernel *(1.0 / gridOS2);
        // //???
        pDst[2 * common_index] =
            (pSrc[2 * common_index]) / gridKernel * ((T1)1.0 / gridOS2); // Real
        pDst[2 * common_index + 1] = (pSrc[2 * common_index + 1]) / gridKernel *
                                     ((T1)1.0 / gridOS2); // Imaginary
      }
    }
  }
}

// Deapodizes 3d data by FT of the Kasier-Bessel kernel
template <typename T1>
void deapodization3d(T1 *__restrict pDst, T1 *__restrict pSrc, int imageX,
                     int imageY, int imageZ, T1 kernelWidth, T1 beta,
                     T1 gridOS) {

  int Z;
  int Y;
  int X;

  // int imageNumElems = imageX * imageY * imageZ;

  T1 gridKernelZ;
  T1 gridKernelX;
  T1 gridKernelY;

  T1 common_exprZ;
  T1 common_exprY;
  T1 common_exprX;

  T1 common_exprX1;
  T1 common_exprY1;
  T1 common_exprZ1;

  T1 gridOS3;
  // int common_index;
  int destSize = imageX * imageY * imageZ;
#pragma acc kernels present(pDst[0 : 2 * imageX *imageY *imageZ])
  {
#pragma acc loop
    for (int ii = 0; ii < 2 * destSize; ii++) {
      pDst[ii] = (T1)0.0;
    }
  }
#pragma acc parallel loop collapse(3)                                          \
    independent present(pSrc[0 : 2 * imageX *imageY *imageZ])                  \
                            present(pDst[0 : 2 * imageX *imageY *imageZ])
  for (Z = 0; Z < imageZ; Z++) {
    for (X = 0; X < imageX; X++) {
      for (Y = 0; Y < imageY; Y++) {

        gridKernelZ = (T1)((Z) - ((T1)imageZ / 2.0)) / (T1)imageZ;
        gridKernelY = (T1)((Y) - ((T1)imageY / 2.0)) / (T1)imageY;
        gridKernelX = (T1)((X) - ((T1)imageX / 2.0)) / (T1)imageX;

        common_exprX = (MRI_PI * MRI_PI * kernelWidth * kernelWidth *
                            gridKernelX * gridKernelX -
                        beta * beta);
        common_exprY = (MRI_PI * MRI_PI * kernelWidth * kernelWidth *
                            gridKernelY * gridKernelY -
                        beta * beta);
        common_exprZ = (MRI_PI * MRI_PI * kernelWidth * kernelWidth *
                            gridKernelZ * gridKernelZ -
                        beta * beta);

        if (common_exprX >= 0)
          common_exprX1 =
              (SIN(std::sqrt(common_exprX)) / std::sqrt(common_exprX));
        else
          common_exprX1 = (SINH(std::sqrt(-1.0 * common_exprX)) /
                           std::sqrt(-1.0 * common_exprX));

        if (common_exprY >= 0)
          common_exprY1 =
              (SIN(std::sqrt(common_exprY)) / std::sqrt(common_exprY));
        else
          common_exprY1 = (SINH(std::sqrt(-1.0 * common_exprY)) /
                           std::sqrt(-1.0 * common_exprY));

        if (common_exprZ >= 0)
          common_exprZ1 =
              (SIN(std::sqrt(common_exprZ)) / std::sqrt(common_exprZ));
        else
          common_exprZ1 = (SINH(std::sqrt(-1.0 * common_exprZ)) /
                           std::sqrt(-1.0 * common_exprZ));

        T1 gridKernel = common_exprX1 * common_exprY1 * common_exprZ1;

        if (gridKernel == gridKernel) // Checking for NaN
        {
          int common_index = Z * imageY * imageX + X * imageY + Y;
          gridOS3 = gridOS * gridOS * gridOS;
          // dst[common_index] =
          //(src[common_index]) / gridKernel *(1.0 / gridOS3);
          pDst[2 * common_index] =
              (pSrc[2 * common_index]) / gridKernel * (1.0 / gridOS3); // Real
          pDst[2 * common_index + 1] = (pSrc[2 * common_index + 1]) /
                                       gridKernel *
                                       (1.0 / gridOS3); // Imaginary
        }
      }
    }
  }
}

// We oversample the FFT by zeropadding so now we need to crop
template <typename T1>
void crop_center_region2d(T1 *__restrict pDst, T1 *__restrict pSrc,
                          int imageSizeX, int imageSizeY, int gridSizeX,
                          int gridSizeY) {
  // (gridSizeX,gridSizeY) is the size of 'src'

  int offsetY;
  int offsetX;
  int dY_src;
  int dX_src;
  int common_index_dst;
  int common_index_src;

#pragma acc parallel loop collapse(2)                                          \
    independent present(pSrc[0 : 2 * gridSizeX *gridSizeY])                    \
                            present(pDst[0 : 2 * imageSizeX *imageSizeY])
  for (int dX_dst = 0; dX_dst < imageSizeX; dX_dst++) {
    for (int dY_dst = 0; dY_dst < imageSizeY; dY_dst++) {

      offsetY = (int)(((T1)gridSizeY / 2.0f) - ((T1)imageSizeY / 2.0f));
      offsetX = (int)(((T1)gridSizeX / 2.0f) - ((T1)imageSizeX / 2.0f));

      dY_src = dY_dst + offsetY;
      dX_src = dX_dst + offsetX;

      common_index_dst = dX_dst * imageSizeY + dY_dst;
      common_index_src = dX_src * gridSizeY + dY_src;

      pDst[2 * common_index_dst] = pSrc[2 * common_index_src];
      pDst[2 * common_index_dst + 1] = pSrc[2 * common_index_src + 1];
    }
  }
}

// We oversample the FFT so now we need to crop
template <typename T1>
void crop_center_region3d(T1 *__restrict pDst, T1 *__restrict pSrc,
                          int imageSizeX, int imageSizeY, int imageSizeZ,
                          int gridSizeX, int gridSizeY, int gridSizeZ) {
  // (gridSizeX,gridSizeY,gridSizeZ) is the size of 'src'
  //    assert( (!(gridSizeX%2) && !(gridSizeY%2) && !(gridSizeZ%2)) );
  //    assert( (!(imageSizeX%2) && !(imageSizeY%2) && !(imageSizeZ%2)) );

  int offsetY;
  int offsetX;
  int offsetZ;

  int dY_src;
  int dX_src;
  int dZ_src;
  int common_index_dst;
  int common_index_src;

#pragma acc parallel loop collapse(3) independent present(                     \
    pSrc[0 : 2 * gridSizeX *gridSizeY *gridSizeZ])                             \
        present(pDst[0 : 2 * imageSizeX *imageSizeY *imageSizeZ])
  for (int dZ_dst = 0; dZ_dst < imageSizeZ; dZ_dst++) {
    for (int dX_dst = 0; dX_dst < imageSizeX; dX_dst++) {
      for (int dY_dst = 0; dY_dst < imageSizeY; dY_dst++) {

        offsetY = (int)(((T1)gridSizeY / (T1)2.0) - ((T1)imageSizeY / (T1)2.0));
        offsetX = (int)(((T1)gridSizeX / (T1)2.0) - ((T1)imageSizeX / (T1)2.0));
        offsetZ = (int)(((T1)gridSizeZ / (T1)2.0) - ((T1)imageSizeZ / (T1)2.0));

        dY_src = dY_dst + offsetY;
        dX_src = dX_dst + offsetX;
        dZ_src = dZ_dst + offsetZ;

        common_index_dst =
            dZ_dst * imageSizeX * imageSizeY + dX_dst * imageSizeY + dY_dst;
        common_index_src =
            dZ_src * gridSizeX * gridSizeY + dX_src * gridSizeY + dY_src;

        // dst[common_index_dst] = src[common_index_src];

        // dst[common_index_dst].real() = src[common_index_src].real();
        // dst[common_index_dst].imag() = src[common_index_src].imag();
        pDst[2 * common_index_dst] = pSrc[2 * common_index_src];
        pDst[2 * common_index_dst + 1] = pSrc[2 * common_index_src + 1];
      }
    }
  }
}

// We need to oversample the FFT for gridding to work so we need to zero pad
template <typename T1>
void zero_pad2d(T1 *__restrict pDst, T1 *__restrict pSrc, int imageSizeX,
                int imageSizeY, T1 gridOS) {
  // (gridSizeX,gridSizeY) is the size of 'src'

  int offsetY;
  int offsetX;
  int dX_dst;
  int dY_dst;
  int common_index_dst;
  int common_index_src;

  offsetY =
      (int)(((T1)imageSizeY * gridOS / (T1)2.0) - ((T1)imageSizeY / (T1)2.0));
  offsetX =
      (int)(((T1)imageSizeX * gridOS / (T1)2.0) - ((T1)imageSizeX / (T1)2.0));
  int destSize = imageSizeX * gridOS * imageSizeY * gridOS;

#pragma acc kernels present(pDst[0 : 2 * destSize])                            \
                                present(pSrc[0 : 2 * imageSizeX *imageSizeY])
  {
#pragma acc loop
    for (int jj = 0; jj < 2 * destSize; jj++) {
      pDst[jj] = 0.0;
    }

#pragma acc loop collapse(2) independent
    for (int dY_src = 0; dY_src < imageSizeY; dY_src++) {
      for (int dX_src = 0; dX_src < imageSizeX; dX_src++) {

        dY_dst = dY_src + offsetY;
        dX_dst = dX_src + offsetX;

        common_index_dst = dY_dst * imageSizeX * gridOS + dX_dst;
        common_index_src = dY_src * imageSizeX + dX_src;

        // dst[common_index_dst] = src[common_index_src];
        pDst[2 * common_index_dst] = (pSrc[2 * common_index_src]); // Real
        pDst[2 * common_index_dst + 1] =
            (pSrc[2 * common_index_src + 1]); // Imaginary
      }
    }
  }
}

// We need to oversample the FFT for gridding to work so we need to zero pad
template <typename T1>
void zero_pad3d(T1 *__restrict pDst, T1 *__restrict pSrc, int imageSizeX,
                int imageSizeY, int imageSizeZ, T1 gridOS) {
  // (gridSizeX,gridSizeY,gridSizeZ) is the size of 'src'

  int offsetY;
  int offsetX;
  int offsetZ;

  int dY_dst;
  int dX_dst;
  int dZ_dst;
  int common_index_dst;
  int common_index_src;

  offsetY =
      (int)(((T1)imageSizeY * gridOS / (T1)2.0) - ((T1)imageSizeY / (T1)2.0));
  offsetX =
      (int)(((T1)imageSizeX * gridOS / (T1)2.0) - ((T1)imageSizeX / (T1)2.0));
  offsetZ =
      (int)(((T1)imageSizeZ * gridOS / (T1)2.0) - ((T1)imageSizeZ / (T1)2.0));
  int destSize =
      imageSizeX * gridOS * imageSizeY * gridOS * imageSizeZ * gridOS;

#pragma acc kernels present(pDst[0 : 2 * destSize]) present(                   \
    pSrc[0 : 2 * imageSizeX *imageSizeY *imageSizeZ])
  {
#pragma acc loop
    for (int jj = 0; jj < 2 * destSize; jj++) {
      pDst[jj] = 0.0;
    }

#pragma acc loop collapse(3) independent
    for (int dZ_src = 0; dZ_src < imageSizeZ; dZ_src++) {
      for (int dY_src = 0; dY_src < imageSizeY; dY_src++) {
        for (int dX_src = 0; dX_src < imageSizeX; dX_src++) {
          dY_dst = dY_src + offsetY;
          dX_dst = dX_src + offsetX;
          dZ_dst = dZ_src + offsetZ;

          common_index_dst =
              dZ_dst * imageSizeX * imageSizeY * gridOS * gridOS +
              dX_dst * imageSizeY * gridOS + dY_dst;
          common_index_src =
              dZ_src * imageSizeX * imageSizeY + dX_src * imageSizeY + dY_src;

          // dst[common_index_dst] = src[common_index_src];
          pDst[2 * common_index_dst] = (pSrc[2 * common_index_src]); // Real
          pDst[2 * common_index_dst + 1] =
              (pSrc[2 * common_index_src + 1]); // Imaginary
        }
      }
    }
  }
}

// Circshift routines for fftshift.
// Are the X and Y references correct?
template <typename T>
void circshift2(T *__restrict pDst, const T *__restrict pSrc, int xdim,
                int ydim, int xshift, int yshift) {
  int ii, jj;

#pragma acc parallel loop collapse(2) independent present(                     \
    pSrc[0 : 2 * xdim *ydim]) present(pDst[0 : 2 * xdim *ydim])
  for (int x = 0; x < xdim; x++) {
    ii = (x + xshift) % xdim;
    for (int y = 0; y < ydim; y++) {
      jj = (y + yshift) % ydim;
      // out[ii+jj*xdim] = in[x+y*xdim];
      pDst[2 * (ii + jj * xdim)] = pSrc[2 * (x + y * xdim)];
      pDst[2 * (ii + jj * xdim) + 1] = pSrc[2 * (x + y * xdim) + 1];
    }
  }
}

template <typename T>
void circshift3(T *__restrict pDst, const T *__restrict pSrc, int xdim,
                int ydim, int zdim, int xshift, int yshift, int zshift) {
  // cout << "Entering circshift3 " << endl;
  int ii, jj, kk;

#pragma acc parallel loop collapse(3) independent present(                     \
    pSrc[0 : 2 * xdim *ydim *zdim]) present(pDst[0 : 2 * xdim *ydim *zdim])
  for (int x = 0; x < xdim; x++) {
    ii = (x + xshift) % xdim;
    for (int y = 0; y < ydim; y++) {
      jj = (y + yshift) % ydim;
      for (int z = 0; z < zdim; z++) {
        kk = (z + zshift) % zdim;
        // out[jj+ii*ydim+kk*xdim*ydim] = in[y+x*ydim+z*xdim*ydim];
        pDst[2 * (jj + ii * ydim + kk * xdim * ydim)] =
            pSrc[2 * (y + x * ydim + z * xdim * ydim)];
        pDst[2 * (jj + ii * ydim + kk * xdim * ydim) + 1] =
            pSrc[2 * (y + x * ydim + z * xdim * ydim) + 1];
      }
    }
  }
}

// We need to fftshift to move the DC point to the center of the array. We could
// use (-1)^(i+j+k)
// but then we have a problem with the fftshift and ifftshift for odd matrix
// sizes
template <typename T>
void fftshift2(T *__restrict out, const T *__restrict in, int xdim, int ydim) {
  circshift2<T>(out, in, xdim, ydim, std::floor(xdim / 2.0),
                std::floor(ydim / 2.0));
}

template <typename T>
void ifftshift2(T *__restrict out, const T *__restrict in, int xdim, int ydim) {
  circshift2<T>(out, in, xdim, ydim, std::ceil(xdim / 2.0),
                std::ceil(ydim / 2.0));
}

template <typename T>
void fftshift3(T *__restrict out, const T *__restrict in, int xdim, int ydim,
               int zdim) {
  circshift3<T>(out, in, xdim, ydim, zdim, std::floor(xdim / 2.0),
                std::floor(ydim / 2.0), std::floor(zdim / 2.0));
}

template <typename T>
void ifftshift3(T *__restrict out, const T *__restrict in, int xdim, int ydim,
                int zdim) {
  circshift3<T>(out, in, xdim, ydim, zdim, std::ceil(xdim / 2.0),
                std::ceil(ydim / 2.0), std::floor(zdim / 2.0));
}

template <typename T1>
void fft2shift_grid(std::complex<T1> *__restrict src, int dimY, int dimX) {
  //(dimX,dimY) is the size of 'src'
  int common_index_dst;
#pragma acc kernels loop
  for (int dY_dst = 0; dY_dst < dimY; dY_dst++) {
    for (int dX_dst = 0; dX_dst < dimX; dX_dst++) {
      common_index_dst = dY_dst * dimX + dX_dst;
      src[common_index_dst] *= std::pow(-1.0, (T1)(dY_dst + dX_dst));
    }
  }
}

template <typename T1>
void fft3shift_grid(std::complex<T1> *__restrict src, int dimY, int dimX,
                    int dimZ) {
  //(dimX,dimY,dimZ) is the size of 'src'
  int common_index_dst;

#pragma acc kernels loop
  for (int dY_dst = 0; dY_dst < dimY; dY_dst++) {
    for (int dX_dst = 0; dX_dst < dimX; dX_dst++) {
      for (int dZ_dst = 0; dZ_dst < dimZ; dZ_dst++) {
        common_index_dst = dY_dst * dimX * dimZ + dX_dst * dimZ + dZ_dst;
        src[common_index_dst] *= std::pow(-1.0, (T1)(dY_dst + dX_dst + dZ_dst));
      }
    }
  }
}
// Explicit Instantiation

template float bessi0<float>(float x);
template double bessi0<double>(double x);

template void calculateLUT<float>(float beta, float width, float *&LUT,
                                  uword &sizeLUT);
template void calculateLUT<double>(double beta, double width, double *&LUT,
                                   uword &sizeLUT);
/*
template float kernel_value_LUT(float dist, const float *LUT, uword sizeLUT,
                                float width);
template double kernel_value_LUT(double dist, const double *LUT, uword sizeLUT,
                                 double width);
*/
template void deinterleave_data2d<float>(float *__restrict pSrc,
                                         float *__restrict outR_d,
                                         float *__restrict outI_d, int imageX,
                                         int imageY);
template void deinterleave_data2d<double>(double *__restrict pSrc,
                                          double *__restrict outR_d,
                                          double *__restrict outI_d, int imageX,
                                          int imageY);
template void deinterleave_data3d<float>(float *__restrict pSrc,
                                         float *__restrict outR_d,
                                         float *__restrict outI_d, int imageX,
                                         int imageY, int imageZ);
template void deinterleave_data3d<double>(double *__restrict pSrc,
                                          double *__restrict outR_d,
                                          double *__restrict outI_d, int imageX,
                                          int imageY, int imageZ);
template void deapodization2d<float>(float *__restrict pDst,
                                     float *__restrict pSrc, int imageX,
                                     int imageY, float kernelWidth, float beta,
                                     float gridOS);
template void deapodization2d<double>(double *__restrict pDst,
                                      double *__restrict pSrc, int imageX,
                                      int imageY, double kernelWidth,
                                      double beta, double gridOS);
template void deapodization3d<float>(float *__restrict pDst,
                                     float *__restrict pSrc, int imageX,
                                     int imageY, int imageZ, float kernelWidth,
                                     float beta, float gridOS);
template void deapodization3d<double>(double *__restrict pDst,
                                      double *__restrict pSrc, int imageX,
                                      int imageY, int imageZ,
                                      double kernelWidth, double beta,
                                      double gridOS);

template void crop_center_region2d<float>(float *__restrict pDst,
                                          float *__restrict pSrc,
                                          int imageSizeX, int imageSizeY,
                                          int gridSizeX, int gridSizeY);
template void crop_center_region2d<double>(double *__restrict pDst,
                                           double *__restrict pSrc,
                                           int imageSizeX, int imageSizeY,
                                           int gridSizeX, int gridSizeY);

template void crop_center_region3d<float>(float *__restrict pDst,
                                          float *__restrict pSrc,
                                          int imageSizeX, int imageSizeY,
                                          int imageSizeZ, int gridSizeX,
                                          int gridSizeY, int gridSizeZ);
template void crop_center_region3d<double>(double *__restrict pDst,
                                           double *__restrict pSrc,
                                           int imageSizeX, int imageSizeY,
                                           int imageSizeZ, int gridSizeX,
                                           int gridSizeY, int gridSizeZ);

template void zero_pad2d<float>(float *__restrict pDst, float *__restrict pSrc,
                                int imageSizeX, int imageSizeY, float gridOS);
template void zero_pad2d<double>(double *__restrict pDst,
                                 double *__restrict pSrc, int imageSizeX,
                                 int imageSizeY, double gridOS);
template void zero_pad3d<float>(float *__restrict pDst, float *__restrict pSrc,
                                int imageSizeX, int imageSizeY, int imageSizeZ,
                                float gridOS);
template void zero_pad3d<double>(double *__restrict pDst,
                                 double *__restrict pSrc, int imageSizeX,
                                 int imageSizeY, int imageSizeZ, double gridOS);
template void circshift2<float>(float *__restrict pDst,
                                const float *__restrict pSrc, int xdim,
                                int ydim, int xshift, int yshift);
template void circshift2<double>(double *__restrict pDst,
                                 const double *__restrict pSrc, int xdim,
                                 int ydim, int xshift, int yshift);
template void fftshift2<float>(float *__restrict out,
                               const float *__restrict in, int xdim, int ydim);
template void fftshift2<double>(double *__restrict out,
                                const double *__restrict in, int xdim,
                                int ydim);
template void ifftshift2<float>(float *__restrict out,
                                const float *__restrict in, int xdim, int ydim);
template void ifftshift2<double>(double *__restrict out,
                                 const double *__restrict in, int xdim,
                                 int ydim);
template void fftshift3<float>(float *__restrict out,
                               const float *__restrict in, int xdim, int ydim,
                               int zdim);
template void fftshift3<double>(double *__restrict out,
                                const double *__restrict in, int xdim, int ydim,
                                int zdim);
template void ifftshift3<float>(float *__restrict out,
                                const float *__restrict in, int xdim, int ydim,
                                int zdim);
template void ifftshift3<double>(double *__restrict out,
                                 const double *__restrict in, int xdim,
                                 int ydim, int zdim);
template void fft2shift_grid<float>(std::complex<float> *__restrict src,
                                    int dimY, int dimX);
template void fft2shift_grid<double>(std::complex<double> *__restrict src,
                                     int dimY, int dimX);
template void fft3shift_grid<float>(std::complex<float> *__restrict src,
                                    int dimY, int dimX, int dimZ);
template void fft3shift_grid<double>(std::complex<double> *__restrict src,
                                     int dimY, int dimX, int dimZ);
