/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [griddingSupport.h]

    Synopsis    [Support functions for the forward and adjoint non-uniform
                    Fast Fourier Transform (NUFFT) on CPU and GPU via
                    OpenACC.]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/2/2016]

*****************************************************************************/

#ifdef _OPENACC
#include "accelmath.h"
#include "openacc.h"
#define COS(a) cos(a)
#define SIN(a) sin(a)
#define SINH(a) sinh(a)
#else // On CPU
#define COS(a) std::cos(a)
#define SIN(a) std::sin(a)
#define SINH(a) std::sinh(a)
#endif

#include "PGIncludes.h"

using namespace std; // where to put?

// From Numerical Recipes in C, 2nd Edition
// Just a vanilla I(0,x) function approximation
#pragma acc routine seq
template <typename T1> T1 bessi0(T1 x);

template <typename T1>
void calculateLUT(T1 beta, T1 width, T1 *&LUT, uword &sizeLUT);

#pragma acc routine seq
template <typename T1>
T1 kernel_value_LUT(T1 dist, const T1 *LUT, uword sizeLUT, T1 width);

template <typename T1>
void deinterleave_data2d(T1 *__restrict pSrc, T1 *__restrict outR_d,
                         T1 *__restrict outI_d, int imageX, int imageY);

template <typename T1>
void deinterleave_data3d(T1 *__restrict pSrc, T1 *__restrict outR_d,
                         T1 *__restrict outI_d, int imageX, int imageY,
                         int imageZ);

template <typename T1>
void deapodization2d(T1 *__restrict pDst, T1 *__restrict pSrc, int imageX,
                     int imageY, T1 kernelWidth, T1 beta, T1 gridOS);

template <typename T1>
void deapodization3d(T1 *__restrict pDst, T1 *__restrict pSrc, int imageX,
                     int imageY, int imageZ, T1 kernelWidth, T1 beta,
                     T1 gridOS);

template <typename T1>
void crop_center_region2d(T1 *__restrict pDst, T1 *__restrict pSrc,
                          int imageSizeX, int imageSizeY, int gridSizeX,
                          int gridSizeY);

template <typename T1>
void crop_center_region3d(T1 *__restrict pDst, T1 *__restrict pSrc,
                          int imageSizeX, int imageSizeY, int imageSizeZ,
                          int gridSizeX, int gridSizeY, int gridSizeZ);

template <typename T1>
void zero_pad2d(T1 *__restrict pDst, T1 *__restrict pSrc, int imageSizeX,
                int imageSizeY, T1 gridOS);

template <typename T1>
void zero_pad3d(T1 *__restrict pDst, T1 *__restrict pSrc, int imageSizeX,
                int imageSizeY, int imageSizeZ, T1 gridOS);

template <typename T>
void circshift2(T *__restrict pDst, const T *__restrict pSrc, int xdim,
                int ydim, int xshift, int yshift);

template <typename T>
void fftshift2(T *__restrict out, const T *__restrict in, int xdim, int ydim);

template <typename T>
void fftshift3(T *__restrict out, const T *__restrict in, int xdim, int ydim,
               int zdim);

template <typename T>
void ifftshift2(T *__restrict out, const T *__restrict in, int xdim, int ydim);

template <typename T>
void ifftshift3(T *__restrict out, const T *__restrict in, int xdim, int ydim,
                int zdim);

template <typename T1>
void fft2shift_grid(std::complex<T1> *__restrict src, int dimY, int dimX);

template <typename T1>
void fft3shift_grid(std::complex<T1> *__restrict src, int dimY, int dimX,
                    int dimZ);

template <typename T1>
void normalize_fft2d(T1 *__restrict pDst, T1 *__restrict pSrc, int gridSizeX,
		int gridSizeY);

template <typename T1>
void normalize_fft3d(T1 *__restrict pDst, T1 *__restrict pSrc,
		int gridSizeX, int gridSizeY, int gridSizeZ);


// Explicit Instantiation
extern template float bessi0<float>(float x);
extern template double bessi0<double>(double x);
extern template void calculateLUT<float>(float beta, float width, float *&LUT,
                                         uword &sizeLUT);
extern template void calculateLUT<double>(double beta, double width,
                                          double *&LUT, uword &sizeLUT);

extern template
float kernel_value_LUT(float dist, const float *LUT, uword sizeLUT, float width);
extern template
double kernel_value_LUT(double dist, const double *LUT, uword sizeLUT, double width);
extern template void deinterleave_data2d<float>(float *__restrict pSrc,
                                                float *__restrict outR_d,
                                                float *__restrict outI_d,
                                                int imageX, int imageY);
extern template void deinterleave_data2d<double>(double *__restrict pSrc,
                                                 double *__restrict outR_d,
                                                 double *__restrict outI_d,
                                                 int imageX, int imageY);
extern template void deinterleave_data3d<float>(float *__restrict pSrc,
                                                float *__restrict outR_d,
                                                float *__restrict outI_d,
                                                int imageX, int imageY,
                                                int imageZ);
extern template void deinterleave_data3d<double>(double *__restrict pSrc,
                                                 double *__restrict outR_d,
                                                 double *__restrict outI_d,
                                                 int imageX, int imageY,
                                                 int imageZ);
extern template void deapodization2d<float>(float *__restrict pDst,
                                            float *__restrict pSrc, int imageX,
                                            int imageY, float kernelWidth,
                                            float beta, float gridOS);
extern template void deapodization2d<double>(double *__restrict pDst,
                                             double *__restrict pSrc,
                                             int imageX, int imageY,
                                             double kernelWidth, double beta,
                                             double gridOS);
extern template void deapodization3d<float>(float *__restrict pDst,
                                            float *__restrict pSrc, int imageX,
                                            int imageY, int imageZ,
                                            float kernelWidth, float beta,
                                            float gridOS);
extern template void deapodization3d<double>(double *__restrict pDst,
                                             double *__restrict pSrc,
                                             int imageX, int imageY, int imageZ,
                                             double kernelWidth, double beta,
                                             double gridOS);

extern template void crop_center_region2d<float>(float *__restrict pDst,
                                                 float *__restrict pSrc,
                                                 int imageSizeX, int imageSizeY,
                                                 int gridSizeX, int gridSizeY);
extern template void crop_center_region2d<double>(double *__restrict pDst,
                                                  double *__restrict pSrc,
                                                  int imageSizeX,
                                                  int imageSizeY, int gridSizeX,
                                                  int gridSizeY);

extern template void crop_center_region3d<float>(float *__restrict pDst,
                                                 float *__restrict pSrc,
                                                 int imageSizeX, int imageSizeY,
                                                 int imageSizeZ, int gridSizeX,
                                                 int gridSizeY, int gridSizeZ);
extern template void
crop_center_region3d<double>(double *__restrict pDst, double *__restrict pSrc,
                             int imageSizeX, int imageSizeY, int imageSizeZ,
                             int gridSizeX, int gridSizeY, int gridSizeZ);

extern template void zero_pad2d<float>(float *__restrict pDst,
                                       float *__restrict pSrc, int imageSizeX,
                                       int imageSizeY, float gridOS);
extern template void zero_pad2d<double>(double *__restrict pDst,
                                        double *__restrict pSrc, int imageSizeX,
                                        int imageSizeY, double gridOS);
extern template void zero_pad3d<float>(float *__restrict pDst,
                                       float *__restrict pSrc, int imageSizeX,
                                       int imageSizeY, int imageSizeZ,
                                       float gridOS);
extern template void zero_pad3d<double>(double *__restrict pDst,
                                        double *__restrict pSrc, int imageSizeX,
                                        int imageSizeY, int imageSizeZ,
                                        double gridOS);
extern template void circshift2<float>(float *__restrict pDst,
                                       const float *__restrict pSrc, int xdim,
                                       int ydim, int xshift, int yshift);
extern template void circshift2<double>(double *__restrict pDst,
                                        const double *__restrict pSrc, int xdim,
                                        int ydim, int xshift, int yshift);
extern template void fftshift2<float>(float *__restrict out,
                                      const float *__restrict in, int xdim,
                                      int ydim);
extern template void fftshift2<double>(double *__restrict out,
                                       const double *__restrict in, int xdim,
                                       int ydim);
extern template void ifftshift2<float>(float *__restrict out,
                                       const float *__restrict in, int xdim,
                                       int ydim);
extern template void ifftshift2<double>(double *__restrict out,
                                        const double *__restrict in, int xdim,
                                        int ydim);
extern template void fftshift3<float>(float *__restrict out,
                                      const float *__restrict in, int xdim,
                                      int ydim, int zdim);
extern template void fftshift3<double>(double *__restrict out,
                                       const double *__restrict in, int xdim,
                                       int ydim, int zdim);
extern template void ifftshift3<float>(float *__restrict out,
                                       const float *__restrict in, int xdim,
                                       int ydim, int zdim);
extern template void ifftshift3<double>(double *__restrict out,
                                        const double *__restrict in, int xdim,
                                        int ydim, int zdim);
extern template void fft2shift_grid<float>(std::complex<float> *__restrict src,
                                           int dimY, int dimX);
extern template void
fft2shift_grid<double>(std::complex<double> *__restrict src, int dimY,
                       int dimX);
extern template void fft3shift_grid<float>(std::complex<float> *__restrict src,
                                           int dimY, int dimX, int dimZ);
extern template void
fft3shift_grid<double>(std::complex<double> *__restrict src, int dimY, int dimX,
                       int dimZ);

extern template void
normalize_fft2d<float>(float *__restrict pDst, float *__restrict pSrc, int gridSizeX,
		int gridSizeY);

extern template void
normalize_fft2d<double>(double *__restrict pDst, double *__restrict pSrc, int gridSizeX,
		int gridSizeY);

extern template void
normalize_fft3d<double>(double *__restrict pDst, double *__restrict pSrc,
		int gridSizeX, int gridSizeY, int gridSizeZ);

extern template void
normalize_fft3d<float>(float *__restrict pDst, float *__restrict pSrc,
		int gridSizeX, int gridSizeY, int gridSizeZ);