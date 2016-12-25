/*
(C) Copyright 2010-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [ftCpu.h]

    Synopsis    [The CPU and OpenACC annotated version of the discrete Fourier
                transform and inverse discrete Fourier transform.]

    Description []

    Revision    [0.1; Initial build; Yue Zhuo, BIOE UIUC]
    Revision    [0.1.1; Add OpenMP, Code cleaning; Xiao-Long Wu, ECE UIUC]
    Revision    [1.0a; Further optimization, Code cleaning, Adding more
 comments;
                 Xiao-Long Wu, ECE UIUC, Jiading Gai, Beckman Institute]
    Revision    [1.1; Remove OpenMP and add OpenACC annotations for GPU
                  acceleration]
    Date        [4/19/2016]

 *****************************************************************************/

#ifndef FT_CPU_H
#define FT_CPU_H

/*---------------------------------------------------------------------------*/
/*  Included library headers                                                 */
/*---------------------------------------------------------------------------*/
// Numeric constants according to the precision type.
#ifdef ENABLE_DOUBLE_PRECISION
#define MRI_PI 3.1415926535897932384626433832795029
#define MRI_NN 64
#define MRI_DELTAZ 0.003
#define MRI_ZERO 0.0
#define MRI_ONE 1.0
#define MRI_NEG_ONE -1.0
#define MRI_POINT_FIVE 0.5
#define MRI_SMOOTH_FACTOR 0.0000001
#else
#define MRI_PI 3.1415926535897932384626433832795029f
#define MRI_NN 64
#define MRI_DELTAZ 0.003f
#define MRI_ZERO 0.0f
#define MRI_ONE 1.0f
#define MRI_NEG_ONE -1.0f
#define MRI_POINT_FIVE 0.5f
#define MRI_SMOOTH_FACTOR 0.000001f
#endif
/*---------------------------------------------------------------------------*/
/*  Namespace declared - begin                                               */
/*---------------------------------------------------------------------------*/

// namespace uiuc_mri {

/*---------------------------------------------------------------------------*/
/*  Function prototypes                                                      */
/*---------------------------------------------------------------------------*/
/*
    void
ftCpu(T1 *kdata_r, T1 *kdata_i,
      const T1 *idata_r, const T1 *idata_i,
      const DataTraj *ktraj, const DataTraj *itraj,
      const T1 *fm, const T1 *t,
      const int num_k, const int num_i
      );

    void
iftCpu(T1 *idata_r, T1 *idata_i,
       const T1 *kdata_r, const T1 *kdata_i,
       const DataTraj *ktraj, const DataTraj *itraj,
       const T1 *fm, const T1 *t,
       const int num_k, const int num_i
       );
*/
/*---------------------------------------------------------------------------*/
/*  Included library headers                                                 */
/*---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENACC
#include "accelmath.h"
#include "openacc.h"
#define COS(a) std::cos(a)
#define SIN(a) std::sin(a)
#define SQRT(a) std::sqrt(a)
#else
#define COS(a) std::cos(a)
#define SIN(a) std::sin(a)
#define SQRT(a) std::sqrt(a)
#endif

//#include <tools.h>
//#include <structures.h>

/*---------------------------------------------------------------------------*/
/*  Namespace declared - begin                                               */
/*---------------------------------------------------------------------------*/

// namespace uiuc_mri {

/*---------------------------------------------------------------------------*/
/*  Function definitions                                                     */
/*---------------------------------------------------------------------------*/

/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU version of the sin function.]                           */
/*                                                                           */
/*  Description [This function is used to avoid additional computations when */
/*      the data values are too small.]                                      */
/*                                                                           */
/*===========================================================================*/

template <typename T1> T1 sinc_cpu(T1 x);

/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU kernel of the Fourier Transformation (FT).]             */
/*                                                                           */
/*  Description []                                                           */
/*                                                                           */
/*===========================================================================*/
template <typename T1>
void ftCpu(T1 *kdata_r, T1 *kdata_i, const T1 *idata_r, const T1 *idata_i,
           const T1 *kx, const T1 *ky, const T1 *kz, const T1 *ix, const T1 *iy,
           const T1 *iz, const T1 *FM, const T1 *t, const int num_k,
           const int num_i);

// Explicit Instantiations
template void ftCpu<float>(float *, float *, const float *, const float *,
                           const float *, const float *, const float *,
                           const float *, const float *, const float *,
                           const float *, const float *, const int, const int);
template void ftCpu<double>(double *, double *, const double *, const double *,
                            const double *, const double *, const double *,
                            const double *, const double *, const double *,
                            const double *, const double *, const int,
                            const int);
/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU kernel of the Inverse Fourier Transformation (IFT).] */
/*                                                                           */
/*  Description [] */
/*                                                                           */
/*===========================================================================*/
template <typename T1>
void iftCpu(T1 *idata_r, T1 *idata_i, const T1 *kdata_r, const T1 *kdata_i,
            const T1 *kx, const T1 *ky, const T1 *kz, const T1 *ix,
            const T1 *iy, const T1 *iz, const T1 *FM, const T1 *t,
            const int num_k, const int num_i);

// Explicit Instantiations
template void iftCpu<float>(float *, float *, const float *, const float *,
                            const float *, const float *, const float *,
                            const float *, const float *, const float *,
                            const float *, const float *, const int, const int);
template void iftCpu<double>(double *, double *, const double *, const double *,
                             const double *, const double *, const double *,
                             const double *, const double *, const double *,
                             const double *, const double *, const int,
                             const int);
/*---------------------------------------------------------------------------*/
/*  Namespace declared - end                                                 */
/*---------------------------------------------------------------------------*/

//}
//}

#endif // FT_CPU_H
