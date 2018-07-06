/*
(C) Copyright 2010-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [ftCpuWithGrads.h]

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

#ifndef FTWithGrads_CPU_H
#define FTWithGrads_CPU_H

/*---------------------------------------------------------------------------*/
/*  Included library headers                                                 */
/*---------------------------------------------------------------------------*/
#include "PGIncludes.h"
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
    #include "openacc.h"
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

template <typename T1> T1 sinc(T1 x);

/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU kernel of the Fourier Transformation (FT).]             */
/*                                                                           */
/*  Description []                                                           */
/*                                                                           */
/*===========================================================================*/
template <typename T1>
void ftCpuWithGrads(T1 *kdata_r, T1 *kdata_i, const T1 *idata_r, const T1 *idata_i,
           const T1 *kx, const T1 *ky, const T1 *kz, const T1 *ix, const T1 *iy,
           const T1 *iz, const T1 *FM, const T1 *Gx, const T1 *Gy, const T1 *Gz,
           const T1 *t, const int num_k,
           const int num_i, const int num_x, const int num_y, const int num_z);

// Explicit Instantiations
extern template void ftCpuWithGrads<float>(float *, float *, const float *,
                                  const float *, const float *, const float *,
                                  const float *, const float *, const float *,
                                  const float *, const float *, const float *,
                                  const float *, const float *, const float *,
                                  const int, const int, const int, const int,
                                  const int);
extern template void ftCpuWithGrads<double>(double *, double *, const double *,
                                   const double *, const double *,
                                   const double *, const double *,
                                   const double *, const double *,
                                   const double *, const double *,
                                   const double *, 
                                   const double *, const double *,
                                   const double *, const int, const int,
                                   const int, const int, const int);
/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU kernel of the Inverse Fourier Transformation (IFT).] */
/*                                                                           */
/*  Description [] */
/*                                                                           */
/*===========================================================================*/
template <typename T1>
void iftCpuWithGrads(T1 *idata_r, T1 *idata_i, const T1 *kdata_r, const T1 *kdata_i,
            const T1 *kx, const T1 *ky, const T1 *kz, const T1 *ix,
            const T1 *iy, const T1 *iz, const T1 *FM, const T1 *Gx, 
            const T1 *Gy, const T1 *Gz, const T1 *t,
            const int num_k, const int num_i, const int num_x, const int num_y,
            const int num_z);

// Explicit Instantiations
extern template void iftCpuWithGrads<float>(float *, float *, const float *,
                                   const float *, const float *, const float *,
                                   const float *, const float *, const float *,
                                   const float *, const float *, const float *,
                                   const float *, const float *, const float *,
                                   const int, const int, const int, const int,
                                   const int);
extern template void iftCpuWithGrads<double>(double *, double *, const double *,
                                    const double *, const double *,
                                    const double *, const double *,
                                    const double *, const double *,
                                    const double *, const double *,
                                    const double *,
                                    const double *, const double *,
                                    const double *, const int, const int,
                                    const int, const int, const int);
/*---------------------------------------------------------------------------*/
/*  Namespace declared - end                                                 */
/*---------------------------------------------------------------------------*/

//}
//}

#endif // FTWithGrads_CPU_H
