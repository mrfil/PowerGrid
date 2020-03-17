/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [reconSolve.h]

    Synopsis    [Helper functions to setup iterative image reconstructions.]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]
                [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/17/2016]

*****************************************************************************/

#ifndef POWERGRID_RECONSOLVE_H
#define POWERGRID_RECONSOLVE_H

#include "PGIncludes.h"
#include "Gdft.h"
#include "GdftR2.h"
#include "Gnufft.h"
#include "SENSE.h"
#include "pcSENSE.h"
#include "TimeSegmentation.h"
#include "QuadPenalty.h"
#include "TVPenalty.h"
#include "solve_pwls_pcg.hpp"

//#ifdef PowerGridMPI
//#include "MPI/mpipcSENSE.h"
//#endif

using namespace arma;


template <typename T1>
void initImageSpaceCoords(Col<T1> &ix, Col<T1> &iy, Col<T1> &iz, uword Nx,
                          uword Ny, uword Nz);

/*
// Template parameters are  T1: data precision (double, float, FP16 etc...),
// TObj: Transform Object, RObj is regularization object
template <typename T1, typename TObj, typename RObj>
Col<complex<T1> > reconSolve(Col<complex<T1> > data, TObj &Sg, RObj R,
                             Col<T1> kx, Col<T1> ky, Col<T1> kz, uword Nx,
                             uword Ny, uword Nz, Col<T1> tvec, uword niter);
*/

// Template parameters are  T1: data precision (double, float, FP16 etc...),
// TObj: Transform Object, RObj is regularization object
template <typename T1, typename TObj, typename RObj>
Col<complex<T1>> reconSolve(Col<complex<T1>> data, TObj& Sg, RObj R,
                            Col<T1> kx, Col<T1> ky, Col<T1> kz, uword Nx,
                            uword Ny, uword Nz, Col<T1> tvec, uword niter) {
  typedef std::complex<T1> CxT1;

  // Col<T1> ix, iy, iz;
  // initImageSpaceCoords(ix,iy,iz,Nz,Ny,Nz);

  // Data weighting term - use ones unless we have something better to use
  Col<T1> W;
  W.ones(data.n_rows);
  Col<std::complex<T1>> xinit;
  xinit.zeros(Nx * Ny * Nz);

  Col<CxT1> imageOut;
  imageOut = solve_pwls_pcg<T1, TObj, RObj>(xinit, Sg, W, data, R, niter);

  return imageOut;
}

// Explicit Instantiation
extern template void initImageSpaceCoords<float>(Col<float> &, Col<float> &,
                                                 Col<float> &, uword Nx, uword Ny,
                                                 uword Nz);
extern template void initImageSpaceCoords<double>(Col<double> &, Col<double> &,
                                                  Col<double> &, uword Nx, uword Ny,
                                                  uword Nz);
// Explicit Instantiation
extern template
Col<complex<float>>
reconSolve(Col<complex<float>>, SENSE<float, Gnufft<float>>&, QuadPenalty<float>, Col<float>, Col<float>, Col<float>, uword,
              uword, uword, Col<float>, uword);

extern template
Col<complex<float> > reconSolve<float, SENSE<float, Gdft<float> >, QuadPenalty<float>>(Col<complex<float>>, SENSE<float, Gdft<float>>&,
                                                                                           QuadPenalty<float>, Col<float>, Col<float>,
                                                                                           Col<float>, uword, uword, uword, Col<float>,
                                                                                           uword);

                                                                                           extern template
Col<complex<float> > reconSolve<float, SENSE<float, GdftR2<float> >, QuadPenalty<float>>(Col<complex<float>>, SENSE<float, GdftR2<float>>&,
                                                                                           QuadPenalty<float>, Col<float>, Col<float>,
                                                                                           Col<float>, uword, uword, uword, Col<float>,
                                                                                           uword);

extern template
Col<complex<float>>
reconSolve(Col<complex<float>>, pcSENSE<float>&, QuadPenalty<float>, Col<float>, Col<float>, Col<float>, uword,
		uword, uword, Col<float>, uword);


extern template  
Col<complex<float>> reconSolve(Col<complex<float>>, LRobj<float, pcSenseTimeSeg<float>>&,
		R_lowRank<float, QuadPenalty<float>>, Col<float>, Col<float>,
		Col<float>, uword, uword, uword, Col<float>,
		uword);

extern template
Col<complex<double>>
reconSolve(Col<complex<double>>, SENSE<double, Gnufft<double>>&, QuadPenalty<double>, Col<double>, Col<double>, Col<double>, uword,
		uword, uword, Col<double>, uword);

extern template
Col<complex<double> > reconSolve<double, SENSE<double, Gdft<double> >, QuadPenalty<double>>(Col<complex<double>>, SENSE<double, Gdft<double>>&,
		QuadPenalty<double>, Col<double>, Col<double>,
		Col<double>, uword, uword, uword, Col<double>,
		uword);

extern template
Col<complex<double> > reconSolve<double, SENSE<double, GdftR2<double> >, QuadPenalty<double>>(Col<complex<double>>, SENSE<double, GdftR2<double>>&,
		QuadPenalty<double>, Col<double>, Col<double>,
		Col<double>, uword, uword, uword, Col<double>,
		uword);

extern template
Col<complex<double>>
reconSolve(Col<complex<double>>, pcSENSE<double>&, QuadPenalty<double>, Col<double>, Col<double>, Col<double>, uword,
		uword, uword, Col<double>, uword);

/*
extern template
Col<complex<double> >
reconSolve<double, SENSE<double, Gnufft<double> > &,
           QuadPenalty<double> &>(Col<complex<double> >, SENSE<double, Gnufft<double> > &,
                                  QuadPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
                                  uword, uword, Col<double>, uword);
extern template
Col<complex<double> >
reconSolve<double, SENSE<double, Gdft<double> > &, QuadPenalty<double> &>(Col<complex<double> >, SENSE<double, Gdft<double> > &,
                                                                          QuadPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
                                                                          uword, uword, Col<double>, uword);
extern template
Col<complex<float> >
reconSolve<float, SENSE<float, Gnufft<float> > &, TVPenalty<float> &>(Col<complex<float> >, SENSE<float, Gnufft<float> > &,
                                                                      TVPenalty<float> &, Col<float>, Col<float>, Col<float>, uword, uword,
                                                                      uword, Col<float>, uword);
extern template
Col<complex<float> > reconSolve<float, SENSE<float, Gdft<float> > &, TVPenalty<float> &>(Col<complex<float> >, SENSE<float, Gdft<float> > &,
                                                                                         TVPenalty<float> &, Col<float>, Col<float>,
                                                                                         Col<float>, uword, uword, uword, Col<float>,
                                                                                         uword);
extern template
Col<complex<double> >
reconSolve<double, SENSE<double, Gnufft<double> > &, TVPenalty<double> &>(Col<complex<double> >, SENSE<double, Gnufft<double> > &,
                                                                          TVPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
                                                                          uword, uword, Col<double>, uword);
extern template
Col<complex<double> >
reconSolve<double, SENSE<double, Gdft<double> > &, TVPenalty<double> &>(Col<complex<double> >, SENSE<double, Gdft<double> > &,
                                                                        TVPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
                                                                        uword, uword, Col<double>, uword);


*/
#endif // POWERGRID_RECONSOLVE_H
