/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [reconSolve.cpp]

    Synopsis    [Helper functions to setup iterative image reconstructions.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#include "reconSolve.h"
#include "pcSENSE.h"

template <typename T1>
void initImageSpaceCoords(Col<T1> &ix, Col<T1> &iy, Col<T1> &iz, uword Nx,
                          uword Ny, uword Nz) {
  // generate the image space coordinates of the voxels we want to reconstruct
  // after vectorizing ix and iy the image coordinates must match the Field and
  // SENSE map image coordinates
  Cube<T1> ixTemp(Nx, Ny, Nz);
  Cube<T1> iyTemp(Nx, Ny, Nz);
  Cube<T1> izTemp(Nx, Ny, Nz);

  for (uword ii = 0; ii < Ny; ii++) {     // y
    for (uword jj = 0; jj < Nx; jj++) {   // x
      for (uword kk = 0; kk < Nz; kk++) { // z

        ixTemp(ii, jj, kk) = ((T1)jj - (T1)Nx / 2.0) / ((T1)Nx);
        iyTemp(ii, jj, kk) = ((T1)ii - (T1)Ny / 2.0) / ((T1)Ny);
        izTemp(ii, jj, kk) = ((T1)kk - (T1)Nz / 2.0) / ((T1)Nz);
      }
    }
  }

  ix = vectorise(ixTemp);
  iy = vectorise(iyTemp);
  iz = vectorise(izTemp);
}

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
template void initImageSpaceCoords<float>(Col<float> &, Col<float> &,
                                          Col<float> &, uword Nx, uword Ny,
                                          uword Nz);
template void initImageSpaceCoords<double>(Col<double> &, Col<double> &,
                                           Col<double> &, uword Nx, uword Ny,
                                           uword Nz);
/*
template <float, SENSE<float, Gnufft<float>>, QuadPenalty<float>>
Col<complex<float>>
reconSolve(Col<complex<float>>, SENSE<float, Gnufft<float>>,
           QuadPenalty<float>, Col<float>, Col<float>, Col<float>, uword,
           uword, uword, Col<float>, uword);
*/

template
Col<complex<float>>
reconSolve(Col<complex<float>>, SENSE<float, Gnufft<float>>&,
		QuadPenalty<float>, Col<float>, Col<float>, Col<float>, uword,
		uword, uword, Col<float>, uword);

template  Col<complex<float>> reconSolve(Col<complex<float>>, SENSE<float, Gdft<float>>&,
                               QuadPenalty<float>, Col<float>, Col<float>,
                               Col<float>, uword, uword, uword, Col<float>,
                               uword);

template  Col<complex<float>> reconSolve(Col<complex<float>>, SENSE<float, TimeSegmentation<float,Gnufft<float>>>&,
		QuadPenalty<float>, Col<float>, Col<float>,
		Col<float>, uword, uword, uword, Col<float>,
		uword);

template  Col<complex<float>> reconSolve(Col<complex<float>>, pcSENSE<float>&,
		QuadPenalty<float>, Col<float>, Col<float>,
		Col<float>, uword, uword, uword, Col<float>,
		uword);

template
Col<complex<double>>
reconSolve(Col<complex<double>>, SENSE<double, Gnufft<double>>&,
		QuadPenalty<double>, Col<double>, Col<double>, Col<double>, uword,
		uword, uword, Col<double>, uword);

template  Col<complex<double>> reconSolve(Col<complex<double>>, SENSE<double, Gdft<double>>&,
		QuadPenalty<double>, Col<double>, Col<double>,
		Col<double>, uword, uword, uword, Col<double>,
		uword);

template  Col<complex<double>> reconSolve(Col<complex<double>>, SENSE<double, TimeSegmentation<double,Gnufft<double>>>&,
		QuadPenalty<double>, Col<double>, Col<double>,
		Col<double>, uword, uword, uword, Col<double>,
		uword);

template  Col<complex<double>> reconSolve(Col<complex<double>>, pcSENSE<double>&,
		QuadPenalty<double>, Col<double>, Col<double>,
		Col<double>, uword, uword, uword, Col<double>,
		uword);

#ifdef PowerGridMPI

template  Col<complex<float>> reconSolve(Col<complex<float>>, mpipcSENSE<float>&,
                               QuadPenalty<float>, Col<float>, Col<float>,
                               Col<float>, uword, uword, uword, Col<float>,
                               uword);
#endif
/*
template <double &, SENSE<double, Gnufft<double>> &, QuadPenalty<double> &> Col<complex<double>> reconSolve(Col<complex<double>>, SENSE<double, Gnufft<double>> &,
           QuadPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
           uword, uword, Col<double>, uword);
template Col<complex<double>> reconSolve<double, SENSE<double, Gdft<double>>, QuadPenalty<double>> (Col<complex<double>>, SENSE<double, Gdft<double>> &,
           QuadPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
           uword, uword, Col<double>, uword);
template <float &, SENSE<float, Gnufft<float>> &, TVPenalty<float> &> Col<complex<float>> reconSolve(Col<complex<float>>, SENSE<float, Gnufft<float>> &,
           TVPenalty<float> &, Col<float>, Col<float>, Col<float>, uword, uword,
           uword, Col<float>, uword);
template <float &, SENSE<float, Gdft<float>> &, TVPenalty<float> &> Col<complex<float>> reconSolve(Col<complex<float>>, SENSE<float, Gdft<float>> &,
                               TVPenalty<float> &, Col<float>, Col<float>,
                               Col<float>, uword, uword, uword, Col<float>,
                               uword);
template <double &, SENSE<double, Gnufft<double>> &, TVPenalty<double> &> Col<complex<double>> reconSolve(Col<complex<double>>, SENSE<double, Gnufft<double>> &,
           TVPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
           uword, uword, Col<double>, uword);
template <double &, SENSE<double, Gdft<double>> &, TVPenalty<double> &>  Col<complex<double>> reconSolve(Col<complex<double>>, SENSE<double, Gdft<double>> &,
           TVPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
           uword, uword, Col<double>, uword);
*/