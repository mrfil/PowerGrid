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
#include "pcSenseTimeSeg.h"

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

template  Col<complex<float>> reconSolve(Col<complex<float>>, SENSE<float, GdftR2<float>>&,
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

template  Col<complex<float>> reconSolve(Col<complex<float>>, pcSenseTimeSeg<float>&,
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

template  Col<complex<double>> reconSolve(Col<complex<double>>, SENSE<double, GdftR2<double>>&,
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

template  Col<complex<double>> reconSolve(Col<complex<double>>, pcSenseTimeSeg<double>&,
		QuadPenalty<double>, Col<double>, Col<double>,
		Col<double>, uword, uword, uword, Col<double>,
		uword);


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
