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

using namespace arma;
using namespace PowerGrid;

template <typename T1>
void initImageSpaceCoords(Col<T1> &ix, Col<T1> &iy, Col<T1> &iz, uword Nx,
                          uword Ny, uword Nz);

// Template parameters are  T1: data precision (double, float, FP16 etc...),
// TObj: Transform Object, RObj is regularization object
template <typename T1, typename TObj, typename RObj>
Col<complex<T1>> reconSolve(Col<complex<T1>> data, TObj &Sg, RObj &R,
                            Col<T1> kx, Col<T1> ky, Col<T1> kz, uword Nx,
                            uword Ny, uword Nz, Col<T1> tvec, uword niter);

// Explicit Instantiation
template void initImageSpaceCoords<float>(Col<float> &, Col<float> &,
                                          Col<float> &, uword Nx, uword Ny,
                                          uword Nz);
template void initImageSpaceCoords<double>(Col<double> &, Col<double> &,
                                           Col<double> &, uword Nx, uword Ny,
                                           uword Nz);

template <float, SENSE<float, Gnufft<float>>, QuadPenalty<float>>
reconSolve(Col<complex<float>>, SENSE<float, Gnufft<float>> &,
           QuadPenalty<float> &, Col<float>, Col<float>, Col<float>, uword,
           uword, uword, Col<float>, uword);
template <float, SENSE<float, Gdft<float>>, QuadPenalty<float>>
reconSolve(Col<complex<float>>, SENSE<float, Gdft<float>> &,
           QuadPenalty<float> &, Col<float>, Col<float>, Col<float>, uword,
           uword, uword, Col<float>, uword);
template <double, SENSE<double, Gnufft<double>>, QuadPenalty<double>>
reconSolve(Col<complex<double>>, SENSE<double, Gnufft<double>> &,
           QuadPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
           uword, uword, Col<double>, uword);
template <double, SENSE<double, Gdft<double>>, QuadPenalty<double>>
reconSolve(Col<complex<double>>, SENSE<double, Gdft<double>> &,
           QuadPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
           uword, uword, Col<double>, uword);
template <float, SENSE<float, Gnufft<float>>, TVPenalty<float>>
reconSolve(Col<complex<float>>, SENSE<float, Gnufft<float>> &,
           TVPenalty<float> &, Col<float>, Col<float>, Col<float>, uword, uword,
           uword, Col<float>, uword);
template <float, SENSE<float, Gdft<float>>, TVPenalty<float>>
reconSolve(Col<complex<float>>, SENSE<float, Gdft<float>> &, TVPenalty<float> &,
           Col<float>, Col<float>, Col<float>, uword, uword, uword, Col<float>,
           uword);
template <double, SENSE<double, Gnufft<double>>, TVPenalty<double>>
reconSolve(Col<complex<double>>, SENSE<double, Gnufft<double>> &,
           TVPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
           uword, uword, Col<double>, uword);
template <double, SENSE<double, Gdft<double>>, TVPenalty<double>>
reconSolve(Col<complex<double>>, SENSE<double, Gdft<double>> &,
           TVPenalty<double> &, Col<double>, Col<double>, Col<double>, uword,
           uword, uword, Col<double>, uword);

#endif // POWERGRID_RECONSOLVE_H
