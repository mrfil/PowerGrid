/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [directRecon.h]

    Synopsis    [Code for calculating sum-of-squares images from fully samled 
                    non-cartesian data]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [11/10/2019]

 *****************************************************************************/

#ifndef POWERGRID_DIRECTRECON_H
#define POWERGRID_DIRECTRECON_H

#include "../PowerGrid.h"
#include "../processIsmrmrd.hpp"
using namespace arma;

template<typename T>
Col<T> calc3DDensityCompensation(Col<T> kx, Col<T> ky, Col<T> kz, uword Ninplane, uword Nz);

template<typename T>
Col<complex<T>> gridCoilImages(uword Ninplane, uword Nz, ISMRMRD::Dataset *d, ISMRMRD::IsmrmrdHeader *hdr, 
                            acqTracking *acqTrack, uword NPhase, uword NEcho, uword NAvg, uword NRep);
template<typename T>
Col<T> calcSumOfSquaresImage(Col<complex<T>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword NCoils);

extern template Col<float> calc3DDensityCompensation<float>(Col<float>, Col<float>, Col<float>, uword, uword);
extern template Col<double> calc3DDensityCompensation<double>(Col<double>, Col<double>, Col<double>, uword, uword);

extern template Col<complex<float>> gridCoilImages<float>( uword, uword, ISMRMRD::Dataset*, ISMRMRD::IsmrmrdHeader*, acqTracking*,  uword,  uword,  uword,  uword);
extern template Col<complex<double>> gridCoilImages<double>( uword, uword, ISMRMRD::Dataset*, ISMRMRD::IsmrmrdHeader*, acqTracking*,  uword,  uword,  uword,  uword);

extern template Col<float> calcSumOfSquaresImage<float>(Col<complex<float>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword NCoils);
extern template Col<double> calcSumOfSquaresImage<double>(Col<complex<double>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword NCoils);
#endif // DIRECTRECON

