/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [directRecon.cpp]

    Synopsis    [Code for calculating sum-of-squares images from fully samled 
                    non-cartesian data]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [11/10/2019]

 *****************************************************************************/
#include "directRecon.h"

template<typename T>
Col<T> calc3DDensityCompensation(Col<T> kx, Col<T> ky, Col<T> kz, uword Ninplane, uword Nz) {
    
    uword Niter = 2;
    Col<T> ix, iy, iz;
    uword overSampling = 3.0;
    initImageSpaceCoords(ix, iy, iz, Ninplane, Ninplane, Nz);
    Gnufft<T> G(kx.n_rows, (T) overSampling, Ninplane, Ninplane, Nz, kx, ky, kz, ix, iy, iz);

    Col<complex<T>> dcf(kx.n_rows);
    dcf.ones();
    dcf /= ((T)dcf.n_rows * (T)dcf.n_rows);
    Col<complex<T>> C1(overSampling * Ninplane * overSampling * Ninplane * overSampling * Nz);
    C1.ones();
    
    //Col<T> C1abs(overSampling * Ninplane * overSampling * Ninplane * overSampling * Nz);
    //C1abs.ones();
    
    Col<complex<T>> C2(kx.n_rows);
    C2.ones();

    for(int ii = 0; ii < Niter; ii++) {

        C1 = G.adjointSpatialInterp(dcf);
        C2 = G.forwardSpatialInterp(C1);
        dcf = dcf / real(C2);

    }

    dcf = normalise(dcf,2);

    return abs(dcf);

}

template<typename T>
Col<T> calcSumOfSquaresImage(Col<complex<T>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword Ncoils)
{


    // Let's do this via armadillo
    Mat<complex<T>> inputData = reshape(coilImages, Nx * Nx * Nz * NSliceMax, Ncoils);

    Col<T> imSOS = sqrt(sum(abs(inputData),1)); // Remember indices are uniformly zero indexed in armadillo/C++)

    return imSOS;
}

template<typename T>
Col<complex<T>> gridCoilImages(uword Ninplane, uword Nz, ISMRMRD::Dataset *d, ISMRMRD::IsmrmrdHeader *hdr, 
                            acqTracking *acqTrack, uword NPhase, uword NEcho, uword NAvg, uword NRep) 
{

    // Grab first acquisition to get parameters (We assume all subsequent
    // acquisitions will be similar).
    ISMRMRD::Acquisition acq;
    d->readAcquisition(0, acq);
    //uword nro = acq.number_of_samples();
    
    uword Ncoils = acq.active_channels();

	//uword NSet = 0; //Set is only used for arrayed ADCs
	//uword NSeg = 0;

    int NSlices = hdr->encoding[0].encodingLimits.slice->maximum + 1;

    Col<std::complex<T>> ImageTemp(Ninplane * Ninplane * Nz * NSlices * Ncoils);
	Col<T> kx, ky, kz, tvec;
    Col<std::complex<T>> data; 
    
    Col<T> ix, iy, iz;
    initImageSpaceCoords(ix, iy, iz, Ninplane, Ninplane, Nz);

    getCompleteISMRMRDAcqData<T>(d, acqTrack, 0, 0, 0, 0, 0, data, kx, ky,
	    kz, tvec);
    
    Col<T> dcf;
    dcf = calc3DDensityCompensation(kx, ky, kz, Ninplane, Nz);
    dcf.save("dcf.dat", raw_ascii);

    uword nro = dcf.n_rows;
    
    for (uword NSlice = 0; NSlice < NSlices; NSlice++) {
        Col<T> kx(nro), ky(nro), kz(nro), tvec(nro);
        Col<std::complex<T>> data(nro * Ncoils);

	    getCompleteISMRMRDAcqData<T>(d, acqTrack, NSlice, NRep, NAvg, NEcho, NPhase, data, kx, ky,
		    kz, tvec);

	    std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
	    std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
	    std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
	    std::cout << "Number of rows in data = " << data.n_rows << std::endl;
	    std::cout << "Number of columns in data = " << data.n_cols << std::endl;
                   
	    Gnufft<T> G(kx.n_rows, (T) 2.0, Ninplane, Ninplane, Nz, kx, ky, kz, ix, iy, iz);

        Col<complex<T>> tempData;
        Col<complex<T>> colImageTemp;
        for (int ii = 0; ii < Ncoils; ii++) {
            tempData = data.subvec(nro * ii, nro * (ii + 1) - 1);
            colImageTemp = G / ( dcf % tempData);
            ImageTemp.subvec((Ninplane * Ninplane * Nz * NSlices) * ii + (Ninplane * Ninplane * Nz * NSlice) , (Ninplane * Ninplane * Nz * NSlices) * ii + (Ninplane * Ninplane * Nz * (NSlice + 1)) - 1) 
                            =  colImageTemp;
        }     
    }
    return ImageTemp;
    //writeNiftiMagPhsImage<float>(filename,ImageTemp,Nx,Ny,Nz,NSlices,Ncoils);

}

template Col<float> calc3DDensityCompensation<float>(Col<float>, Col<float>, Col<float>,  uword,  uword);
template Col<double> calc3DDensityCompensation<double>(Col<double>, Col<double>, Col<double>,  uword,  uword);

template Col<complex<float>> gridCoilImages<float>( uword, uword, ISMRMRD::Dataset*, ISMRMRD::IsmrmrdHeader*, acqTracking*, uword, uword, uword, uword);
template Col<complex<double>> gridCoilImages<double>( uword, uword, ISMRMRD::Dataset*, ISMRMRD::IsmrmrdHeader*, acqTracking*, uword, uword, uword, uword);

template Col<float> calcSumOfSquaresImage<float>(Col<complex<float>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword NCoils);
template Col<double> calcSumOfSquaresImage<double>(Col<complex<double>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword NCoils);


