/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [processIsmrmrd.hpp]

    Synopsis    [Helper functions for working with ISMRMRD files.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef POWERGRID_PROCESSISMRMRD_HPP
#define POWERGRID_PROCESSISMRMRD_HPP

#include "PowerGrid.h"

using namespace arma;

void openISMRMRDData(std::string inputDataFile, ISMRMRD::Dataset *&d, ISMRMRD::IsmrmrdHeader &hdr, acqTracking *&acqTrack) {
    std::cout << "trying to create an ISMRMD::Dataset object" << std::endl;
    d = new ISMRMRD::Dataset(inputDataFile.c_str(), "dataset", false);
    //std::cout << "address of the  ISMRMD::Dataset object = " << d << std::endl;
    std::string xml;
    std::cout << "trying to read the header from the ISMRMD::Dataset object" << std::endl;
    d->readHeader(xml);
    std::cout << "read the header from the ISMRMD::Dataset object" << std::endl;
    //ISMRMRD::IsmrmrdHeader hdr;
    ISMRMRD::deserialize(xml.c_str(), hdr);
    std::cout << "trying to deserialze the xml header from the string" << std::endl;

	//Intitalize the acqTracking object to handle bookkeeping for the file
	acqTrack = new acqTracking(d,hdr);

}

void closeISMRMRDData( ISMRMRD::Dataset *&d, ISMRMRD::IsmrmrdHeader &hdr, acqTracking *&acqTrack) {

	//delete hdr;
	delete acqTrack;
	delete d;

}


//Write conversion from Image format to Armadillo matrix format for further use.
template<typename T1>
arma::Col<T1> convertFromNDArrayToArma(ISMRMRD::NDArray<T1> &inArray) {
    std::cout << "Converting NDArray to Arma Column" << std::endl;
    arma::Col<T1> temp;
    arma::uword numElems = inArray.getNumberOfElements();
    std::cout << "Num Elements to be converted = " << numElems << std::endl;
    const T1 *auxData = NULL;
    auxData = inArray.getDataPtr();
    temp = arma::Col<T1>(auxData, numElems);

    return temp;
}

template<typename T1>
arma::Col<T1> getISMRMRDFieldMap(ISMRMRD::Dataset *d) {
    const std::string fieldMap = "FieldMap";
    arma::Col<double> FM_temp;
    arma::Col<T1> FM;
    if (d->getNumberOfNDArrays(fieldMap) > 1) {
        //Throw error here
    }
    ISMRMRD::NDArray<double> tempArray;
    d->readNDArray(fieldMap, 0, tempArray);
    FM_temp = convertFromNDArrayToArma<double>(tempArray);
    FM = conv_to<arma::Col<T1>>::from(FM_temp);
    return FM;
}

template<typename T1>
arma::Col<T1> getISMRMRDSenseMap(ISMRMRD::Dataset *d) {
    const std::string senseMap = "SENSEMap";
    arma::Col<std::complex<double>> sen_temp;
    arma::Col<T1> sen;
    std::cout << "About to get number of NDArrays (SENSEMap)" << std::endl;
    if (d->getNumberOfNDArrays(senseMap) > 1) {
        //Throw error here
        std::cout << "OH NO!!!! SEGV APPROACHING!" << std::endl;
    }
    std::cout << "Got number of NDArrays (SENSEMap)" << std::endl;
    ISMRMRD::NDArray<std::complex<double>> tempArray;
    d->readNDArray(senseMap, 0, tempArray);
    sen_temp = convertFromNDArrayToArma(tempArray);
    sen = conv_to<arma::Col<T1>>::from(sen_temp);
    return sen;
}


template<typename T1>
arma::Col<T1> getISMRMRDPhaseMaps(ISMRMRD::Dataset *d) {
	const std::string phaseMaps = "PhaseMaps";
	arma::Col<double> pMaps_temp;
	arma::Col<T1> pMaps;
	std::cout << "About to get number of NDArrays (PhaseMaps)" << std::endl;
	if (d->getNumberOfNDArrays(phaseMaps) > 1) {
		//Throw error here
		std::cout << "OH NO!!!! SEGV APPROACHING!" << std::endl;
	}
	std::cout << "Got number of NDArrays (PhaseMaps)" << std::endl;
	ISMRMRD::NDArray<double> tempArray;
	d->readNDArray(phaseMaps, 0, tempArray);
	pMaps_temp = convertFromNDArrayToArma(tempArray);
	pMaps = conv_to<arma::Col<T1>>::from(pMaps_temp);
	return pMaps;
}


template<typename T1>
arma::Col<T1> getISMRMRDCompletePhaseMap(ISMRMRD::Dataset *d, uword NSlice, uword NSet, uword NRep, uword NAvg, uword NPhase, uword NEcho, uword NSeg, uword imageSize)
{
	arma::Col<T1> pMaps = getISMRMRDPhaseMaps<T1>(d);

	std::string xml;
	std::cout << "trying to read the header from the ISMRMD::Dataset object" << std::endl;
	d->readHeader(xml);
	std::cout << "read the header from the ISMRMD::Dataset object" << std::endl;
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(), hdr);

	uword NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum;
	uword NSetMax   = hdr.encoding[0].encodingLimits.set->maximum;
	uword NRepMax   = hdr.encoding[0].encodingLimits.repetition->maximum;
	uword NAvgMax   = hdr.encoding[0].encodingLimits.average->maximum;
	uword NSegMax   = hdr.encoding[0].encodingLimits.segment->maximum;
	uword NEchoMax  = hdr.encoding[0].encodingLimits.contrast->maximum;
	uword NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum;

	uword NShotMax  = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum;
	uword NParMax   = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum;

	uword PMapSize   = imageSize*(NShotMax+1)*(NParMax+1);
	uword startIndex = PMapSize*NSlice + PMapSize*NSliceMax*NAvg + PMapSize*NSliceMax*NAvgMax*NPhase + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEcho + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEchoMax*NRep + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEchoMax*NRepMax*NSeg;

	std::cout << "PMap slicing startIndex = " << startIndex << std::endl;

	uword endIndex = PMapSize*NSlice + PMapSize*NSliceMax*NAvg + PMapSize*NSliceMax*NAvgMax*NPhase + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEcho + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEchoMax*NRep + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEchoMax*NRepMax*NSeg + PMapSize - 1;

	std::cout << "PMap slicing endIndex = " << endIndex << std::endl;

	arma::Col<T1> pMapOut = pMaps.subvec(startIndex, endIndex);
	return pMapOut;
}

template<typename T1>
arma::Col<T1> getISMRMRDCompleteSENSEMap(ISMRMRD::Dataset *d, uword NSlice, uword imageSize)
{
	arma::Col<T1> SENSEMaps = getISMRMRDSenseMap<T1>(d);

	std::string xml;
	std::cout << "trying to read the header from the ISMRMD::Dataset object" << std::endl;
	d->readHeader(xml);
	std::cout << "read the header from the ISMRMD::Dataset object" << std::endl;
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(), hdr);


	// Grab first acquisition to get parameters (We assume all subsequent
	// acquisitions will be similar).
	ISMRMRD::Acquisition acq;
	d->readAcquisition(0, acq);
	uword nCoils = acq.active_channels();

	uword NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum;

	uword SENSEMapSize   = imageSize*nCoils;
	uword startIndex = SENSEMapSize*NSlice;

	std::cout << "SENSEMap slicing startIndex = " << startIndex << std::endl;

	uword endIndex = SENSEMapSize*(NSlice+1) - 1;

	std::cout << "SENSEMap slicing endIndex = " << endIndex << std::endl;
	std::cout << "SENSEMap number of rows = " << SENSEMaps.n_rows << std::endl;
	arma::Col<T1> SENSEMapOut = SENSEMaps.subvec(startIndex, endIndex);

	return SENSEMapOut;
}

template<typename T1>
arma::Col<T1> getISMRMRDCompleteFieldMap(ISMRMRD::Dataset *d, uword NSlice, uword imageSize)
{
	arma::Col<T1> FieldMaps = getISMRMRDFieldMap<T1>(d);

	std::string xml;
	std::cout << "trying to read the header from the ISMRMD::Dataset object" << std::endl;
	d->readHeader(xml);
	std::cout << "read the header from the ISMRMD::Dataset object" << std::endl;
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(), hdr);


	// Grab first acquisition to get parameters (We assume all subsequent
	// acquisitions will be similar).
	ISMRMRD::Acquisition acq;
	d->readAcquisition(0, acq);
	uword nCoils = acq.active_channels();

	uword NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum;

	uword startIndex = imageSize*NSlice;

	std::cout << "FieldMap slicing startIndex = " << startIndex << std::endl;

	uword endIndex = imageSize*(NSlice+1) - 1;

	std::cout << "FieldMap slicing endIndex = " << endIndex << std::endl;

	arma::Col<T1> FieldMapOut = FieldMaps.subvec(startIndex, endIndex);
	return FieldMapOut;
}

template<typename T1>
void processISMRMRDInput(std::string inputDataFile, ISMRMRD::Dataset *&d, ISMRMRD::IsmrmrdHeader &hdr,
                         arma::Col<T1> &FM, arma::Col<std::complex<T1>> &sen, acqTracking *&acqTrack) {
    std::cout << "About to open ISMRMRD file for input" << std::endl;
	openISMRMRDData(inputDataFile, d, hdr, acqTrack);
    std::cout << "Opened ISMRMRD file for input " << std::endl;
    std::cout << "About to get the Field map" << std::endl;
    FM = getISMRMRDFieldMap<T1>(d);
    std::cout << "About to get the SENSE map" << std::endl;
    sen = getISMRMRDSenseMap<std::complex<T1>>(d);

    return;
}

template<typename T1>
void getISMRMRDAcqData(ISMRMRD::Dataset *d, uword Nacq, Col<std::complex<T1>> &data, Col<T1> &kx, Col<T1> &ky,
                       Col<T1> &kz, Col<T1> &tvec) {
    ISMRMRD::Acquisition acq;
    d->readAcquisition(Nacq, acq);
    uword nro = acq.number_of_samples();
    uword nc = acq.active_channels();

    Mat<std::complex<T1>> dataTemp(nro, nc);
    for (uword jj = 0; jj < nc; jj++) {
        for (uword kk = 0; kk < nro; kk++) {
            dataTemp(kk, jj) = acq.data(kk, jj);
        }
    }

    //Deal with trajectories
    for (uword ii = 0; ii < acq.number_of_samples(); ii++) {
        kx(ii) = acq.traj(0, ii);
        ky(ii) = acq.traj(1, ii);
        kz(ii) = acq.traj(2, ii);
        tvec(ii) = acq.traj(3, ii);
    }
    data = vectorise(dataTemp);
    return;
}
/*
template<typename T1>
ISMRMRD::Acquisition getISMRMRDAcq(ISMRMRD::Dataset *d, uword Nacq) {
    ISMRMRD::Acquisition acq;
    d->readAcquisition(Nacq, acq);
    uword nro = acq.number_of_samples();
    uword nc = acq.active_channels();

    Mat<std::complex<T1>> dataTemp(nro, nc);
    for (uword jj = 0; jj < nc; jj++) {
        for (uword kk = 0; kk < nro; kk++) {
            dataTemp(kk, jj) = acq.data(kk, jj);
        }
    }

    //Deal with trajectories
    for (uword ii = 0; ii < acq.number_of_samples(); ii++) {
        kx(ii) = acq.traj(0, ii);
        ky(ii) = acq.traj(1, ii);
        kz(ii) = acq.traj(2, ii);
        tvec(ii) = acq.traj(3, ii);
    }
    data = vectorise(dataTemp);
    return acq;
}
*/
template<typename T1>
void writeISMRMRDImageData(ISMRMRD::Dataset *d, Col<std::complex<T1>> &image, uword Nx, uword Ny, uword Nz) {
    ISMRMRD::Image<std::complex<T1>> img_out(Nx, Ny, Nz, 1);
    savemat("testImage.mat", "img", image);
    for (int ii = 0; ii < Ny; ii++) {
        for (int jj = 0; jj < Nx; jj++) {
            for (int kk = 0; kk < Nz; kk++) {
                img_out(jj, ii, kk) = image(ii + jj * Ny + kk * Ny * Nx);
            }
        }
    }

    //Let's set some header details
    img_out.setImageType(ISMRMRD::ISMRMRD_IMTYPE_COMPLEX);

    //Write out the image
    d->appendImage("image", img_out);

}

template<typename T1>
void getCompleteISMRMRDAcqData(ISMRMRD::Dataset *d, acqTracking *acqTrack, uword NSlice, uword NRep, uword NAve, uword NEcho, uword NPhase, Col<std::complex<T1>> &data,
                               Col<T1> &kx, Col<T1> &ky, Col<T1> &kz, Col<T1> &tvec)
{

	//Initialization
	Mat<std::complex<T1>> dataTemp, acqWork;

	Mat<T1> kxTemp, kyTemp, kzTemp, tvecTemp;
	Mat<T1> kxWork, kyWork, kzWork, tvecWork;
	uword numAcq = d->getNumberOfAcquisitions();
	bool firstData = true;
	ISMRMRD::Acquisition acq;
	std::cout << "Num of acquisitions in dataset = " << numAcq << std::endl;
	int acqIndx = -1;
	//
	//for (uword acqIndx = 0; acqIndx < numAcq; acqIndx++) {
	for (uword NPar = 0; NPar < acqTrack->NParMax; NPar++) {
		for (uword NShot = 0; NShot < acqTrack->NShotMax; NShot++) {

			acqIndx = acqTrack->acqArray(NShot, NPar, NSlice, NRep, NAve, NEcho, NPhase);
			if (acqIndx != -1) {
				d->readAcquisition(acqIndx, acq);
				uword nro = acq.number_of_samples();
				uword nc = acq.active_channels();
				ISMRMRD::EncodingCounters encIdx = acq.idx();

				std::cout << "Grabbing acq index #" << acqIndx << std::endl;
				acqWork.zeros(nro, nc);
				kxWork.zeros(nro,1);
				kyWork.zeros(nro,1);
				kzWork.zeros(nro,1);
				tvecWork.zeros(nro,1);

				for (uword jj = 0; jj<nc; jj++) {
					for (uword kk = 0; kk<nro; kk++) {
						acqWork(kk, jj) = static_cast<std::complex<T1>>(acq.data(kk, jj));
					}
				}

				std::cout << "Get Number of Trajectory entries = " << acq.getNumberOfTrajElements() << std::endl;

				int NSlice = encIdx.slice;
				int NRep   = encIdx.repetition;
				int NAvg   = encIdx.average;
				int NEcho  = encIdx.contrast;
				int NPhase = encIdx.phase;


				std::cout << "NPar = "   << NPar   << std::endl;
				std::cout << "NShot = "  << NShot  << std::endl;
				std::cout << "NSlice = " << NSlice << std::endl;
				std::cout << "NRep = "   << NRep   << std::endl;
				std::cout << "NAvg = "   << NAvg   << std::endl;
				std::cout << "NEcho = "  << NEcho  << std::endl;
				std::cout << "NPhase = " << NPhase << std::endl;

				//Deal with trajectories
				for (uword ii = 0; ii<nro; ii++) {
					kxWork(ii)   = static_cast<T1>(acq.traj(0, ii));
					kyWork(ii)   = static_cast<T1>(acq.traj(1, ii));
					kzWork(ii)   = static_cast<T1>(acq.traj(2, ii));
					//std::cout << "kxWork(" << ii << ")= " << kxWork(ii) << std::endl;
					tvecWork(ii) = static_cast<T1>(acq.traj(3, ii));
				}

				//Append Data points to vectors
				if (firstData) {
					firstData = false; //Sentinel
					dataTemp = acqWork;
					kxTemp   = kxWork;
					kyTemp   = kyWork;
					kzTemp   = kzWork;
					tvecTemp = tvecWork;
				} else {
					std::cout << "Concatenating arrays" << std::endl;
					dataTemp = join_cols(dataTemp, acqWork);
					kxTemp   = join_cols(kxTemp, kxWork).eval();
					for (int test = 0; test < kxTemp.n_rows; test++) {
						std::cout << "kxTemp(" << test << ") = " << kxTemp(test,0) << std::endl;
					}
					kyTemp   = join_cols(kyTemp, kyWork).eval();
					kzTemp   = join_cols(kzTemp, kzWork).eval();
					tvecTemp = join_cols(tvecTemp, tvecWork).eval();
				}

			}
		}
	}


    //Vectorise coils from matrix to column vector
    data = vectorise(dataTemp);
	kx   = vectorise(kxTemp);
	ky   = vectorise(kyTemp);
	kz   = vectorise(kzTemp);
	tvec = vectorise(tvecTemp);

	kx.save("kx.dat", raw_ascii);
	ky.save("ky.dat", raw_ascii);
	kz.save("kz.dat", raw_ascii);
	tvec.save("tvec.dat", raw_ascii);
    //savemat("kxOut.mat", "kx", kx);
    //savemat("kyOut.mat", "ky", ky);
    //savemat("kzOut.mat", "kz", kz);
    //savemat("tvecOut.mat", "tvec", tvec);
    savemat("dataOut.mat", "dataOut", data);

    return;
}

#endif //POWERGRID_PROCESSISMRMRD_HPP
