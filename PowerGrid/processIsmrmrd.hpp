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
typedef std::tuple<std::size_t,std::size_t,std::size_t> D3tuple;

void openISMRMRDData(std::string inputDataFile, ISMRMRD::Dataset *&d, ISMRMRD::IsmrmrdHeader &hdr, acqTracking *&acqTrack) {
	RANGE()
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
	RANGE()
    std::cout << "Converting NDArray to Arma Column" << std::endl;
    arma::Col<T1> temp;
    arma::uword numElems = inArray.getNumberOfElements();
    std::cout << "Num Elements to be converted = " << numElems << std::endl;
    T1 *auxData = NULL;
    auxData = inArray.getDataPtr();
    temp = arma::Col<T1>(auxData, numElems, false, false);

    return temp;
}

template<typename T1>
arma::Col<T1> getISMRMRDFieldMap(ISMRMRD::Dataset *d) {
	RANGE()
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
	RANGE()
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
	RANGE()
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
arma::Col<complex<T1>> getISMRMRDTemporalBasis(ISMRMRD::Dataset *d) {
	RANGE()
    const std::string tempBasis = "v";
    arma::Col<complex<double>> vBasis_temp;
    arma::Col<complex<T1>> vBasis;
    if (d->getNumberOfNDArrays(tempBasis) > 1) {
        //Throw error here
    }
    ISMRMRD::NDArray<complex<double>> tempArray;
    d->readNDArray(tempBasis, 0, tempArray);
    vBasis_temp = convertFromNDArrayToArma<complex<double>>(tempBasis);
    vBasis = conv_to<arma::Col<complex<T1>>>::from(vBasis_temp);
    return vBasis;
}

template<typename T1>
arma::Col<T1> getISMRMRDCompletePhaseMap(ISMRMRD::Dataset *d, uword NSlice, uword NSet, uword NRep, uword NAvg, uword NPhase, uword NEcho, uword NSeg, uword imageSize)
{
	RANGE()
	arma::Col<T1> pMaps = getISMRMRDPhaseMaps<T1>(d);

	std::string xml;
	std::cout << "trying to read the header from the ISMRMD::Dataset object" << std::endl;
	d->readHeader(xml);
	std::cout << "read the header from the ISMRMD::Dataset object" << std::endl;
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(), hdr);

	uword NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum + 1;
	uword NSetMax   = hdr.encoding[0].encodingLimits.set->maximum + 1;
	uword NRepMax   = hdr.encoding[0].encodingLimits.repetition->maximum + 1;
	uword NAvgMax   = hdr.encoding[0].encodingLimits.average->maximum + 1;
	uword NSegMax   = hdr.encoding[0].encodingLimits.segment->maximum + 1;
	uword NEchoMax  = hdr.encoding[0].encodingLimits.contrast->maximum  + 1;
	uword NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum  + 1;

	uword NShotMax  = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum + 1;
	uword NParMax   = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum + 1;

    std::cout << "NSliceMax = " << NSliceMax << std::endl;
    std::cout << "NSetMax = " << NSetMax << std::endl;
    std::cout << "NRepMax = " << NRepMax << std::endl;
    std::cout << "NAvgMax = " << NAvgMax << std::endl;
    std::cout << "NEchoMax = " << NEchoMax << std::endl;
    std::cout << "NPhaseMax = " << NPhaseMax << std::endl;
    std::cout << "NSegMax = " << NSegMax << std::endl;

	uword PMapSize   = imageSize*(NShotMax)*(NParMax);
	uword startIndex = PMapSize*NSlice + PMapSize*NSliceMax*NAvg + PMapSize*NSliceMax*NAvgMax*NPhase + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEcho + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEchoMax*NRep + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEchoMax*NRepMax*NSeg;

	std::cout << "PMap slicing startIndex = " << startIndex << std::endl;

	uword endIndex = PMapSize*NSlice + PMapSize*NSliceMax*NAvg + PMapSize*NSliceMax*NAvgMax*NPhase + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEcho + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEchoMax*NRep + PMapSize*NSliceMax*NAvgMax*NPhaseMax*NEchoMax*NRepMax*NSeg + PMapSize - 1;

	std::cout << "PMap slicing endIndex = " << endIndex << std::endl;
  std::cout << "pMaps length = " << pMaps.n_rows << std::endl;
  std::cout << "PMapSize length = " << PMapSize << std::endl;

	arma::Col<T1> pMapOut = pMaps.subvec(startIndex, endIndex);

	return pMapOut;
}

template<typename T1>
arma::Col<T1> getISMRMRDCompleteSENSEMap(ISMRMRD::Dataset *d, arma::Col<T1> &SENSEMaps, uword NSlice, uword imageSize)
{
	RANGE()
	//arma::Col<T1> SENSEMaps = getISMRMRDSenseMap<T1>(d);

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
arma::Col<T1> getISMRMRDCompleteFieldMap(ISMRMRD::Dataset *d, arma::Col<T1> &FieldMaps, uword NSlice, uword imageSize)
{
	RANGE()
	//arma::Col<T1> FieldMaps = getISMRMRDFieldMap<T1>(d);

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

	//FieldMapOut.save("FieldMap.dat", raw_ascii);

	return FieldMapOut;
}

template<typename T1>
void processISMRMRDInput(std::string inputDataFile, ISMRMRD::Dataset *&d, ISMRMRD::IsmrmrdHeader &hdr,
                         arma::Col<T1> &FM, arma::Col<std::complex<T1>> &sen, acqTracking *&acqTrack) {
	RANGE()
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
	RANGE()
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
	RANGE()
    ISMRMRD::Image<std::complex<T1>> img_out(Nx, Ny, Nz, 1);
	//image.save("testImage.dat", raw_ascii);
	Col<T1> realImage = real(image).eval();
	Col<T1> imagImage = imag(image).eval();
	//realImage.save("testImageReal.dat", raw_ascii);
	//imagImage.save("testImageImag.dat", raw_ascii);

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
	RANGE()
	//Initialization
	Mat<std::complex<T1>> acqWork;
  	Cube<std::complex<T1>> dataWork;
	Mat<T1> kxWork, kyWork, kzWork, tvecWork;
	uword numAcqTotal = d->getNumberOfAcquisitions();
	bool firstData = true;
	ISMRMRD::Acquisition acq;
	std::cout << "Num of acquisitions in dataset = " << numAcqTotal << std::endl;
	int acqIndx = -1;
	int numAcqs = 0;
  	int nro = -1, nc = -1;
	//for (uword acqIndx = 0; acqIndx < numAcq; acqIndx++) {
  	for (uword NPar = 0; NPar < acqTrack->NParMax; NPar++) {
		for (uword NShot = 0; NShot < acqTrack->NShotMax; NShot++) {

			acqIndx = acqTrack->acqArray(NShot, NPar, NSlice, NRep, NAve, NEcho, NPhase);
      if (acqIndx != -1) {
        numAcqs++;
        if(firstData) {
          d->readAcquisition(acqIndx, acq);
          nro = acq.number_of_samples();
          nc = acq.active_channels();
          //std::cout << "Nro = " << nro << std::endl;
          //std::cout << "Nc = " << nc << std::endl;
          firstData = false;
        }
      }
    }
  }

  // Preallocating storage for all of the data and trajectories.
  dataWork.zeros(nro,nc,numAcqs);
  kxWork.zeros(nro,numAcqs);
  kyWork.zeros(nro,numAcqs);
  kzWork.zeros(nro,numAcqs);
  tvecWork.zeros(nro,numAcqs);
  acqWork.zeros(nro,nc);
  uword curAcq = 0;
	for (uword NPar = 0; NPar < acqTrack->NParMax; NPar++) {
		for (uword NShot = 0; NShot < acqTrack->NShotMax; NShot++) {

			acqIndx = acqTrack->acqArray(NShot, NPar, NSlice, NRep, NAve, NEcho, NPhase);
			if (acqIndx != -1) {
				d->readAcquisition(acqIndx, acq);
				nro = acq.number_of_samples();
				nc = acq.active_channels();

				ISMRMRD::EncodingCounters encIdx = acq.idx();

				//std::cout << "Grabbing acq index #" << acqIndx << std::endl;

				for (uword jj = 0; jj<nc; jj++) {
					for (uword kk = 0; kk<nro; kk++) {
						acqWork(kk, jj) = static_cast<std::complex<T1>>(acq.data(kk, jj));
					}
				}
				/*
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
				*/
				//Deal with trajectories
				for (uword ii = 0; ii<nro; ii++) {
					kxWork(ii,curAcq)   = static_cast<T1>(acq.traj(0, ii));
					kyWork(ii,curAcq)   = static_cast<T1>(acq.traj(1, ii));
					kzWork(ii,curAcq)   = static_cast<T1>(acq.traj(2, ii));
					tvecWork(ii,curAcq) = static_cast<T1>(acq.traj(3, ii));
				}

				dataWork.slice(curAcq) = acqWork;

        curAcq++;


			}
		}
	}

  // Need to permute the dataWork.
  dataWork = permute(dataWork,D3tuple(1,3,2));

  //Vectorise coils from matrix to column vector
  data = vectorise(dataWork);
	kx   = vectorise(kxWork);
	ky   = vectorise(kyWork);
	kz   = vectorise(kzWork);
	tvec = vectorise(tvecWork);
  /*
	kx.save("kx.dat", raw_ascii);
	ky.save("ky.dat", raw_ascii);
	kz.save("kz.dat", raw_ascii);
	tvec.save("tvec.dat", raw_ascii);
  */

  return;
}

#endif //POWERGRID_PROCESSISMRMRD_HPP
