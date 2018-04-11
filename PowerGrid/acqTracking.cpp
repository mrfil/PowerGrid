//
// Created by acerja2 on 10/22/17.
//

#include "acqTracking.h"

acqTracking::acqTracking(ISMRMRD::Dataset *dataSet, ISMRMRD::IsmrmrdHeader &hdr) {

	std::cout << "Entering acqTracking constructor." << std::endl;
	d = dataSet;

	uword numAcq = d->getNumberOfAcquisitions();
	//bool firstData = true;
	ISMRMRD::Acquisition acq;


	// Get bound of acquisitions headers for initializing the matrix
	NShotMax  = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum + 1;
	NParMax   = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum + 1;
	NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum + 1;
	//int NSetMax   = hdr.encoding[0].encodingLimits.set->maximum  + 1;
	NRepMax   = hdr.encoding[0].encodingLimits.repetition->maximum + 1 ;
	NAvgMax   = hdr.encoding[0].encodingLimits.average->maximum + 1;
	//int NSegMax   = hdr.encoding[0].encodingLimits.segment->maximum + 1;
	NEchoMax  = hdr.encoding[0].encodingLimits.contrast->maximum + 1;
	NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum + 1;


	int NShot, NPar, NSlice, NRep, NAvg, NEcho, NPhase;



	// Need to initialize a N-dimensional matrix or array to track which acquisition index to go to for a specific encoding
	// Can't use the ISMRMRD NDArray because it is only 7 dimensions long and we have more dimensions than that.
	// Buuuuutttt, we're going to stick with it for now because I don't want to figure out another nice ND array solution.

	std::vector<size_t> dims;
	dims.push_back(NShotMax);
	dims.push_back(NParMax);
	dims.push_back(NSliceMax);
	dims.push_back(NRepMax);
	dims.push_back(NAvgMax);
	dims.push_back(NEchoMax);
	dims.push_back(NPhaseMax);
	// Notice that we are skipping NSetMax and NSegMax because we have chosen to avoid those indexes as a lab.
	// If ISMRMRD gives us more than 7 dimensions in NDArray we can put them in!

	this->acqArray.resize(dims);

	//Set the entire thing to -1 to signify a nonexistant acquisition. Can't use zero because acquisitions are zero indexed.
	memset(this->acqArray.getDataPtr(),-1,sizeof(int)*NShotMax*NParMax*NSliceMax*NRepMax*NAvgMax*NEchoMax*NPhaseMax);


	// Now we can sort the acquisitions
	for (uword acqIndx = 0; acqIndx < numAcq; acqIndx++) {
		std::cout << "Scanning through acquisition # " << acqIndx << std::endl;
		// Scanning through the file.
		d->readAcquisition(acqIndx, acq);
		//uword nro = acq.number_of_samples();
		//uword nc = acq.active_channels();
		ISMRMRD::EncodingCounters encIdx = acq.idx();

		NShot = encIdx.kspace_encode_step_1;
		NPar = encIdx.kspace_encode_step_2;
		NSlice = encIdx.slice;
		NRep = encIdx.repetition;
		NAvg = encIdx.average;
		NEcho = encIdx.contrast;
		NPhase = encIdx.phase;



		std::cout << "Writing back to the NDArray" << std::endl;
		std::cout << "NPar = "   << NPar << std::endl;
		std::cout << "NShot = "  << NShot << std::endl;
		std::cout << "NSlice = " << NSlice << std::endl;
		std::cout << "NRep = "   << NRep << std::endl;
		std::cout << "NAvg = "   << NAvg << std::endl;
		std::cout << "NEcho = "  << NEcho << std::endl;
		std::cout << "NPhase = " << NPhase << std::endl;

		// Sort the indexes as we need to.
		this->acqArray(NShot,NPar,NSlice,NRep,NAvg,NEcho,NPhase) = acqIndx;
	}


	return;

}
