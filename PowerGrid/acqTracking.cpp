//
// Created by acerja2 on 10/22/17.
//

#include "acqTracking.h"

acqTracking::acqTracking(ISMRMRD::Dataset &dataSet, ISMRMRD::IsmrmrdHeader &hdr) {


	d = &dataSet;

	uword numAcq = d->getNumberOfAcquisitions();
	bool firstData = true;
	ISMRMRD::Acquisition acq;


	// Get bound of acquisitions headers for initializing the matrix
	int NShotMax  = hdr.encoding[0].encodingLimits.kspace_encoding_step_0->maximum;
	int NParMax   = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum;
	int NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum;
	int NSetMax   = hdr.encoding[0].encodingLimits.set->maximum;
	int NRepMax   = hdr.encoding[0].encodingLimits.repetition->maximum;
	int NAvgMax   = hdr.encoding[0].encodingLimits.average->maximum;
	int NSegMax   = hdr.encoding[0].encodingLimits.segment->maximum;
	int NEchoMax  = hdr.encoding[0].encodingLimits.contrast->maximum;
	int NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum;


	int NShot, NPar, NSlice, NRep, NAvg, NEcho, NPhase;



	// Need to initialize a N-dimensional matrix or array to track which acquisition index to go to for a specific encoding
	// Can't use the ISMRMRD NDArray because it is only 7 dimensions long and we have more dimensions than that.
	// Buuuuutttt, we're going to stick with it for now because I don't want to figure out another nice ND array solution.

	std::vector<size_t> dims;
	dims.push_back(NShotMax  + 1);
	dims.push_back(NParMax   + 1);
	dims.push_back(NSliceMax + 1);
	dims.push_back(NRepMax   + 1);
	dims.push_back(NAvgMax   + 1);
	dims.push_back(NEchoMax  + 1);
	dims.push_back(NPhaseMax + 1);
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
		uword nro = acq.number_of_samples();
		uword nc = acq.active_channels();
		ISMRMRD::EncodingCounters encIdx = acq.idx();

		// Sort the indexes as we need to.
		this->acqArray(NShot,NPar,NSlice,NRep,NAvg,NEcho,NPhase) = acqIndx;
	}


	return;

}