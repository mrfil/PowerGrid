//
// Created by acerja2 on 10/22/17.
//

#ifndef POWERGRID_ACQTRACKING_H
#define POWERGRID_ACQTRACKING_H

#include "PowerGrid.h"

class acqTracking {

public:
	~acqTracking() {};
	acqTracking(ISMRMRD::Dataset *d, ISMRMRD::IsmrmrdHeader &hdr);

	int NShotMax, NParMax, NSliceMax, NRepMax, NAvgMax, NEchoMax, NPhaseMax;

	ISMRMRD::NDArray<int> acqArray;
private:
	ISMRMRD::Dataset *d;



};

#endif //POWERGRID_ACQTRACKING_H
