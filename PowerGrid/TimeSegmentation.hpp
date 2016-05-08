/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [TimeSegmentation.hpp]

    Synopsis    [Wrappers to the cuFFT library supporting single and double
    				precision for GPU accelerated FFTs]

    Description []

    Revision    [0.1.0; Giang-Chau Ngo, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

// This using field correction by time segmentation
// The data is corrected to time 0 with reference to the time vector passed
//
//

#ifndef PowerGrid_TimeSegmentation_hpp
#define PowerGrid_TimeSegmentation_hpp

// We are using two template types at the moment. One for the type of data to be processed (ie Col<cx_double>) and one for the type of G object (ie Gfft<Col<cx_double>>
// T1 is the data type for complex, T2 is the data type for real data
template<typename T1, typename Tobj>
class TimeSegmentation {
    typedef complex <T1> CxT1;
public:
	TimeSegmentation();

    //Class variables go here
    uword n1; //Data size
    uword n2; //Image size
    uword L; //number of time segments
    uword type; //type of time segmentation
    uword Nshots; //Number of shots, used to reduce complexity of calculating interpolator
    T1 tau;        //time segment length
    T1 T_min;   // minimum time in the time vector (i.e. TE for spiral out)
    Tobj* obj;
    Col <T1> fieldMap; //Field map (in radians per second)
    Col <T1> timeVec;  //timing vector of when each data point was collected relative to the echo time (in seconds)
    Mat <CxT1> AA;        //interpolator coefficients for the different time segments
    CxT1 i = CxT1(0., 1.);

    //Class constructor
	TimeSegmentation(Tobj &G, Col <T1> map_in, Col <T1> timeVec_in, uword a, uword b, uword c, uword interptype = 1,
					 uword shots = 1)
    {
		cout << "Entering Class constructor" << endl;
	    n1 = a; //Data size
	    n2 = b;//Image size
	    L = c; //number of time segments
	    type = interptype; // type of time segmentation performed
	    Nshots = shots; // number of shots
	    obj = &G;
	    fieldMap = map_in;
    	cout << "N1 = " << n1 << endl;
		cout << "N2 = " << n2 << endl;
		cout << "L = " << L << endl;


		AA.set_size(n1, L); //time segments weights
	    timeVec = timeVec_in;
	    T_min =timeVec.min();
	    T1 rangt = timeVec.max()-T_min;
		tau = (rangt + datum::eps) / (L - 1); // it was L-1 before
	    timeVec = timeVec-T_min;

	    uword NOneShot = n1/Nshots;
	    if (L==1) {
		    tau = 0;
		    AA.ones();
	    }
	    else {
			Mat <CxT1> tempAA(NOneShot, L);
		    if (type==1) {// Hanning interpolator
			    cout << "Hanning interpolation" << endl;
			    for (unsigned int ii = 0; ii<L; ii++) {
				    for (unsigned int jj = 0; jj<NOneShot; jj++) {
					    if ((std::abs(timeVec(jj)-((ii)*tau)))<=tau) {
						    tempAA(jj, ii) = 0.5+0.5*std::cos((datum::pi)*(timeVec(jj)-((ii)*tau))/tau);
					    }
					    else {
						    tempAA(jj, ii) = 0.0;
					    }
				    }
			    }
			    AA = repmat(tempAA, Nshots, 1);
		    }
		    else if (type==2) { // Min-max interpolator: Exact LS interpolator

			    cout << "Min Max time segmentation" << endl;

			    Mat <CxT1> Ltp;
			    Ltp.ones(1, L);
			    Col <CxT1> ggtp;
			    ggtp.ones(n2, 1);
			    Mat <CxT1> gg;
			    gg = exp(i*fieldMap*tau)*Ltp;
			    Mat <CxT1> iGTGGT;
			    iGTGGT.set_size(L+1, n2);
			    Mat <CxT1> gl;
			    gl.zeros(n2, L);


			    for (unsigned int ii = 0; ii<L; ii++) {
				    for (unsigned int jj = 0; jj<n2; jj++) {
					    gl(jj, ii) = pow(gg(jj, ii), (T1) (ii+1));
				    }
			    }

			    Mat <CxT1> G;
			    G.set_size(n2, L);

			    for (unsigned int jj = 0; jj<L; jj++) {
				    if (jj==0) {
					    G.col(jj) = ggtp;
				    }
				    else {
					    G.col(jj) = gl.col(jj-1);
				    }
			    }

			    Col <CxT1> glsum;
			    Mat <CxT1> GTG;
			    GTG.zeros(L, L);
			    GTG.diag(0) += n2;
			    glsum = sum(gl.t(), 1);
				Mat <CxT1> GTGtp(L, L);
				for (unsigned int ii = 0; ii < (L - 1); ii++) {
					GTGtp.zeros();
				    GTGtp.diag(-(T1) (ii+1)) += glsum(ii);
				    GTGtp.diag((T1) (ii+1)) += std::conj(glsum(ii));
				    GTG = GTG+GTGtp;
			    }

			    T1 rcn = 1/cond(GTG);
			    if (rcn>10*2e-16) { //condition number of GTG
				    iGTGGT = inv(GTG)*G.t();

			    }
			    else {
				    iGTGGT = pinv(GTG)*G.t(); // pseudo inverse
			    }


			    Mat <CxT1> iGTGGTtp;
			    Mat <CxT1> ftp;
			    Col <CxT1> res, temp;

			    for (unsigned int ii = 0; ii<NOneShot; ii++) {
				    ftp = exp(i*fieldMap*timeVec(ii));
				    res = iGTGGT*ftp;
				    tempAA.row(ii) = res.t();
			    }
			    AA = repmat(tempAA, Nshots, 1);
		    }
	    }
		savemat("aamat.mat", "AA", vectorise(AA));
		cout << "Exiting class constructor." << endl;
    }




//Overloaded operators go here

//Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col <CxT1> operator*(const Col <CxT1>& d) const
    {
	    Tobj* G = this->obj;
	    //output is the size of the kspace data
	    Col <CxT1> outData = zeros<Col<CxT1 >> (this->n1);
		//cout << "OutData size = " << this->n1 << endl;
	    Col <CxT1> Wo;
		//uvec dataMaskTrimmed;
	    //loop through time segments
	    for (unsigned int ii = 0; ii<this->L; ii++) {
		    //cout << "Entering time segmentation loop" << endl;
		    //apply a phase to each time segment
		    Wo = exp(-i*(this->fieldMap)*((ii)*this->tau+this->T_min));

		    //perform multiplication by the object and sum up the time segments
			//outData += (this->AA.col(ii))%(*G*(Wo%d));

			//dataMaskTrimmed = find(abs(this->AA.col(ii)) > 0);
			//std::cout << "Length dataMaskTrimmed = " << dataMaskTrimmed.n_rows << std::endl;

			outData += (this->AA.col(ii)) % ((*G).trimmedForwardOp(Wo % d, this->AA.col(ii)));

	    }
	    return outData;
    }

    Col <CxT1> operator/(const Col <CxT1>& d) const
    {
		Tobj *G = this->obj;
	    //output is the size of the image
	    Col <CxT1> outData = zeros<Col<CxT1 >> (this->n2);
	    Col <CxT1> Wo;
	    //loop through the time segemtns
	    for (unsigned int ii = 0; ii<this->L; ii++) {

		    //create the phase map for the Lth time segment
		    Wo = exp(i*(this->fieldMap)*((ii)*this->tau+this->T_min));

		    //perform adjoint operation by the object and sum up the time segments
			outData += Wo % ((*G).trimmedAdjointOp((AA.col(ii) % d), AA.col(ii)));

	    }

	    return outData;

    }



};

#endif
