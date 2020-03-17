/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [oscillateRecon.cpp]

    Synopsis    [PowerGrid reconstruction executable supporting ISMRMRD format
                                as input and Phase Corrected Sense Algorithm for
								use with Oscillate Low Rank MRE recon.]

    Description [This reconstruction supports 2D and 3D reconstructions with
                                time segmentation for field correction. This
 support is experimental.]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [03/16/2020]

 *****************************************************************************/

// //Project headers.
#include "PowerGrid.h"
#include "PGIncludes.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/version.h"
#include "ismrmrd/xml.h"
#include "processIsmrmrd.hpp"
#include "processNIFTI.hpp"
#include <boost/program_options.hpp>
#include "R_lowRank.h"
#include "LRobj.h"

namespace po = boost::program_options;
//using namespace PowerGrid;

int main(int argc, char **argv)
{
	std::string rawDataFilePath, outputImageFilePath, TimeSegmentationInterp;
	uword Nx, Ny, Nz, NShots = 1, type = 1, L = 0, NIter = 10, rank = 1;
	//uword type, L = 0;
	double beta = 0.0;
	uword dims2penalize = 3;
	po::options_description desc("Allowed options");
	desc.add_options()("help,h", "produce help message")(
			"inputData,i", po::value<std::string>(&rawDataFilePath)->required(),
			"input ISMRMRD Raw Data file")
			("outputImage,o", po::value<std::string>(&outputImageFilePath)->required(), "output file path for NIFTIimages")
			("Nx,x", po::value<uword>(&Nx), "Image size in X")
			("Ny,y", po::value<uword>(&Ny), "Image size in Y")
			("Nz,z", po::value<uword>(&Nz), "Image size in Z")
			("Rank,L", po::value<uword>(&rank), "Rank (model order) of the temporal basis")
            ("TimeSegmentationInterp,I", po::value<std::string>(&TimeSegmentationInterp)->required(), "Field Correction Interpolator (Required)")
            ("TimeSegments,t", po::value<uword>(&L)->required(), "Number of time segments (Required)")
			("NShots,s", po::value<uword>(&NShots), "Number of shots per image")
			("Beta,B", po::value<double>(&beta), "Spatial regularization penalty weight")
			("Dims2Penalize,D", po::value<uword>(&dims2penalize), "Dimensions to apply regularization to (2 or 3)")
			("CGIterations,n", po::value<uword>(&NIter),
					"Number of preconditioned conjugate gradient interations for main solver");

	po::variables_map vm;

	try {

		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << desc << std::endl;
			return 1;
		}

	}
	catch (boost::program_options::error& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cout << desc << std::endl;
		return 1;
	}

	arma::Col<float> FM;
	arma::Col<std::complex<float>> sen;
	
	arma::Col<float> FMap;
	arma::Col<std::complex<float>> SMap;
	ISMRMRD::Dataset* d;
	ISMRMRD::IsmrmrdHeader hdr;
	acqTracking* acqTrack;
	processISMRMRDInput<float>(rawDataFilePath, d, hdr, FM, sen, acqTrack);

	std::cout << "Number of elements in SENSE Map = " << sen.n_rows << std::endl;
	std::cout << "Number of elements in Field Map = " << FM.n_rows << std::endl;

	uword numAcq = d->getNumberOfAcquisitions();

	// Grab first acquisition to get parameters (We assume all subsequent
	// acquisitions will be similar).
	ISMRMRD::Acquisition acq;
	d->readAcquisition(0, acq);
	uword nro = acq.number_of_samples();
	uword nc = acq.active_channels();

	// Handle Nx, Ny, Nz
	if(!vm.count("Nx")) {
		Nx = hdr.encoding[0].encodedSpace.matrixSize.x;
	} 
	if(!vm.count("Ny")) {
		Ny = hdr.encoding[0].encodedSpace.matrixSize.y;
	} 
	if(!vm.count("Nz")) {
		Nz = hdr.encoding[0].encodedSpace.matrixSize.z;
	} 

	Col<float> ix, iy, iz;
	initImageSpaceCoords(ix, iy, iz, Nx, Ny, Nz);
	Col<float> kxTemp(nro), kyTemp(nro), kzTemp(nro), tvecTemp(nro);


	// Check and abort if we have more than one encoding space (No Navigators for
	// now).
	if (hdr.encoding.size()!=1) {
		std::cout << "There are " << hdr.encoding.size()
		          << " encoding spaces in this file" << std::endl;
		std::cout
				<< "This recon does not handle more than one encoding space. Aborting."
				<< std::endl;
		return -1;
	}

	uword NShotMax = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum + 1;
	uword NParMax = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum + 1;
	uword NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum + 1;
	uword NSetMax = hdr.encoding[0].encodingLimits.set->maximum + 1;
	uword NRepMax = hdr.encoding[0].encodingLimits.repetition->maximum + 1;
	uword NAvgMax = hdr.encoding[0].encodingLimits.average->maximum + 1;
	uword NSegMax = hdr.encoding[0].encodingLimits.segment->maximum + 1;
	uword NEchoMax = hdr.encoding[0].encodingLimits.contrast->maximum + 1;
	uword NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum + 1;

	std::cout << "NParMax = " << NParMax << std::endl;
	std::cout << "NShotMax = " << NShotMax << std::endl;
	std::cout << "NSliceMax = " << NSliceMax << std::endl;
	std::cout << "NSetMax = " << NSetMax << std::endl;
	std::cout << "NRepMax = " << NRepMax << std::endl;
	std::cout << "NAvgMax = " << NAvgMax << std::endl;
	std::cout << "NEchoMax = " << NEchoMax << std::endl;
	std::cout << "NPhaseMax = " << NPhaseMax << std::endl;
	std::cout << "NSegMax = " << NSegMax << std::endl;
	std::cout << "About to loop through the counters and scan the file"
	          << std::endl;
  	std::string baseFilename = "img";
	std::string filename;
	if (!outputImageFilePath.empty() && *outputImageFilePath.rbegin() != '/') {
    	outputImageFilePath += '/';
	}

	Mat<std::complex<float>> data(nro * nc * NShotMax * NParMax, NRepMax);
	Col<std::complex<float>> dataTemp;
	Mat<float> kx(nro * NShotMax * NParMax, NRepMax), 
			   ky(nro * NShotMax * NParMax, NRepMax), 
			   kz(nro * NShotMax * NParMax, NRepMax), 
			   tvec(nro * NShotMax * NParMax, NRepMax);

	Mat<float> PMap(Nx * Ny * Nz * NShotMax * NParMax, NRepMax);
	Col<complex<float>> vTemp = getISMRMRDTemporalBasis<float>(d);
	Mat<complex<float>> v = reshape(vTemp, NRepMax, NRepMax);
	
	Mat<complex<float>> vLR = trans(v.cols(0,rank-1)).eval();
	// Declare array of pcSENSE pointers to generate obejcts for the low rank object.
	pcSenseTimeSeg<float>** A_list = NULL; 

	Col<std::complex<float>> ImageTemp(Nx * Ny * Nz * NRepMax);

	uword NSet = 0; //Set is only used for arrayed ADCs
	uword NSeg = 0;
	for (uword NPhase = 0; NPhase < NPhaseMax; NPhase++) {
		for (uword NEcho = 0; NEcho < NEchoMax; NEcho++) {
			for (uword NAvg = 0; NAvg < NAvgMax; NAvg++) {
				for (uword NSlice = 0; NSlice < NSliceMax; NSlice++) {

					filename = outputImageFilePath + baseFilename + "_" + "Slice" + std::to_string(NSlice) +
						"_" + "Avg" + std::to_string(NAvg) +
						"_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);

					SMap = getISMRMRDCompleteSENSEMap<std::complex<float>>(d, sen, NSlice, (uword) (Nx*Ny*Nz));
					FMap = getISMRMRDCompleteFieldMap<float>(d, FM, NSlice, (uword) (Nx*Ny*Nz));

					//Allocate array of pcSenseTimeSeg objects
					A_list = new pcSenseTimeSeg<float> *[NRepMax];

					for (uword NRep = 0; NRep < NRepMax; NRep++) {

						getCompleteISMRMRDAcqData<float>(d, acqTrack, NSlice, NRep, NAvg, NEcho, NPhase, dataTemp, kxTemp, kyTemp,
								kzTemp, tvecTemp);
						cout << "about to copy in rep data" << endl;
						// Collect and assemble a complete data vector (once the matrix is flattened with vectorise() )
						cout << "data size = " << data.n_rows << " and dataTemp size = " << dataTemp.n_rows << endl;
						data.col(NRep) = dataTemp;
						cout << "kx size = " << kx.n_rows << " and kxTemp size = " << kxTemp.n_rows << endl;
						kx.col(NRep)   = kxTemp;
						cout << "ky size = " << ky.n_rows << " and kyTemp size = " << kyTemp.n_rows << endl;
						ky.col(NRep)   = kyTemp;
						cout << "kz size = " << kz.n_rows << " and kzTemp size = " << kzTemp.n_rows << endl;
						kz.col(NRep)   = kzTemp;
						cout << "tvec size = " << tvec.n_rows << " and tvecTemp size = " << tvecTemp.n_rows << endl;
						tvec.col(NRep) = tvecTemp;
						
						// Deal with the number of time segments
						if(!vm.count("TimeSegments")) {
							switch(type) {
								case 1:
									L = ceil((arma::max(tvecTemp) - arma::min(tvecTemp))/2E-3);
								break;
								case 2:
									L = ceil((arma::max(tvecTemp) - arma::min(tvecTemp))/3E-3);
								break;
								case 3:
									L = ceil((arma::max(tvecTemp) - arma::min(tvecTemp))/3E-3);
								break;
								default: 
									L = 0;
							}
							std::cout << "Info: Setting L = " << L << " by default." << std::endl; 
						}

						PMap.col(NRep) = getISMRMRDCompletePhaseMap<float>(d, NSlice, NSet, NRep, NAvg, NPhase, NEcho, NSeg,
																(uword) (Nx*Ny*Nz));

						cout << "Returned from getISMRMRDCompletePhaseMap" << endl;

						A_list[NRep] = new pcSenseTimeSeg<float>(kx.col(NRep), ky.col(NRep), kz.col(NRep), Nx, Ny, Nz, nc, tvec.col(NRep), L, type, SMap, FMap,
															   0-PMap.col(NRep));

					} // End of rep loop
					
					cout << "End of rep loop" << endl;

					std::cout << "Number of elements in SMap = " << SMap.n_rows << std::endl;
					std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
					std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
					std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
					std::cout << "Number of rows in phase map = " << PMap.n_rows << std::endl;
					std::cout << "Number of rows in data = " << data.n_rows << std::endl;

					std::cout << "Number of columns in data = " << data.n_cols << std::endl;

					LRobj<float, pcSenseTimeSeg<float>> A_lr(Nx * Ny * Nz, nro * nc * NShotMax * NParMax, rank, NRepMax, vLR, A_list);
					QuadPenalty<float> R(Nx, Ny, Nz, beta, dims2penalize);
					R_lowRank<float, QuadPenalty<float>> R_lr(R, rank, NRepMax, v);

					Col<float> W;
  					W.ones(data.n_elem);
  					Col<complex<float>> xinit;
  					xinit.zeros(Nx * Ny * Nz * rank);

  					Col<complex<float>> imageOut;
  					imageOut = solve_pwls_pcg<float, LRobj<float, pcSenseTimeSeg<float>>, R_lowRank<float, QuadPenalty<float>>>(xinit, A_lr, W, vectorise(data), R_lr, NIter);
					
					Mat<complex<float>> imgB = reshape(imageOut, Nx * Ny * Nz, rank);
					
					// Form complete image set
					ImageTemp = vectorise(imgB * vLR);
					
					writeNiftiMagPhsImage<float>(filename, ImageTemp, Nx, Ny, Nz, NRepMax);
					 
				}
			}
		}
	}

  // Close ISMRMRD::Dataset
  delete d;

  return 0;
}
