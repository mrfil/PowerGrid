/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [PowerGridPcSenseMPI_TS.cpp]

    Synopsis    [PowerGrid reconstruction executable supporting ISMRMRD format
                                as input and Phase Corrected Sense Algorithm
                                with MPI support.]

    Description [This reconstruction supports 2D and 3D reconstructions with
                                time segmentation for field correction. This
 support is experimental.]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [06/25/2017]

 *****************************************************************************/

// //Project headers.
#include "PowerGrid.h"
#include "processIsmrmrd.hpp"
#include "processNIFTI.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;
namespace bmpi = boost::mpi;

//using namespace PowerGrid;

int main(int argc, char** argv)
{
    //required for boost::mpi to work.
    bmpi::environment env(argc, argv, true);
    bmpi::communicator world;

    std::string rawDataFilePath, outputImageFilePath, TimeSegmentationInterp;
    uword Nx, Ny, Nz, NShots = 1, NIter = 10, L = 0, type = 1;
    //uword  type, L = 0;
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
        ("NShots,s", po::value<uword>(&NShots), "Number of shots per image")
        ("TimeSegmentationInterp,I", po::value<std::string>(&TimeSegmentationInterp), "Field Correction Interpolator (Required)")
        ("TimeSegments,t", po::value<uword>(&L), "Number of time segments (Required)")
        ("Beta,B", po::value<double>(&beta), "Spatial regularization penalty weight")("Dims2Penalize,D", po::value<uword>(&dims2penalize), "Dimensions to apply regularization to (2 or 3)")
        ("CGIterations,n", po::value<uword>(&NIter), "Number of preconditioned conjugate gradient interations for main solver");

    po::variables_map vm;

    try {

        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 1;
        }

    } catch (boost::program_options::error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }

    arma::Col<float> FM;
    arma::Col<std::complex<float> > sen;
    arma::Col<float> PMap;
    arma::Col<float> FMap;
    arma::Col<std::complex<float> > SMap;

    ISMRMRD::Dataset* d;
    ISMRMRD::IsmrmrdHeader hdr;
    acqTracking* acqTrack;
    uword nro, nc;

    processISMRMRDInput<float>(rawDataFilePath, d, hdr, FM, sen, acqTrack);
    std::cout << "Number of elements in SENSE Map = " << sen.n_rows << std::endl;
    std::cout << "Number of elements in Field Map = " << FM.n_rows << std::endl;

    uword numAcq = d->getNumberOfAcquisitions();

    // Grab first acquisition to get parameters (We assume all subsequent
    // acquisitions will be similar).
    ISMRMRD::Acquisition acq;
    d->readAcquisition(0, acq);
    nro = acq.number_of_samples();
    nc = acq.active_channels();

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
    Col<float> kx(nro), ky(nro), kz(nro), tvec(nro);
    Col<std::complex<float> > data(nro * nc);
    Col<std::complex<float> > ImageTemp(Nx * Ny * Nz);

    if (hdr.encoding.size() != 1) {

        std::cout << "There are " << hdr.encoding.size()
                  << " encoding spaces in this file" << std::endl;
        std::cout
            << "This recon does not handle more than one encoding space. Aborting."
            << std::endl;

        return -1;
    }

    uword NSliceMax, NSetMax, NRepMax, NAvgMax, NSegMax, NEchoMax, NPhaseMax, NParMax, NShotMax;

    NShotMax = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum + 1;
    NParMax = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum + 1;
    NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum + 1;
    NSetMax = hdr.encoding[0].encodingLimits.set->maximum + 1;
    NRepMax = hdr.encoding[0].encodingLimits.repetition->maximum + 1;
    NAvgMax = hdr.encoding[0].encodingLimits.average->maximum + 1;
    NSegMax = hdr.encoding[0].encodingLimits.segment->maximum + 1;
    NEchoMax = hdr.encoding[0].encodingLimits.contrast->maximum + 1;
    NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum + 1;

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

    //uword dataNrows;
    //uword dataNcols;
    //uword PMapNrows;
    //uword kxNrows;

    // Figure out what slices to work on in this rank
    Col<uword> sliceList;
    Col<uword> phaseList;
    Col<uword> echoList;
    Col<uword> avgList;
    Col<uword> repList;

    uword numTasks = NSliceMax * NPhaseMax * NEchoMax * NAvgMax * NRepMax;

    sliceList.resize(numTasks);
    phaseList.resize(numTasks);
    echoList.resize(numTasks);
    avgList.resize(numTasks);
    repList.resize(numTasks);

    // Now we make lists of the tasks such that we are mapping out all the tasks
    for (uword ii = 0; ii < NSliceMax; ii++) { // slice loop
        for (uword jj = 0; jj < NPhaseMax; jj++) { // phase Loop
            for (uword kk = 0; kk < NEchoMax; kk++) { // echo loop
                for (uword ll = 0; ll < NAvgMax; ll++) { // avg Loop
                    for (uword mm = 0; mm < NRepMax; mm++) { // rep loop
                        sliceList(ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = ii;
                        phaseList(ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = jj;
                        echoList(ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = kk;
                        avgList(ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = ll;
                        repList(ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = mm;
                    }
                }
            }
        }
    }

    // Start by creating a 2D vector (vector of vectors) to hold the MPI rank ->
    // task mapping
    vector<uword> init(0);
    std::vector<std::vector<uword> >* taskList = new std::vector<std::vector<uword> >(world.size(), init);
    uword process = 0;
    for (uword task = 0; task < numTasks; task++) {
        (*taskList)[process].push_back(task);

        process++;
        if (process == world.size()) {
            process = 0;
        }
    }
    uword taskIndex;
    std::cout << " Number of tasks = " << numTasks << std::endl;

#ifdef OPENACC_GPU
    // Spread out our tasks across GPUs.
    uword gpunum;
    uword ngpus = acc_get_num_devices(acc_device_nvidia);
    if (ngpus) {
        gpunum = world.rank() % ngpus;
        acc_set_device_num(gpunum, acc_device_nvidia);
    } else {
        acc_set_device_type(acc_device_host);
    }
#endif

    std::cout << "Rank = " << world.rank() << std::endl;

    uword NSlice, NRep, NAvg, NEcho, NPhase;
    uword NSet = 0, NSeg = 0;
    for (uword ii = 0; ii < (*taskList)[world.rank()].size(); ii++) {

        taskIndex = (*taskList)[world.rank()].at(ii);

        NSlice = sliceList(taskIndex);
        NRep = repList(taskIndex);
        NAvg = avgList(taskIndex);
        NEcho = echoList(taskIndex);
        NPhase = phaseList(taskIndex);

        std::cout << "NSlice = " << NSlice << std::endl;
        std::cout << "NRep   = " << NRep << std::endl;
        std::cout << "NAvg   = " << NAvg << std::endl;
        std::cout << "NEcho  = " << NEcho << std::endl;
        std::cout << "NPhase = " << NPhase << std::endl;

        filename = baseFilename + "_" + "Slice" + std::to_string(NSlice) + "_" + "Rep" + std::to_string(NRep) + "_" + "Avg" + std::to_string(NAvg) + "_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);

        getCompleteISMRMRDAcqData<float>(d, acqTrack, NSlice, NRep, NAvg, NEcho, NPhase, data, kx, ky,
            kz, tvec);


       // Deal with the number of time segments
		if(!vm.count("TimeSegments")) {
			switch(type) {
		    	case 1:
			    	L = ceil((arma::max(tvec) - arma::min(tvec))/2E-3);
				break;
				case 2:
					L = ceil((arma::max(tvec) - arma::min(tvec))/3E-3);
				break;
				case 3:
			    	L = ceil((arma::max(tvec) - arma::min(tvec))/3E-3);
				break;
				default: 
					L = 0;
			}
			std::cout << "Info: Setting L = " << L << " by default." << std::endl; 
		}

        PMap = getISMRMRDCompletePhaseMap<float>(d, NSlice, NSet, NRep, NAvg, NPhase, NEcho, NSeg,
            (uword)(Nx * Ny * Nz));
        SMap = getISMRMRDCompleteSENSEMap<std::complex<float> >(d, sen, NSlice, (uword)(Nx * Ny * Nz));
        FMap = getISMRMRDCompleteFieldMap<float>(d, FM, NSlice, (uword)(Nx * Ny * Nz));

        std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
        std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
        std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
        std::cout << "Number of rows in phase map = " << PMap.n_rows << std::endl;
        std::cout << "Number of rows in data = " << data.n_rows << std::endl;
        std::cout << "Number of columns in data = " << data.n_cols << std::endl;
        QuadPenalty<float> R(Nx, Ny, Nz, beta);

        
        pcSenseTimeSeg<float> S_DWI(kx, ky, kz, Nx, Ny, Nz, nc, tvec, L, 1, SMap, FMap,
            0 - PMap);

        ImageTemp = reconSolve<float, pcSenseTimeSeg<float>, QuadPenalty<float> >(data, S_DWI, R, kx, ky, kz, Nx,
            Ny, Nz, tvec, NIter);
        
        writeNiftiMagPhsImage<float>(filename, ImageTemp, Nx, Ny, Nz);
    }

    // Close ISMRMRD::Dataset, hdr, and acqTrack
    closeISMRMRDData(d, hdr, acqTrack);

    return 0;
}
