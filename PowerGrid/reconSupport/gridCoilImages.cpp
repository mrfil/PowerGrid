
/*
(C) Copyright 2015-2019 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [gridCoilImages.cpp]

    Synopsis    [GPU accelerated gridding of fully-sampled coil images.]

    Description [This reconstruction supports 2D and 3D reconstructions for fully
                 sampled data with sum of squares (SoS) coil combination.]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [2019/12/8]

 *****************************************************************************/

//Project headers.

#include "../PowerGrid.h"
#include "../processIsmrmrd.hpp"
#include "../processNIFTI.hpp"
#include <boost/program_options.hpp>
#include <chrono>
#include "../reconSupport/directRecon.h"

namespace po = boost::program_options;

int main(int argc, char** argv)
{
    std::string rawDataFilePath, outputImageFilePath, precisionString;

    uword Nx, Ny, Nz, NShots = 1, type = 1;

    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")(
        "inputData,i", po::value<std::string>(&rawDataFilePath)->required(),
        "input ISMRMRD Raw Data file")("outputImage,o", po::value<std::string>(&outputImageFilePath)->required(), "output file path for NIFTIimages")("Nx,x", po::value<uword>(&Nx), "Image size in X")("Ny,y", po::value<uword>(&Ny), "Image size in Y")("Nz,z", po::value<uword>(&Nz), "Image size in Z");

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

    ISMRMRD::Dataset* d;
    ISMRMRD::IsmrmrdHeader hdr;
    acqTracking* acqTrack;
    Col<float> FM;
    Col<std::complex<float> > sen;
    processISMRMRDInput<float>(rawDataFilePath, d, hdr, FM, sen, acqTrack);

    uword numAcq = d->getNumberOfAcquisitions();

    // Grab first acquisition to get parameters (We assume all subsequent
    // acquisitions will be similar).
    ISMRMRD::Acquisition acq;
    d->readAcquisition(0, acq);
    uword nro = acq.number_of_samples();
    uword Ncoils = acq.active_channels();

    // Handle Nx, Ny, Nz
    if (!vm.count("Nx")) {
        Nx = hdr.encoding[0].encodedSpace.matrixSize.x;
    }
    if (!vm.count("Ny")) {
        Ny = hdr.encoding[0].encodedSpace.matrixSize.y;
    }
    if (!vm.count("Nz")) {
        Nz = hdr.encoding[0].encodedSpace.matrixSize.z;
    }

    Col<float> ix, iy, iz;
    initImageSpaceCoords(ix, iy, iz, Nx, Ny, Nz);

    // Check and abort if we have more than one encoding space (No Navigators for
    // now).
    std::cout << "hdr.encoding.size() = " << hdr.encoding.size() << std::endl;
    if (hdr.encoding.size() != 1) {
        std::cout << "There are " << hdr.encoding.size()
                  << " encoding spaces in this file" << std::endl;
        std::cout
            << "This recon does not handle more than one encoding space. Aborting."
            << std::endl;
        return -1;
    }

    int NShotMax = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum + 1;
    int NParMax = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum + 1;
    int NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum + 1;
    int NSetMax = hdr.encoding[0].encodingLimits.set->maximum + 1;
    int NRepMax = hdr.encoding[0].encodingLimits.repetition->maximum + 1;
    int NAvgMax = hdr.encoding[0].encodingLimits.average->maximum + 1;
    int NSegMax = hdr.encoding[0].encodingLimits.segment->maximum + 1;
    int NEchoMax = hdr.encoding[0].encodingLimits.contrast->maximum + 1;
    int NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum + 1;

    std::cout << "NParMax = " << NParMax << std::endl;
    std::cout << "NShotMax = " << NShotMax << std::endl;
    std::cout << "NSliceMax = " << NSliceMax << std::endl;
    std::cout << "NSetMax = " << NSetMax << std::endl;
    std::cout << "NRepMax = " << NRepMax << std::endl;
    std::cout << "NAvgMax = " << NAvgMax << std::endl;
    std::cout << "NEchoMax = " << NEchoMax << std::endl;
    std::cout << "NPhaseMax = " << NPhaseMax << std::endl;
    std::cout << "NSegMax = " << NSegMax << std::endl;

    std::cout << "About to loop through the counters and the file"
              << std::endl;

    std::string baseFilename = "img";
    std::string filename, filenameSOS;
    if (!outputImageFilePath.empty() && *outputImageFilePath.rbegin() != '/') {
        outputImageFilePath += '/';
    }

    //uword NSet = 0; //Set is only used for arrayed ADCs
    //uword NSeg = 0;

    arma::Col<complex<float>> ImageTemp;
    arma::Col<float> imSOS;

    for (uword NPhase = 0; NPhase < NPhaseMax; NPhase++) {
        for (uword NEcho = 0; NEcho < NEchoMax; NEcho++) {
            for (uword NAvg = 0; NAvg < NAvgMax; NAvg++) {
                for (uword NRep = 0; NRep < NRepMax; NRep++) {

                    filename = outputImageFilePath + baseFilename + "_" + "Rep" + std::to_string(NRep) + "_" + "Avg" + std::to_string(NAvg) + "_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);
                    filenameSOS = outputImageFilePath + "imSOS" + "_" + "Rep" + std::to_string(NRep) + "_" + "Avg" + std::to_string(NAvg) + "_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);
                    
                    ImageTemp = gridCoilImages<float>(Nx, Nz, d, &hdr, acqTrack, NPhase, NEcho, NAvg, NRep);
                    // write coil images to file
                    writeNiftiMagPhsImage<float>(filename, ImageTemp, Nx, Ny, Nz, NSliceMax, Ncoils);
                    
                    // calculate the sum of squares
                    imSOS = calcSumOfSquaresImage<float>(ImageTemp, Nx, Nz, NSliceMax, Ncoils);

                    writeNiftiRealImage<float>(filenameSOS, imSOS, Nx, Ny, Nz, NSliceMax);
                }
            }
        }
    }

    //

    // Close ISMRMRD::Dataset, hdr, and acqTrack
    closeISMRMRDData(d, hdr, acqTrack);

return 0;
}
