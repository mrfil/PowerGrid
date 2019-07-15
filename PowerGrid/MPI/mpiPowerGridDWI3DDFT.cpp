//
//  PowerGridDWI3DDFT.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic 11/4/2015.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#include "../PowerGrid.h"
namespace mpi = boost::mpi;

using namespace arma; //Armdillo stuff is in the arma namespace
using namespace std; //complex type comes from the STL
using namespace PowerGrid;

int main(int argc, char **argv) {

    mpi::environment env(argc, argv, true);
    mpi::communicator world;

    arma::Col<float> FM;
    arma::Col<std::complex<float>> sen;
    arma::Col<float> PMap;
    ISMRMRD::Dataset *d;
    ISMRMRD::IsmrmrdHeader hdr;
    processISMRMRDInput<float>(rawDataFilePath, d, hdr, FM, sen);

    //savemat("testFM.mat", "FM", FM);
    //savemat("testSen.mat", "sen", sen);

    std::cout << "Number of elements in SENSE Map = " << sen.n_rows << std::endl;
    std::cout << "Number of elements in Field Map = " << FM.n_rows << std::endl;

    uword numAcq = d->getNumberOfAcquisitions();

    // Grab first acquisition to get parameters (We assume all subsequent
    // acquisitions will be similar).
    ISMRMRD::Acquisition acq;
    d->readAcquisition(0, acq);
    uword nro = acq.number_of_samples();
    uword nc = acq.active_channels();

    Col<float> ix, iy, iz;
    initImageSpaceCoords(ix, iy, iz, Nx, Ny, Nz);
    Col<float> kx(nro), ky(nro), kz(nro), tvec(nro);
    Col<std::complex<float>> data(nro * nc);
    Col<std::complex<float>> ImageTemp(Nx * Ny * Nz);

    // Check and abort if we have more than one encoding space (No Navigators for
    // now).
    if (hdr.encoding.size() != 1) {
        if (world.rank() == 0) {
            std::cout << "There are " << hdr.encoding.size()
                      << " encoding spaces in this file" << std::endl;
            std::cout
                    << "This recon does not handle more than one encoding space. Aborting."
                    << std::endl;
        }
        return -1;
    }


    uword NSliceMax = hdr.encoding[0].encodingLimits.slice;
    uword NSetMax   = hdr.encoding[0].encodingLimits.set;
    uword NRepMax   = hdr.encoding[0].encodingLimits.repetition;
    uword NAvgMax   = hdr.encoding[0].encodingLimits.average;
    uword NSegMax   = hdr.encoding[0].encodingLimits.segment;
    uword NEchoMax  = hdr.encoding[0].encodingLimits.contrast;
    uword NPhaseMax = hdr.encoding[0].encodingLimits.phase;
    if (world.rank() == 0) {
        std::cout << "NSliceMax = " << NSliceMax << std::endl;
        std::cout << "NSetMax = " << NSetMax << std::endl;
        std::cout << "NRepMax = " << NRepMax << std::endl;
        std::cout << "NAvgMax = " << NAvgMax << std::endl;
        std::cout << "NEchoMax = " << NEchoMax << std::endl;
        std::cout << "NPhaseMax = " << NPhaseMax << std::endl;
        std::cout << "NSegMax = " << NSegMax << std::endl;
        std::cout << "About to loop through the counters and scan the file"
                  << std::endl;
    }
    uword NSet = 0; //Set is only used for arrayed ADCs
    uword NSeg = 0;
    for (uword NEcho = 0; NEcho <= NEchoMax; NEcho++) {
        for (uword NSeg = 0; NSeg <= NSegMax; NSeg++) {
            for (uword NRep = 0; NRep < NRepMax; NRep++) {
                for( uword NAvg = 0; NAvg < 1; NAvg++) {
                    for (uword NPhase = 0; NPhase <= NPhaseMax; NPhase++) {
                        //for (uword NSet = 0; NSet < NSetMax + 1; NSet++) {
                        for (uword NSlice = 0; NSlice <= NSliceMax; NSlice++) {
                            getCompleteISMRMRDAcqData<float>(d, NSlice, NSet, NRep, NAvg, data, kx, ky,
                                    kz, tvec);

                            PMap = getISMRMRDCompletePhaseMap<float>(d, NSlice, NSet, NRep, NAvg, NPhase, NEcho, NSeg, (uword)(Nx*Ny*Nz));

                            std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
                            std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
                            std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
                            std::cout << "Number of rows in phase map = " << PMap.n_rows << std::endl;
                            std::cout << "Number of rows in data = " << data.n_rows << std::endl;

                            std::cout << "Number of columns in data = " << data.n_cols << std::endl;

                            //Gnufft<float> A(kx.n_rows, (float) 2.0, Nx, Ny, Nz, kx, ky, kz, ix,
                            // iy, iz);
                            //Gdft<float> A(kx.n_rows, Nx*Ny*Nz,kx,ky,kz,ix,iy,iz,FM,tvec);
                            pcSENSE<float> S_DWI(kx, ky, kz, Nx, Ny, Nz, nc, tvec, sen, FM,
                                    0 - PMap);
                            //pcSENSE<float, Gnufft<float>> Sg(A, sen, kx.n_rows, Nx*Ny*Nz, nc);
                            QuadPenalty<float> R(Nx, Ny, Nz, beta);

                            ImageTemp = reconSolve<float, pcSENSE<float>, QuadPenalty<float>>(data, S_DWI, R, kx, ky, kz, Nx,
                                    Ny, Nz, tvec, NIter);
                            writeISMRMRDImageData<float>(d, ImageTemp, Nx, Ny, Nz);
                        }
                    }
                }
            }
        }
    }
    // Close ISMRMRD::Dataset
    delete d;

    return 0;
}
