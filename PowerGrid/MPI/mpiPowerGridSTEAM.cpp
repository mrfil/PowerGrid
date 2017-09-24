//
//  mpiPowerGridDWI3DDFT.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic 11/4/2015.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#include "../PowerGrid.h"
#include "../../Support/CeempleMatio.h"


namespace mpi = boost::mpi;

using namespace arma; //Armdillo stuff is in the arma namespace
using namespace std; //complex type comes from the STL
//using namespace PowerGrid;

int main(int argc, char **argv) {

        mpi::environment env(argc, argv, true);
        mpi::communicator world;

        uword Nx, Ny, Nz, Niter = 1, NL = 1, Ncoils, Nimages, Nphases, Nslabs;
        uword startIndex, endIndex;

        string testPath, configPath;
        if (argc > 1) {
                testPath = std::string(argv[1]);
                configPath = testPath + "config.xml";
        } else {
                cout << "Enter a path to find test files." << endl;
                return -1;
        }

        if (argc > 2) {
                startIndex = atoi(argv[2]);
                endIndex = atoi(argv[3]);
                cout << "Start Index." << startIndex << "End Index" << endIndex << endl;

                try {
                        auto_ptr <PowerGridConfig_t> cfg(PowerGridConfig(configPath.c_str()));

                        Nx = cfg->Nx();
                        Ny = cfg->Ny();
                        Nz = cfg->Nz();
                        NL = cfg->Ntimeseg();
                        Niter = cfg->Niter();
                        Ncoils = cfg->Ncoils();

                }
                catch (const xml_schema::exception &e) {
                        cerr << e << endl;
                        return 1;
                }

                //int test = reconfMRIGnufft<double>(testPath, Nx, Ny, Nz, NL, Niter, Ncoils, startIndex, endIndex);
        } else {
                std::cout << configPath.c_str() << std::endl;
                try {
                        auto_ptr <PowerGridConfig_t> cfg(PowerGridConfig(configPath.c_str()));

                        Nx = cfg->Nx();
                        Ny = cfg->Ny();
                        Nz = cfg->Nz();
                        NL = cfg->Ntimeseg();
                        Niter = cfg->Niter();
                        Ncoils = cfg->Ncoils();
                        Nimages = cfg->Nimages();
                        Nphases = cfg->Nphases();
                        Nslabs = cfg->Nslabs();

                }
                catch (const xml_schema::exception &e) {
                        cerr << e << endl;
                        return 1;
                }

                //int test = test_SpeedCompareGnufft<double>(testPath, Nx,Ny,Nz,NL,Niter,Ncoils);
                for (uword slab = 0; slab < Nslabs; slab++) {
                        for (uword image = 0; image < Nimages; image++) {
                                for (uword phase = 0; phase < Nphases; phase++) {
                                        int test = DWIRecon<float>(testPath, Nx, Ny, Nz, NL, Niter, Ncoils, image, phase, slab, env,
                                                                      world);
                                }
                        }

                }

        }

        return 0;
}
