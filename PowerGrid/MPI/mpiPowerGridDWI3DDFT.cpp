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

    uword Nx, Ny, Nz, Niter = 1, NL = 1, Ncoils;
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
            auto_ptr<PowerGridConfig_t> cfg(PowerGridConfig(configPath.c_str()));

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

        int test = reconfMRIGgrid<double>(testPath, Nx, Ny, Nz, NL, Niter, Ncoils, startIndex, endIndex);
    } else {

        try {
            auto_ptr<PowerGridConfig_t> cfg(PowerGridConfig(configPath.c_str()));

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

        //int test = test_SpeedCompareGgrid<double>(testPath, Nx,Ny,Nz,NL,Niter,Ncoils);
        int test = test_DWIDft<float>(testPath, Nx, Ny, Nz, NL, Niter, Ncoils, env, world);

    }

    return 0;
}

