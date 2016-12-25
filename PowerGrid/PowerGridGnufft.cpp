/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [PowerGridGnufft.cpp]

    Synopsis    [Field corrected recon using Gdft object.]

    Description []

    Revision    [0.1.0; Joseph Holtrop, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#include "../Support/CeempleMatio.h" //Headers for using savemat and loadmat
#include "PowerGrid.h"               //Project headers.

using namespace arma; // Armdillo stuff is in the arma namespace
using namespace std;  // complex type comes from the STL
using namespace PowerGrid;

int main(int argc, char **argv) {

  uword Nx, Ny, Nz, Niter = 1, NL = 1, Ncoils, Nshots = 1;
  uword startIndex, endIndex;
  double Beta = 0.0;
  string testPath, configPath;
  if (argc > 1) {
    testPath = std::string(argv[1]);
    configPath = testPath + "config.xml";
  } else {
    cout << "Enter a path to find test files." << endl;
    return -1;
  }

  try {

    auto_ptr<PowerGridConfig_t> cfg(PowerGridConfig(configPath.c_str()));

    Nx = cfg->Nx();
    Ny = cfg->Ny();
    Nz = cfg->Nz();
    NL = cfg->Ntimeseg();
    Niter = cfg->Niter();
    Ncoils = cfg->Ncoils();
    Nshots = cfg->Nshots();
    Beta = cfg->Beta();
  } catch (const xml_schema::exception &e) {
    cerr << e << endl;
    return 1;
  }
  int test =
      GnufftRecon<float>(testPath, Nx, Ny, Nz, NL, Niter, Ncoils, Nshots, Beta);

  return 0;
}
