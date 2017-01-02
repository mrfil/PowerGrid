/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [GdftRecon.h]

    Synopsis    [Implements a basic 2D and 3D recon with the DFT and field
                    correction.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/
#ifndef __PowerGrid__GdftRecon
#define __PowerGrid__GdftRecon

#include "PGIncludes.h"
#include "Gdft.h"
#include "SENSE.h"
#include "QuadPenalty.h"
#include "solve_pwls_pcg.hpp"

using namespace arma;
using namespace std;

template <typename T1>
int GdftRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
              uword niter, uword nc, uword nshots, T1 beta);

// Explicit Instantiations
extern template int GdftRecon<float>(string, uword, uword, uword, uword, uword,
                                     uword, uword, float);
extern template int GdftRecon<double>(string, uword, uword, uword, uword, uword,
                                      uword, uword, double);

#endif // GdftRecon
