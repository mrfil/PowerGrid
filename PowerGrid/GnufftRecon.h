/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [GnufftRecon.h]

    Synopsis    [Implements a basic 2D and 3D recon with NUFFT and Time
                                Segmentation for field correction.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/26/2016]

 *****************************************************************************/

#ifndef __PowerGrid__GnufftRecon
#define __PowerGrid__GnufftRecon

#include "PowerGrid.h"

using namespace arma;

template <typename T1>
int GnufftRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
                uword niter, uword nc, uword nshots, T1 beta);

// Explicit Instantiation
extern template int GnufftRecon<float>(string, uword, uword, uword, uword,
                                       uword, uword, uword, float);
extern template int GnufftRecon<double>(string, uword, uword, uword, uword,
                                        uword, uword, uword, double);

#endif //__PowerGrid__GnufftRecon
