/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [DWIRecon.h]

    Synopsis    [Code implementing diffusion weighted image reconstructions
                    using MPI using the DFT object.]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC
                 0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/13/2016]

*****************************************************************************/

#ifndef __PowerGrid__DWIRecon
#define __PowerGrid__DWIRecon

#include "PGIncludes.h"
#include "Gdft.h"
#include "pcSENSE.h"
#include "mpipcSENSE.h"
#include "QuadPenalty.h"
#include "solve_pwls_pcg.hpp"

using namespace arma;
using namespace std;

template <typename T1>
int DWIRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
             uword niter, uword nc);

// Explicit Instantiations
extern template int DWIRecon<float>(string, uword, uword, uword, uword, uword,
                                    uword);
extern template int DWIRecon<double>(string, uword, uword, uword, uword, uword,
                                     uword);

#ifdef PowerGridMPI
// boost::mpi version of DWIDft reconstruction routine.
template <typename T1>
int DWIRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
             uword niter, uword nc, bmpi::environment &env,
             bmpi::communicator &world);

// Explicit Instantiations
extern template int DWIRecon<float>(string dataPath, uword Nx, uword Ny,
                                    uword Nz, uword L, uword niter, uword nc,
                                    bmpi::environment &env,
                                    bmpi::communicator &world);
extern template int DWIRecon<double>(string dataPath, uword Nx, uword Ny,
                                     uword Nz, uword L, uword niter, uword nc,
                                     bmpi::environment &env,
                                     bmpi::communicator &world);

template <typename T1>
int DWIRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
             uword niter, uword nc, uword nimage, uword nphase, uword nslab,
             bmpi::environment &env, bmpi::communicator &world);

// Explicit Instantiations
extern template int DWIRecon<float>(string dataPath, uword Nx, uword Ny,
                                    uword Nz, uword L, uword niter, uword nc,
                                    uword nimage, uword nphase, uword nslab,
                                    bmpi::environment &env,
                                    bmpi::communicator &world);
extern template int DWIRecon<double>(string dataPath, uword Nx, uword Ny,
                                     uword Nz, uword L, uword niter, uword nc,
                                     uword nimage, uword nphase, uword nslab,
                                     bmpi::environment &env,
                                     bmpi::communicator &world);

#endif // PowerGridMPI
template <typename T1>
int DWIRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
             uword niter, uword nc, uword nimage, uword nphase, uword nslab);

// Explicit Instantiations
extern template int DWIRecon<float>(string dataPath, uword Nx, uword Ny,
                                    uword Nz, uword L, uword niter, uword nc,
                                    uword nimage, uword nphase, uword nslab);
extern template int DWIRecon<double>(string dataPath, uword Nx, uword Ny,
                                     uword Nz, uword L, uword niter, uword nc,
                                     uword nimage, uword nphase, uword nslab);

#endif // Include Guard __PowerGrid__DWIRecon
