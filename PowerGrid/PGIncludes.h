
/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [PGIncludes.h]

    Synopsis    [Main header file for internal use.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/27/2016]

*****************************************************************************/
#ifndef PGINCLUDES_H
#define PGINCLUDES_H

#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>
#include "../Support/ArmaExtensions/permute.hpp"

#ifdef PowerGridMPI
#include "../Support/ArmaExtensions/arma_extensions.h"
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#endif // PowerGridMPI

// Support Headers for making it easier to work with Armadillo and Matio.
//#include "../Support/CeempleComplex.h"
#ifdef _OPENACC
#include "nvToolsExt.h"
#endif

#include "../Support/CeempleArmadillo.h"
#include "../Support/CeempleMatio.h"

#include "config.hxx"
#include "griddingTypes.h"
#include <memory>

// Numeric constants according to the precision type.
#ifdef ENABLE_DOUBLE_PRECISION
#define MRI_PI 3.1415926535897932384626433832795029
#define MRI_NN 64
#define MRI_DELTAZ 0.003
#define MRI_ZERO 0.0
#define MRI_ONE 1.0
#define MRI_NEG_ONE -1.0
#define MRI_POINT_FIVE 0.5
#define MRI_SMOOTH_FACTOR 0.0000001
#else
#define MRI_PI 3.1415926535897932384626433832795029f
#define MRI_NN 64
#define MRI_DELTAZ 0.003f
#define MRI_ZERO 0.0f
#define MRI_ONE 1.0f
#define MRI_NEG_ONE -1.0f
#define MRI_POINT_FIVE 0.5f
#define MRI_SMOOTH_FACTOR 0.000001f
#endif

#endif // PGINCLUDES_H
