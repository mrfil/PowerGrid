/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [PowerGrid.h]

    Synopsis    [Main header file for the project.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef PowerGrid_PowerGrid_h
#define PowerGrid_PowerGrid_h
//#define ARMA_NO_DEBUG // Disable this comment only for release.
#include <cmath>
#include <iostream>
#include <string>

#ifdef PowerGridMPI
#include "../Support/ArmaExtensions/arma_extensions.h"
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#endif // PowerGridMPI

// Support Headers for making it easier to work with Armadillo and Matio.
//#include "../Support/CeempleComplex.h"

#include "../Support/CeempleArmadillo.h"
#include "../Support/CeempleMatio.h"

#include "config.hxx"
#include "griddingTypes.h"
#include <memory>

// Headers for ISMRMRD Support
//#include "ismrmrd/ismrmrd.h"
//#include "ismrmrd/xml.h"
//#include "ismrmrd/dataset.h"
//#include "ismrmrd/version.h"

namespace PowerGrid {

#include "Gdft.h"
#include "Gfft.h"
#include "Gnufft.h"
#include "QuadPenalty.h"
#include "Robject.h"
#include "SENSE.h"
#include "TVPenalty.h"
#include "TimeSegmentation.h"
#include "fftGPU.h"
#include "fftshift.hpp"
#include "ftCpu.h"
#include "gridding.h"
#include "griddingSupport.hpp"
#include "pcSENSE.h"
#include "reconSolve.hpp"
#include "solve_pwls_pcg.hpp"

#ifdef PowerGridMPI
#include "mpipcSENSE.h"
#endif // PowerGridMPI

#include "DWIRecon.h"
#include "GdftRecon.h"
#include "GnufftRecon.h"
}

#endif
