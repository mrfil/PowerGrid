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
#include <iostream>
#include <cmath>
#include <string>

#ifdef PowerGridMPI
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include "../Support/ArmaExtensions/arma_extensions.h"
#endif //PowerGridMPI

//Support Headers for making it easier to work with Armadillo and Matio.
//#include "../Support/CeempleComplex.h"


#include "../Support/CeempleArmadillo.h"
#include "../Support/CeempleMatio.h"

#include "griddingTypes.h"
#include <memory>
#include "config.hxx"


//Headers for ISMRMRD Support
//#include "ismrmrd/ismrmrd.h"
//#include "ismrmrd/xml.h"
//#include "ismrmrd/dataset.h"
//#include "ismrmrd/version.h"


namespace PowerGrid {

#include "fftshift.hpp"
#include "ftCpu.hpp"
#include "Gdft.hpp"
#include "fftGPU.hpp"
#include "griddingSupport.hpp"
#include "gridding.hpp"
#include "Gnufft.hpp"
#include "Gfft.hpp"
#include "SENSE.hpp"
#include "Robject.hpp"
#include "QuadPenalty.hpp"
#include "TVPenalty.hpp"
#include "solve_pwls_pcg.hpp"
#include "reconSolve.hpp"
#include "TimeSegmentation.hpp"
#include "pcSENSE.hpp"

#ifdef PowerGridMPI
#include "mpipcSENSE.hpp"
#endif //PowerGridMPI

#include "DWIRecon.hpp"
#include "GdftRecon.hpp"
#include "GnufftRecon.hpp"
}

#endif
