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

#include "PGIncludes.h"

// Headers for ISMRMRD Support
//#include "ismrmrd/ismrmrd.h"
//#include "ismrmrd/xml.h"
//#include "ismrmrd/dataset.h"
//#include "ismrmrd/version.h"

//namespace PowerGrid {
#include "Robject.h"
#include "TVPenalty.h"

#include "Gdft.h"
#include "Gfft.h"
#include "Gnufft.h"

#include "pcSENSE.h"
#include "solve_pwls_pcg.hpp"

#include "SENSE.h"
#include "TimeSegmentation.h"
#include "fftGPU.h"
#include "fftshift.hpp"
#include "ftCpu.h"
#include "gridding.h"
#include "griddingSupport.h"

#include "reconSolve.h"

#ifdef PowerGridMPI
#include "mpipcSENSE.h"
#endif // PowerGridMPI

#include "DWIRecon.h"
#include "GdftRecon.h"
#include "GnufftRecon.h"

//}

#endif
