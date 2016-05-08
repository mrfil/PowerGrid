/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [griddingTypes.hpp]

    Synopsis    [Support times.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef PowerGrid_griddingTypes_h
#define PowerGrid_griddingTypes_h
template<typename T1>
struct parameters{
    int numSamples;
    int imageSize[3];
    int gridSize[3];
    T1 gridOS;
    T1 kernelWidth;
    int binsize;
    int useLUT;
    int sync;
};

template <typename T1>
struct ReconstructionSample{
    T1 real;
    T1 imag;
    T1 kX;
    T1 kY;
    T1 kZ;
    T1 sdc;
    T1 t;
    T1 dummy;
};

#endif
