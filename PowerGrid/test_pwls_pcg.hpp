/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [test_pwls_pcg.hpp]

    Synopsis    []

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/


#ifndef PowerGrid_test_TimeSegmentation_h
#define PowerGrid_test_TimeSegmentation_h


template<typename T1>
Mat<std::complex<T1>> test_pwls_pcg()
{
    typedef std::complex<T1> CxT1;
    std::string testPath = "/shared/mrfil-data/jholtrop/repos/PowerGrid/Resources/";

    // Load real data
    Mat<T1> test;

    //Load our read double matrix object. (You need to match type to avoide mangling the data.)
    loadmat(testPath+"test.mat","test",&test);

    // Formard operator
    Gfft<T1> G(64,64);
    Col<CxT1> TestForward;
    TestForward = G *vectorise(conv_to<Mat<cx_double>>::from(test));

    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(64,64);

    QuadPenalty<T1> R(64, 64, 1, 1, ReconMask);

    uword niter = 2;
    Col<T1> xinit = zeros<Col<CxT1>>(64*64); // initial estimate of xd
    double W;
    W = 1.0;

    Col<CxT1> x_t;
    x_t = solve_pwls_pcg<T1,Gfft<T1>,QuadPenalty<T1>>(xinit, G, W, TestForward, R, niter);

    return x_t;
}

#endif
