//
//  main.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/12/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//


#include "PowerGrid.h" //Project headers.
#include "../Support/CeempleMatio.h" //Headers for using savemat and loadmat

using namespace arma; //Armdillo stuff is in the arma namespace
using namespace std; //complex type comes from the STL


int main(int argc, char** argv)
{
    //Create an empty unsized matrix of type real double
    Mat<double> test;

    //Load our read double matrix object. (You need to match type to avoide mangling the data.)
    loadmat("/Users/alexcerjanic/Developer/PowerGrid/Resources/test.mat","test",&test);

    //Create our Gfft object for type complex double and N1=64, N2=64
    Gfft<Col<cx_double>> G(64,64);

    //Perform the forward transformation and store the result in a column vector of type complex double
    //Note the conv_to<type to convert to>::from(data) command to convert our real double phantom to type complex double
    Col<cx_double> testForward = G*vectorise(conv_to<Mat<cx_double>>::from(test));

    //Now perform the adjoint transform
    Col<cx_double> testAdjoint = G/testForward;

    //Testing the my fftshift function. Armadillo delays evaluations to batch them. You can force the evaluation with data.eval()
    //If we just passed the variable rather than the variable wrapped in a function, we wouldn't need the eval.

    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/out.mat","fftShiftOut", fftshift(test).eval());

    //Writing the transforms to disk for evaluation in MATLAB. Note that we can't read or write complex natively (yet), so lets
    //separate real and imaginary first. Also, we can't put multiple variables in a single file yet. These are TODOs.
    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testForwardReal.mat","testForwardReal", real(testForward).eval());
    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testForwardImag.mat","testForwardImag", imag(testForward).eval());

    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testAdjointReal.mat","testAdjointReal", real(testAdjoint).eval());
    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testAdjointImag.mat","testAdjointImag", imag(testAdjoint).eval());


    return 0;
}
