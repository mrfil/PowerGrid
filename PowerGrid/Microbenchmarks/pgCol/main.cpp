
#include "../../PGIncludes.h"

#include "../../pgCol.hpp"
#include "../../pgMat.hpp"
#include "../../pgComplex.hpp"

#ifdef _OPENACC
#include "openacc.h"
#endif

#include <iostream>
#include <complex>
#include <chrono>

using namespace std;


int main()
{
    typedef std::complex<float> CxT1;
    arma::uword N1 = 128*128;

    arma::uword N2 = 256*256;
    arma::uword N3 = 20;
    // Establish Armadillo timeframes
    auto start = std::chrono::high_resolution_clock::now();
    {
    Col<CxT1> A1 = ones<Col<CxT1>>(N1);
    Col<CxT1> A2 = ones<Col<CxT1>>(N1);
    
    Col<CxT1> A = A1 % A2 + A1;
    
    Col<CxT1> B1 = ones<Col<CxT1>>(N2);
    Col<CxT1> B2 = ones<Col<CxT1>>(N2);

    Col<CxT1> B = B1 % B2 + B1;
    }
    auto finish = std::chrono::high_resolution_clock::now();

    #ifdef _OPENACC
    acc_init(acc_device_default);
    #endif
    std::chrono::duration<float> elapsed = finish - start;
    std::cout << "arma::Col<complex<float>> Test Elapsed time: " << elapsed.count() << " s\n";
  
    start = std::chrono::high_resolution_clock::now();
    {
    pgCol<pgComplex<float>> A1(N1);
    pgCol<pgComplex<float>> A2(N1);

    A1.ones();
    A2.ones();

    pgCol<pgComplex<float>> A = A1 % A2 + A1;
    //A1+=A2;

    pgCol<pgComplex<float>> B1(N2);
    pgCol<pgComplex<float>> B2(N2);

    B1.ones();
    B2.ones();

    pgCol<pgComplex<float>> B = B1 % B2 + B1;
    }
    finish = std::chrono::high_resolution_clock::now();

    elapsed = finish - start;
    std::cout << "pgCol<pgComplex<float>> For Elapsed time: " << elapsed.count() << " s\n";
    
    start = std::chrono::high_resolution_clock::now();
    {
    Mat<CxT1> A1 = ones<Mat<CxT1>>(N1,N3);
    Mat<CxT1> A2 = ones<Mat<CxT1>>(N1,N3);
    
    Mat<CxT1> A = A1 % A2 + A1;
    
    Mat<CxT1> B1 = ones<Mat<CxT1>>(N2,N3);
    Mat<CxT1> B2 = ones<Mat<CxT1>>(N2,N3);

    Mat<CxT1> B = B1 % B2 + B1;
    }
    finish = std::chrono::high_resolution_clock::now();

    #ifdef _OPENACC
    acc_init(acc_device_default);
    #endif
    elapsed = finish - start;
    std::cout << "arma::Mat<complex<float>> Test Elapsed time: " << elapsed.count() << " s\n";
  
    start = std::chrono::high_resolution_clock::now();
    {
    pgMat<pgComplex<float>> A1(N1,N3);
    pgMat<pgComplex<float>> A2(N1,N3);

    A1.ones();
    A2.ones();

    pgMat<pgComplex<float>> A = A1 % A2 + A1;
    //A1+=A2;

    pgMat<pgComplex<float>> B1(N2,N3);
    pgMat<pgComplex<float>> B2(N2,N3);

    B1.ones();
    B2.ones();

    pgMat<pgComplex<float>> B = B1 % B2 + B1;
    }
    finish = std::chrono::high_resolution_clock::now();

    elapsed = finish - start;
    std::cout << "pgMat<pgComplex<float>> For Elapsed time: " << elapsed.count() << " s\n";
    
    
    
    
    
    
    
    
    
    #ifdef _OPENACC
    acc_shutdown(acc_device_default);
    #endif





    return 0;
}