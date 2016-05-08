/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [fftGPU.hpp]

    Synopsis    [Wrappers to the cuFFT library supporting single and double
    				precision for GPU accelerated FFTs]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

// Based on https://www.olcf.ornl.gov/tutorials/mixing-openacc-with-gpu-libraries/


#ifndef PowerGrid_fftGPU_hpp
#define PowerGrid_fftGPU_hpp
#ifdef _OPENACC
#include <type_traits>
#include "cufft.h"
#include <complex>

//Like Armadillo, we're using SFINAE here to choose between float and double. (Maybe FP16 some day in the future)
//We need enable_if to choose which version to run based on the type of the template parameter.
template<typename T1, typename std::enable_if<std::is_same<T1, float>::value, int>::type = 0>
void ifft2dGPU(T1* d_data, int nx, int ny, void* stream)
{
	cufftHandle plan;

	cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL);

	if (cufftPlan2d(&plan, nx, ny, CUFFT_C2C)!=CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT error: Plan creation failed");
	}

	cufftSetStream(plan, (cudaStream_t) stream);
	cufftExecC2C(plan, (cufftComplex*) d_data, (cufftComplex*) d_data, CUFFT_INVERSE);
	cufftDestroy(plan);
}

template<typename T1, typename std::enable_if<std::is_same<T1, float>::value, int>::type = 0>
void fft2dGPU(T1* d_data, int nx, int ny, void* stream)
{
	cufftHandle plan;
	cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL);

	if (cufftPlan2d(&plan, nx, ny, CUFFT_C2C)!=CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT error: Plan creation failed");
	}

	cufftSetStream(plan, (cudaStream_t) stream);
	cufftExecC2C(plan, (cufftComplex*) d_data, (cufftComplex*) d_data, CUFFT_FORWARD);
	cufftDestroy(plan);
}

template<typename T1, typename std::enable_if<std::is_same<T1, float>::value, int>::type = 0>
void ifft3dGPU(T1* d_data, int nx, int ny, int nz, void* stream)
{
	cufftHandle plan;
	cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL);

	if (cufftPlan3d(&plan, nz, ny, nx, CUFFT_C2C)!=CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT error: Plan creation failed");
	}

	cufftSetStream(plan, (cudaStream_t) stream);
	cufftExecC2C(plan, (cufftComplex*) d_data, (cufftComplex*) d_data, CUFFT_INVERSE);
	cufftDestroy(plan);
}

template<typename T1, typename std::enable_if<std::is_same<T1, float>::value, int>::type = 0>
void fft3dGPU(T1* d_data, int nx, int ny, int nz, void* stream)
{
	cufftHandle plan;
	cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL);

	if (cufftPlan3d(&plan, nz, ny, nx, CUFFT_C2C)!=CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT error: Plan creation failed");
	}

	cufftSetStream(plan, (cudaStream_t) stream);
	cufftExecC2C(plan, (cufftComplex*) d_data, (cufftComplex*) d_data, CUFFT_FORWARD);
	cufftDestroy(plan);
}

template<typename T1, typename std::enable_if<std::is_same<T1, double>::value, int>::type = 0>
void ifft2dGPU(T1* d_data, int nx, int ny, void* stream)
{
	//printf("Running 2d inverse xform \n");
	cufftHandle plan;

	cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL);

	if (cufftPlan2d(&plan, ny, nx, CUFFT_Z2Z)!=CUFFT_SUCCESS) {
		printf( "CUFFT error: Plan creation failed\n");
	}
	//printf("Built plan \n");
	cufftSetStream(plan, (cudaStream_t) stream);
	if (cufftExecZ2Z(plan, (cufftDoubleComplex*) d_data, (cufftDoubleComplex*) d_data, CUFFT_INVERSE)!=CUFFT_SUCCESS) {
		printf("CUFFT error: Plan execution failed\n");
	};
	cufftDestroy(plan);
}

template<typename T1, typename std::enable_if<std::is_same<T1, double>::value, int>::type = 0>
void fft2dGPU(T1* d_data, int nx, int ny, void* stream)
{
	//printf("Running 2d forward xform \n");
	cufftHandle plan;

	cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL);

	if (cufftPlan2d(&plan, ny, nx, CUFFT_Z2Z)!=CUFFT_SUCCESS) {
		printf("CUFFT error: Plan creation failed\n");
	}
	//printf("Built plan \n");

	cufftSetStream(plan, (cudaStream_t) stream);
	if (cufftExecZ2Z(plan, (cufftDoubleComplex*) d_data, (cufftDoubleComplex*) d_data, CUFFT_FORWARD)!=CUFFT_SUCCESS) {
		printf("CUFFT error: Plan execution failed\n");
	};
	cufftDestroy(plan);
}

template<typename T1, typename std::enable_if<std::is_same<T1, double>::value, int>::type = 0>
void ifft3dGPU(T1*d_data, int nx, int ny, int nz, void* stream)
{
	//printf("Running 3d inverse xform \n");
	cufftHandle plan;
	cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL);

	if (cufftPlan3d(&plan, nz, ny, nx, CUFFT_Z2Z)!=CUFFT_SUCCESS) {
		printf("CUFFT error: Plan creation failed\n");
	}
	//printf("Built plan \n");

	cufftSetStream(plan, (cudaStream_t) stream);
	if (cufftExecZ2Z(plan, (cufftDoubleComplex*) d_data, (cufftDoubleComplex*) d_data, CUFFT_INVERSE)!=CUFFT_SUCCESS) {
		printf("CUFFT error: Plan execution failed\n");
	};
	cufftDestroy(plan);
}

template<typename T1, typename std::enable_if<std::is_same<T1, double>::value, int>::type = 0>
void fft3dGPU(T1* d_data, int nx, int ny, int nz, void* stream)
{   //printf("Running 3d forward xform \n");
	cufftHandle plan;
	cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_ALL);

	if (cufftPlan3d(&plan, nz, ny, nx, CUFFT_Z2Z)!=CUFFT_SUCCESS) {
		printf("CUFFT error: Plan creation failed\n");
	}
	//printf("Built plan \n");

	cufftSetStream(plan, (cudaStream_t) stream);

	if (cufftExecZ2Z(plan, (cufftDoubleComplex*) d_data, (cufftDoubleComplex*) d_data, CUFFT_FORWARD)!=CUFFT_SUCCESS) {
		printf("CUFFT error: Plan execution failed\n");
	};
	cufftDestroy(plan);
}
#endif //_OPENACC
#endif //PowerGrid_fftGPU_hpp