/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [GnufftRecon.hpp]

    Synopsis    [Implements a basic 2D and 3D recon with NUFFT and Time
    				Segmentation for field correction.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/26/2016]

 *****************************************************************************/

#ifndef __PowerGrid__GnufftRecon
#define __PowerGrid__GnufftRecon

#include "PowerGrid.h"

using namespace arma;


template<typename T1>
int GnufftRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L, uword niter, uword nc, uword nshots,
		T1 beta)
{
	string testPath = dataPath;
	typedef complex <T1> CxT1;

	//Setup image space coordinates/trajectory
	Cube <T1> ix(Nx, Ny, Nz);
	Cube <T1> iy(Nx, Ny, Nz);
	Cube <T1> iz(Nx, Ny, Nz);
	Col <T1> FM;
	Col <CxT1> SMap;

	ix.zeros();
	iy.zeros();
	iz.zeros();

	//generate the image space coordinates of the voxels we want to reconstruct
	// after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
	for (uword ii = 0; ii<Ny; ii++) { //y
		for (uword jj = 0; jj<Nx; jj++) { //x
			for (uword kk = 0; kk<Nz; kk++) { //z

				ix(ii, jj, kk) = ((T1) jj-(T1) Nx/2.0)/((T1) Nx);
				iy(ii, jj, kk) = ((T1) ii-(T1) Ny/2.0)/((T1) Ny);
				iz(ii, jj, kk) = ((T1) kk-(T1) Nz/2.0)/((T1) Nz);
			}
		}
	}

	Col <T1> kx;
	loadmat(testPath+"kx.mat", "kx", &kx);

	Col <T1> ky;
	loadmat(testPath+"ky.mat", "ky", &ky);

	Col <T1> kz;
	loadmat(testPath+"kz.mat", "kz", &kz);



	uword nro;
	nro = kx.n_elem;

	Col <T1> tvec;
	loadmat(testPath+"t.mat", "t", &tvec);

	loadmat(testPath+"FM.mat", "FM", &FM);


	// Fourier transform operator
	Gnufft<T1> Gg(nro, 2.0, Nx, Ny, Nz, kx, ky, kz, vectorise(ix), vectorise(iy), vectorise(iz));

	// Field correction operation
	uword type = 1; // 2 for min max time seg and 1 for Hanning
	TimeSegmentation<T1, Gnufft<T1>> A(Gg, vectorise(FM), vectorise(tvec), nro, Nx*Ny*Nz, L, type, nshots);

	loadmat(testPath+"SMap.mat", "SMap", &SMap);


	// Sense operator
	SENSE<T1, TimeSegmentation<T1, Gnufft<T1>>> Sg(A, SMap, nro, Nx*Ny*Nz, nc);


	Col <CxT1> data;
	loadmat(testPath+"data.mat", "data", &data);


	// Variables needed for the recon: Penalty object, num of iterations

	QuadPenalty<T1> R(Nx, Ny, Nz, beta);

	//uword niter = 10;
	Col <CxT1> xinit(Nx*Ny*Nz); // initial estimate of x
	xinit.zeros();
	Col <T1> W;
	W.ones(nro*nc);


	Col <CxT1> test_pwls;
	test_pwls = solve_pwls_pcg<T1, SENSE<T1, TimeSegmentation<T1, Gnufft<T1>>>, QuadPenalty<T1>>(xinit, Sg, W, data, R,
			niter);
	savemat(testPath+"test_pwls.mat", "img", test_pwls);

	return 0;
}

#endif //__PowerGrid__GnufftRecon
