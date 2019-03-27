/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [pcSenseTimeSeg.cpp]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Joseph Holtrop, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/
#include "pcSenseTimeSeg.h"

template <typename T1> pcSenseTimeSeg<T1>::~pcSenseTimeSeg() {

	for (uword jj = 0; jj < Ns; jj++) {
		delete G[jj];
		delete AObj[jj];
	}
	delete[] G;
	delete[] AObj;

}

// Class constructor
template <typename T1>
pcSenseTimeSeg<T1>::pcSenseTimeSeg(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx, uword ny,
                     uword nz, uword nc, Col<T1> t, uword lSegs, uword intType, Col<complex<T1>> SENSEmap,
                     Col<T1> FieldMap, Col<T1> ShotPhaseMap) {
        Ni = nx * ny * nz;
        Nc = nc;
        Ns = ShotPhaseMap.n_elem / Ni;
        Nd = kx.n_elem / Ns;
        std::cout << "Nd = " << Nd << std::endl;
        std::cout << "Ns = " << Ns << std::endl;
        std::cout << "Nc = " << Nc << std::endl;
        std::cout << "Ni = " << Ni << std::endl;
        SMap = reshape(SENSEmap, Ni, Nc);
        PMap = reshape(ShotPhaseMap, Ni, Ns);
        FMap = FieldMap;
        Kx = reshape(kx, Nd, Ns);
        Ky = reshape(ky, Nd, Ns);
        Kz = reshape(kz, Nd, Ns);
        Nx = nx;
        Ny = ny;
        Nz = nz;
        Tvec = reshape(t, Nd, Ns);
        L = lSegs;
        type = intType;
        Cube<T1> ix;
        ix.zeros(Nx, Ny, Nz);
        Cube<T1> iy;
        iy.zeros(Nx, Ny, Nz);
        Cube<T1> iz;
        iz.zeros(Nx, Ny, Nz);

        // generate the image space coordinates of the voxels we want to reconstruct
        // after vectorizing ix and iy the image coordinates must match the Field and
        // SENSe map image coordinates
        for (uword ii = 0; ii < Ny; ii++) { // y
                for (uword jj = 0; jj < Nx; jj++) { // x
                        for (uword kk = 0; kk < Nz; kk++) { // z
                                ix(ii, jj, kk) = ((T1)jj - (T1)Nx / 2.0) / ((T1)Nx);
                                iy(ii, jj, kk) = ((T1)ii - (T1)Ny / 2.0) / ((T1)Ny);
                                iz(ii, jj, kk) = ((T1)kk - (T1)Nz / 2.0) / ((T1)Nz);
                        }
                }
        }

        Ix = vectorise(ix);
        Iy = vectorise(iy);
        Iz = vectorise(iz);

        //AObj = new Gdft<T1> *[Ns];
		G = new Gnufft<T1> *[Ns];
		AObj = new TimeSegmentation <T1, Gnufft<T1>> *[Ns];

        // Initialize the field correction and G objects we need for this
        // reconstruction
        for (uword jj = 0; jj < Ns; jj++) {

                //AObj[jj] =
                //        new Gdft<T1>(Nd, Nx * Ny * Nz, Kx.col(jj), Ky.col(jj), Kz.col(jj), Ix,
                //                     Iy, Iz, vectorise(FMap), vectorise(Tvec.col(jj)));
	        G[jj] = new Gnufft<T1>(Nd, 2.0, Nx, Ny, Nz, Kx.col(jj), Ky.col(jj), Kz.col(jj), Ix, Iy, Iz);
                AObj[jj] = new TimeSegmentation <T1, Gnufft<T1>>(*G[jj], vectorise(FMap),
                 vectorise(Tvec.col(jj)), (uword) Nd, (uword)(Nx * Ny * Nz), (uword) L, type, Ns);
        }
}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1>
Col<complex<T1> > pcSenseTimeSeg<T1>::operator*(const Col<complex<T1> > &d) const {
        RANGE("pcSenseTimeSeg::operator*")
        Mat<complex<T1> > outData = zeros<Mat<complex<T1> > >(Nd, Ns * Nc);
        Mat<complex<T1> > expiPMap = exp(-i * PMap);
        Mat<T1> temp2;
        // Coil loop. Each coil exists for each shot, so we need to work with these.
        for (unsigned int ii = 0; ii < Nc; ii++) {

                /*
                   Gnufft <T1>*  G = new Gnufft<T1>(Nd, 2.0, Nx, Ny, Nz, Kx.col(jj),
                   Ky.col(jj), Kz.col(jj), Ix, Iy, Iz);
                   TimeSegmentation<T1,Gnufft<T1>>*  AObj = new
                   TimeSegmentation<T1,Gnufft<T1>>(*G, vectorise(FMap),
                   vectorise(Tvec.col(jj)), (uword)Nd, (uword)(Nx * Ny * Nz), (uword)L,
                                                     type , (uword) 1);
                 */
                /*
                   std::cout << "A address = " << &AObj << std::endl;
                   std::cout << "Nd = " << Nd << std::endl;
                   std::cout << "NxNyNz = " << (uword)(Nx * Ny * Nz) << std::endl;
                   std::cout << "A.n1 = " << AObj->n1 << std::endl;
                   std::cout << "A.n2 = " << AObj->n2 << std::endl;
                   std::cout << "A.L = " << AObj->L << std::endl;
                   std::cout << "A.Nshots = " << AObj->Nshots << std::endl;
                 */
                // Shot loop. Each shot has its own kspace trajectory
                for (unsigned int jj = 0; jj < Ns; jj++) {
                        outData.col(jj + ii * Ns) =
                                (*AObj[jj]) * (d % (SMap.col(ii) % expiPMap.col(jj)));
                        //std::cout << "Processed shot # " << jj << " coil # " << ii << std::endl;
                }
                // delete AObj;
                // delete G;
        }
        // equivalent to returning col(output) in MATLAB with IRT
        return vectorise(outData);
}

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
template <typename T1>
Col<complex<T1> > pcSenseTimeSeg<T1>::operator/(const Col<complex<T1> > &d) const {
        RANGE("pcSenseTimeSeg::operator/");
        Mat<complex<T1> > inData = reshape(d, Nd, Ns * Nc);
        Mat<complex<T1> > expiPMap = conj(exp(-i * PMap));
        Mat<complex<T1> > conjSMap = conj(SMap);
        Col<complex<T1> > outData = zeros<Col<complex<T1> > >(Ni);
        // Coil Loop - for each shot we have a full set of coil data.
        for (unsigned int ii = 0; ii < Nc; ii++) {
        
                // Use grid or DFT?
                // Gdft<T1> G(Nd, Ni, Kx.col(jj), Ky.col(jj), Kz.col(jj), Ix, Iy, Iz,FMap,
                // vectorise(Tvec.col(jj)));
                /*
                   Gnufft <T1>* G = new Gnufft<T1>(Nd, 2.0, Nx, Ny, Nz, Kx.col(jj), Ky.col(jj),
                   Kz.col(jj), Ix, Iy, Iz);
                   TimeSegmentation <T1, Gnufft<T1>>* AObj = new
                   TimeSegmentation<T1,Gnufft<T1>>(*G, vectorise(FMap),
                   vectorise(Tvec.col(jj)), Nd, Nx * Ny * Nz, L, type);
                 */

                // Shot Loop. Each shot has it's own k-space trajectory
                for (unsigned int jj = 0; jj < Ns; jj++) {
                        //outData += conj(SMap.col(ii) % exp(-i * (PMap.col(jj)))) %
                        //           ((*AObj[jj]) / inData.col(jj + ii * Ns));
                        outData += (conjSMap.col(ii) % expiPMap.col(jj)) %
                                   ((*AObj[jj]) / inData.col(jj + ii * Ns));
                        //std::cout << "Processed shot # " << jj << " coil # " << ii << std::endl;
                }
                // delete AObj;
                // delete G;
        }

        // equivalent to returning col(output) in MATLAB with IRT
        return vectorise(outData);
}

// Explicit Instantiation
template class pcSenseTimeSeg<float>;
template class pcSenseTimeSeg<double>;
