/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [DWIRecon.cpp]

    Synopsis    [Code implementing diffusion weighted image reconstructions
                    using MPI using the DFT object.]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC
                 0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/13/2016]

 *****************************************************************************/

#include "DWIRecon.h"

using namespace arma;

template <typename T1>
int DWIRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
             uword niter, uword nc) {
  string testPath = dataPath;
  typedef complex<T1> CxT1;
  // uword Nx =64;
  // uword Ny =64;
  // uword Nz =16;

  // Setup image space coordinates/trajectory
  Cube<T1> ix(Nx, Ny, Nz);
  Cube<T1> iy(Nx, Ny, Nz);
  Cube<T1> iz(Nx, Ny, Nz);
  Col<T1> FM;
  Col<CxT1> SMap;
  Col<T1> PMap;

  ix.zeros();
  iy.zeros();
  iz.zeros();

  // generate the image space coordinates of the voxels we want to reconstruct
  // after vectorizing ix and iy the image coordinates must match the Field and
  // SENSe map image coordinates
  for (uword ii = 0; ii < Ny; ii++) {     // y
    for (uword jj = 0; jj < Nx; jj++) {   // x
      for (uword kk = 0; kk < Nz; kk++) { // z

        ix(ii, jj, kk) = ((T1)jj - (T1)Nx / 2.0) / ((T1)Nx);
        iy(ii, jj, kk) = ((T1)ii - (T1)Ny / 2.0) / ((T1)Ny);
        iz(ii, jj, kk) = ((T1)kk - (T1)Nz / 2.0) / ((T1)Nz);
      }
    }
  }

  Col<T1> kx;
  loadmat(testPath + "kx.mat", "kx", &kx);
  Col<T1> ky;
  loadmat(testPath + "ky.mat", "ky", &ky);
  Col<T1> kz;
  loadmat(testPath + "kz.mat", "kz", &kz);

  uword nro;
  nro = kx.n_elem;

  Col<T1> tvec;
  loadmat(testPath + "t.mat", "t", &tvec);

  loadmat(testPath + "FM.mat", "FM", &FM);
  // FM.zeros();

  loadmat(testPath + "SMap.mat", "SMap", &SMap);

  // nonlinear motion correction for DWI
  loadmat(testPath + "PMap.mat", "PMap", &PMap);
  // The navigator phase is being referenced to zero
  pcSENSE<T1> S_DWI(kx, ky, kz, Nx, Ny, Nz, nc, tvec, SMap, vectorise(FM),
                    0 - PMap);

  cout << "loading data" << endl;
  Col<CxT1> data;
  data = 1E-6 * data;
  loadmat(testPath + "data.mat", "data", &data);

  // Variables needed for the recon: Penalty object, num of iterations

  cout << "Iniitalizing QuadPenalty" << endl;
  QuadPenalty<T1> R(Nx, Ny, Nz, 0.0);
  cout << "QuadPenalty setup successfull" << endl;

  // uword niter = 10;
  Col<CxT1> xinit(Nx * Ny * Nz); // initial estimate of x
  xinit.zeros();
  Col<T1> W;
  // W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is
  // it right?
  // W=1.0;
  W.ones(nro * nc);

  cout << "Runing pwls with Gnufft" << endl;
  Col<CxT1> test_DWI_img;
  test_DWI_img = solve_pwls_pcg<T1, pcSENSE<T1>, QuadPenalty<T1>>(
      xinit, S_DWI, W, data, R, niter);
  savemat(testPath + "test_DWICGMC.mat", "img", test_DWI_img);

  return 0;
}

#ifdef PowerGridMPI
// boost::mpi version of DWIDft reconstruction routine.
template <typename T1>
int DWIRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
             uword niter, uword nc, bmpi::environment &env,
             bmpi::communicator &world) {
  string testPath = dataPath;
  typedef complex<T1> CxT1;

  // Setup image space coordinates/trajectory
  Cube<T1> ix(Nx, Ny, Nz);
  Cube<T1> iy(Nx, Ny, Nz);
  Cube<T1> iz(Nx, Ny, Nz);
  Col<T1> FM;
  Col<CxT1> SMap;
  Col<T1> PMap;

  ix.zeros();
  iy.zeros();
  iz.zeros();

  // generate the image space coordinates of the voxels we want to reconstruct
  // after vectorizing ix and iy the image coordinates must match the Field and
  // SENSe map image coordinates
  for (uword ii = 0; ii < Ny; ii++) {     // y
    for (uword jj = 0; jj < Nx; jj++) {   // x
      for (uword kk = 0; kk < Nz; kk++) { // z

        ix(ii, jj, kk) = ((T1)jj - (T1)Nx / 2.0) / ((T1)Nx);
        iy(ii, jj, kk) = ((T1)ii - (T1)Ny / 2.0) / ((T1)Ny);
        iz(ii, jj, kk) = ((T1)kk - (T1)Nz / 2.0) / ((T1)Nz);
      }
    }
  }

  Col<T1> kx;
  loadmat(testPath + "kx.mat", "kx", &kx);
  Col<T1> ky;
  loadmat(testPath + "ky.mat", "ky", &ky);
  Col<T1> kz;
  loadmat(testPath + "kz.mat", "kz", &kz);

  uword nro;
  nro = kx.n_elem;

  Col<T1> tvec;
  loadmat(testPath + "t.mat", "t", &tvec);

  loadmat(testPath + "FM.mat", "FM", &FM);

  loadmat(testPath + "SMap.mat", "SMap", &SMap);

  // nonlinear motion correction for DWI
  loadmat(testPath + "PMap.mat", "PMap", &PMap);
  // The navigator phase is being referenced to zero
  mpipcSENSE<T1> S_DWI(kx, ky, kz, Nx, Ny, Nz, nc, tvec, SMap, vectorise(FM),
                       0 - PMap, env, world);

  cout << "loading data" << endl;
  Col<CxT1> data;
  data = 1E-6 * data;
  loadmat(testPath + "data.mat", "data", &data);

  // Variables needed for the recon: Penalty object, num of iterations

  cout << "Iniitalizing QuadPenalty" << endl;
  // TVPenalty<T1>R(Nx,Ny,Nz,100000000.0,.000001);
  // TVPenalty<T1>R(Nx,Ny,Nz,0.0,.01);
  QuadPenalty<T1> R(Nx, Ny, Nz, 0.0);
  cout << "QuadPenalty setup successfull" << endl;

  // uword niter = 10;
  Col<CxT1> xinit(Nx * Ny * Nz); // initial estimate of x
  xinit.zeros();
  Col<T1> W;
  // W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is
  // it right?
  // W=1.0;
  W.ones(nro * nc);

  cout << "Runing pwls with Gnufft" << endl;
  Col<CxT1> test_DWI_img;
  test_DWI_img = solve_pwls_pcg<T1, mpipcSENSE<T1>, QuadPenalty<T1>>(
      xinit, S_DWI, W, data, R, niter);
  if (world.rank() == 0) {
    savemat(testPath + "test_DWICGMC.mat", "img", test_DWI_img);
  }
  return 0;
}

template <typename T1>
int DWIRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
             uword niter, uword nc, uword nimage, uword nphase, uword nslab,
             bmpi::environment &env, bmpi::communicator &world) {
  string testPath = dataPath;
  typedef complex<T1> CxT1;
  // uword Nx =64;
  // uword Ny =64;
  // uword Nz =16;

  // Shift to index one for files
  nslab++;
  nphase++;
  nimage += 2;

  // Setup image space coordinates/trajectory
  Cube<T1> ix(Nx, Ny, Nz);
  Cube<T1> iy(Nx, Ny, Nz);
  Cube<T1> iz(Nx, Ny, Nz);
  Col<T1> FM;
  Col<CxT1> SMap;
  Col<T1> PMap;

  ix.zeros();
  iy.zeros();
  iz.zeros();

  // generate the image space coordinates of the voxels we want to reconstruct
  // after vectorizing ix and iy the image coordinates must match the Field and
  // SENSe map image coordinates
  for (uword ii = 0; ii < Ny; ii++) {     // y
    for (uword jj = 0; jj < Nx; jj++) {   // x
      for (uword kk = 0; kk < Nz; kk++) { // z

        ix(ii, jj, kk) = ((T1)jj - (T1)Nx / 2.0) / ((T1)Nx);
        iy(ii, jj, kk) = ((T1)ii - (T1)Ny / 2.0) / ((T1)Ny);
        iz(ii, jj, kk) = ((T1)kk - (T1)Nz / 2.0) / ((T1)Nz);
      }
    }
  }

  Col<T1> kx;
  loadmat(testPath + "kx.mat", "kx", &kx);
  Col<T1> ky;
  loadmat(testPath + "ky.mat", "ky", &ky);
  Col<T1> kz;
  loadmat(testPath + "kz.mat", "kz", &kz);

  uword nro;
  nro = kx.n_elem;

  Col<T1> tvec;
  loadmat(testPath + "t.mat", "t", &tvec);

  // uword L = 1;
  loadmat(testPath + "FM_" + std::to_string(nslab) + ".mat", "FM", &FM);
  // FM.zeros();

  // uword nc = 4;
  loadmat(testPath + "sen_" + std::to_string(nslab) + ".mat", "sen", &SMap);

  // nonlinear motion correction for DWI
  loadmat(testPath + "phaseMap_" + std::to_string(nimage) + "_" +
              std::to_string(nslab) + "_" + std::to_string(nphase) + ".mat",
          "PMap", &PMap);
  // PMap.zeros();
  // The navigator phase is being referenced to zero
  mpipcSENSE<T1> S_DWI(kx, ky, kz, Nx, Ny, Nz, nc, tvec, SMap, vectorise(FM),
                       0 - PMap, env, world);

  cout << "loading data" << endl;
  Col<CxT1> data;
  data = 1E-7 * data;
  loadmat(testPath + "data_" + std::to_string(nimage) + "_" +
              std::to_string(nslab) + "_" + std::to_string(nphase) + ".mat",
          "data", &data);

  // Variables needed for the recon: Penalty object, num of iterations

  cout << "Iniitalizing QuadPenalty" << endl;
  // TVPenalty<T1>R(Nx,Ny,Nz,100000000.0,.000001);
  // TVPenalty<T1>R(Nx,Ny,Nz,0.0,.01);
  QuadPenalty<T1> R(Nx, Ny, Nz, 0.0);
  cout << "QuadPenalty setup successfull" << endl;

  // uword niter = 10;
  Col<CxT1> xinit(Nx * Ny * Nz); // initial estimate of x
  xinit.zeros();
  Col<T1> W;
  // W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is
  // it right?
  // W=1.0;
  W.ones(nro * nc);

  cout << "Runing pwls with Gnufft" << endl;
  Col<CxT1> test_DWI_img;
  test_DWI_img = solve_pwls_pcg<T1, mpipcSENSE<T1>, QuadPenalty<T1>>(
      xinit, S_DWI, W, data, R, niter);
  savemat(testPath + "PGout_" + std::to_string(nimage) + "_" +
              std::to_string(nslab) + "_" + std::to_string(nphase) + ".mat",
          "img", test_DWI_img);

  return 0;
}
#endif // PowerGridMPI
template <typename T1>
int DWIRecon(string dataPath, uword Nx, uword Ny, uword Nz, uword L,
             uword niter, uword nc, uword nimage, uword nphase, uword nslab) {
  string testPath = dataPath;
  typedef complex<T1> CxT1;
  // uword Nx =64;
  // uword Ny =64;
  // uword Nz =16;

  // Shift to index one for files
  nslab++;
  nphase++;
  nimage += 2;

  // Setup image space coordinates/trajectory
  Cube<T1> ix(Nx, Ny, Nz);
  Cube<T1> iy(Nx, Ny, Nz);
  Cube<T1> iz(Nx, Ny, Nz);
  Col<T1> FM;
  Col<CxT1> SMap;
  Col<T1> PMap;

  ix.zeros();
  iy.zeros();
  iz.zeros();

  // generate the image space coordinates of the voxels we want to reconstruct
  // after vectorizing ix and iy the image coordinates must match the Field and
  // SENSe map image coordinates
  for (uword ii = 0; ii < Ny; ii++) {     // y
    for (uword jj = 0; jj < Nx; jj++) {   // x
      for (uword kk = 0; kk < Nz; kk++) { // z

        ix(ii, jj, kk) = ((T1)jj - (T1)Nx / 2.0) / ((T1)Nx);
        iy(ii, jj, kk) = ((T1)ii - (T1)Ny / 2.0) / ((T1)Ny);
        iz(ii, jj, kk) = ((T1)kk - (T1)Nz / 2.0) / ((T1)Nz);
      }
    }
  }

  Col<T1> kx;
  loadmat(testPath + "kx.mat", "kx", &kx);
  Col<T1> ky;
  loadmat(testPath + "ky.mat", "ky", &ky);
  Col<T1> kz;
  loadmat(testPath + "kz.mat", "kz", &kz);

  uword nro;
  nro = kx.n_elem;

  Col<T1> tvec;
  loadmat(testPath + "t.mat", "t", &tvec);

  // uword L = 1;
  loadmat(testPath + "FM_" + std::to_string(nslab) + ".mat", "FM", &FM);
  std::cout << testPath + "FM_" + std::to_string(nslab) + ".mat" << std::endl;
  // FM.zeros();

  // uword nc = 4;
  loadmat(testPath + "sen_" + std::to_string(nslab) + ".mat", "sen", &SMap);

  // nonlinear motion correction for DWI
  loadmat(testPath + "phaseMap_" + std::to_string(nimage) + "_" +
              std::to_string(nslab) + "_" + std::to_string(nphase) + ".mat",
          "PMap", &PMap);
  // PMap.zeros();
  // The navigator phase is being referenced to zero
  pcSENSE<T1> S_DWI(kx, ky, kz, Nx, Ny, Nz, nc, tvec, SMap, vectorise(FM),
                    0 - PMap);

  cout << "loading data" << endl;
  Col<CxT1> data;
  data = 1E-9 * data;
  loadmat(testPath + "data_" + std::to_string(nimage) + "_" +
              std::to_string(nslab) + "_" + std::to_string(nphase) + ".mat",
          "data", &data);

  // Variables needed for the recon: Penalty object, num of iterations

  cout << "Iniitalizing QuadPenalty" << endl;
  // TVPenalty<T1>R(Nx,Ny,Nz,100000000.0,.000001);
  // TVPenalty<T1>R(Nx,Ny,Nz,0.0,.01);
  QuadPenalty<T1> R(Nx, Ny, Nz, 0.0);
  cout << "QuadPenalty setup successfull" << endl;

  // uword niter = 10;
  Col<CxT1> xinit(Nx * Ny * Nz); // initial estimate of x
  xinit.zeros();
  Col<T1> W;
  // W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is
  // it right?
  // W=1.0;
  W.ones(nro * nc);

  cout << "Runing pwls with Gnufft" << endl;
  Col<CxT1> test_DWI_img;
  test_DWI_img = solve_pwls_pcg<T1, pcSENSE<T1>, QuadPenalty<T1>>(
      xinit, S_DWI, W, data, R, niter);
  savemat(testPath + "PGout_" + std::to_string(nimage) + "_" +
              std::to_string(nslab) + "_" + std::to_string(nphase) + ".mat",
          "img", test_DWI_img);

  return 0;
}
