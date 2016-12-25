/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [mpipcSENSE.cpp]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/
#include "mpipcSENSE.h"

template <typename T1> mpipcSENSE<T1>::~mpipcSENSE() {

  for (uword jj = 0; jj < (*taskList)[world->rank()].size(); jj++) {
    delete AObj[jj];
  }
  std::cout << "Deleted Objects in AObj array" << std::endl;
  delete[] AObj;
  std::cout << "Deleted Objects AObj array" << std::endl;
  delete taskList;
  std::cout << "Deleted taskList array" << std::endl;
}

template <typename T1>
mpipcSENSE<T1>::mpipcSENSE(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx,
                           uword ny, uword nz, uword nc, Col<T1> t,
                           Col<CxT1> SENSEmap, Col<T1> FieldMap,
                           Col<T1> ShotPhaseMap, bmpi::environment &en,
                           bmpi::communicator &wor) {
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
  env = &en;
  world = &wor;

  Cube<T1> ix;
  ix.zeros(Nx, Ny, Nz);
  Cube<T1> iy;
  iy.zeros(Nx, Ny, Nz);
  Cube<T1> iz;
  iz.zeros(Nx, Ny, Nz);

  // generate the image space coordinates of the voxels we want to reconstruct
  // after vectorizing ix and iy the image coordinates must match the Field
  // and SENSe map image coordinates
  for (uword ii = 0; ii < Ny; ii++) {     // y
    for (uword jj = 0; jj < Nx; jj++) {   // x
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

  // Deal with task allocation to the MPI ranks
  // This gives us the rank (index) of the process in the collective MPI
  // space.
  // The process will only calculate shot jj of the acquisition.
  uword ShotRank = world->rank();
  uword numShots = Ns / world->size();

  shotList.resize(Ns * Nc);
  coilList.resize(Ns * Nc);
  // Need to generate a mapping of task number to coils and to shots
  // separately.
  for (uword ii = 0; ii < Ns; ii++) {   // Shot loop
    for (uword jj = 0; jj < Nc; jj++) { // Coil Loop
      shotList(ii * Nc + jj) = ii;
      coilList(ii * Nc + jj) = jj;
    }
  }

  // Now all we have to do is generate a mapping of tasks to MPI ranks.

  // Start by creating a 2D vector (vector of vectors) to hold the MPI rank ->
  // task mapping
  vector<uword> init(0);
  taskList = new std::vector<std::vector<uword>>(world->size(), init);
  uword remaining = Ns * Nc;
  uword process = 0;
  for (uword task = 0; task < Ns * Nc; task++) {
    (*taskList)[process].push_back(task);
    process++;
    if (process == world->size()) {
      process = 0;
    }
  }

  AObj = new Gdft<T1> *[(*taskList)[world->rank()].size()];
  uword taskIndex;
  uword coilIndex;
  uword shotIndex;
  // Initialize the field correction and G objects we need for this
  // reconstruction
  for (uword jj = 0; jj < (*taskList)[world->rank()].size(); jj++) {
    taskIndex = (*taskList)[world->rank()].at(jj);
    shotIndex = shotList(taskIndex);
    AObj[jj] = new Gdft<T1>(Nd, Nx * Ny * Nz, Kx.col(shotIndex),
                            Ky.col(shotIndex), Kz.col(shotIndex), Ix, Iy, Iz,
                            vectorise(FMap), vectorise(Tvec.col(shotIndex)));
  }
}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1>
Col<CxT1> mpipcSENSE<T1>::operator*(const Col<CxT1> &d) const {
  // uword ShotRank = world->rank();
  Mat<CxT1> outData = zeros<Mat<CxT1>>(Nd, Ns * Nc);
  Mat<CxT1> tempOutData =
      zeros<Mat<CxT1>>(Nd, (*taskList)[world->rank()].size());
  Mat<CxT1> tempOutData2 =
      zeros<Mat<CxT1>>(Nd, (*taskList)[world->rank()].size());
  uword taskIndex;
  uword coilIndex;
  uword shotIndex;
  // Shot loop. Each shot has it's own kspace trajectory
  // for (unsigned int jj = 0; jj < Ns; jj++) {
  std::cout << "Task Number = " << (*taskList)[world->rank()].size()
            << std::endl;
  // Coil loop. Each coil exists for each shot, so we need to work with these.
  for (uword jj = 0; jj < (*taskList)[world->rank()].size(); jj++) {
    taskIndex = (*taskList)[world->rank()].at(jj);
    shotIndex = shotList(taskIndex);
    coilIndex = coilList(taskIndex);
    tempOutData.col(jj) =
        (*AObj[jj]) *
        (d % (SMap.col(coilIndex) % exp(-i * (PMap.col(shotIndex)))));
  }

  //}
  // Now let's do some MPI stuff here.
  if (world->rank() == 0) {
    std::vector<Mat<CxT1>> OutDataGather(world->size());
    // Collect all the data into OutDataGather an std::vector collective
    // std::cout << "Rank #: " << world->rank() << " reached foward xform
    // gather" << std::endl;
    bmpi::gather<Mat<CxT1>>(*world, tempOutData, OutDataGather, 0);
    // std::cout << "Rank #: " << world->rank() << " passed forward xform
    // gather" << std::endl;
    for (uword jj = 0; jj < world->size(); jj++) {
      // std::cout << "About to access element #" << jj << " in gathered data"
      // << std::endl;
      // std::cout << "Size of returned stl::vector<> = " <<
      // OutDataGather.size() << std::endl;
      tempOutData2 = OutDataGather[jj];
      for (uword ii = 0; ii < (*taskList)[jj].size(); ii++) {
        taskIndex = (*taskList)[jj].at(ii);
        shotIndex = shotList(taskIndex);
        coilIndex = coilList(taskIndex);
        outData.col(shotIndex + coilIndex * Ns) = tempOutData2.col(ii);
      }
    }

  } else {
    std::cout << "Rank #: " << world->rank() << " reached foward xform gather"
              << std::endl;
    bmpi::gather<Mat<CxT1>>(*world, tempOutData, 0);
  }
  std::cout << "Rank #: " << world->rank() << " reached foward xform broadcast"
            << std::endl;
  bmpi::broadcast(*world, outData, 0);
  std::cout << "Rank #: " << world->rank() << " passed foward xform broadcast"
            << std::endl;

  // equivalent to returning col(output) in MATLAB with IRT
  return vectorise(outData);
}

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
template <typename T1>
Col<CxT1> mpipcSENSE<T1>::operator/(const Col<CxT1> &d) const {
  // uword ShotRank = world->rank();
  Mat<CxT1> inData = reshape(d, Nd, Ns * Nc);
  uword taskIndex;
  uword coilIndex;
  uword shotIndex;
  // Col <CxT1> tempOutData = zeros < Col < CxT1 >> (Ni);
  Col<CxT1> outData = zeros<Col<CxT1>>(Ni);
  // Shot Loop. Each shot has it's own k-space trajectory
  // for (unsigned int jj = 0; jj < Ns; jj++) {
  // Coil Loop - for each shot we have a full set of coil data.
  for (uword jj = 0; jj < (*taskList)[world->rank()].size(); jj++) {
    taskIndex = (*taskList)[world->rank()].at(jj);
    shotIndex = shotList(taskIndex);
    coilIndex = coilList(taskIndex);
    outData += conj(SMap.col(coilIndex) % exp(-i * (PMap.col(shotIndex)))) %
               ((*AObj[jj]) / inData.col(shotIndex + coilIndex * Ns));
    // std::cout << "Processed shot # " << jj << " coil # " << ii <<
    // std::endl;
  }

  if (world->rank() == 0) {
    std::vector<Col<CxT1>> OutDataGather;
    // Collect all the data into OutDataGather an std::vector collective
    bmpi::gather<Col<CxT1>>(*world, outData, OutDataGather, 0);
    outData.zeros(Ni);
    for (uword jj = 0; jj < world->size(); jj++) {
      // Skip the zeroth rank because that is the outData we started with on
      // this processor.
      // Now we will manually reduce the data by adding across the elements.
      outData += OutDataGather.at(jj);
    }
  } else {
    bmpi::gather<Col<CxT1>>(*world, outData, 0);
  }
  // Broadcast will send the data to all machines if the rank==0 and receive
  // the broadcasted data from rank=0 to
  // all other machines in the communicator, overwriting the outData value
  // already there.
  bmpi::broadcast<Col<CxT1>>(*world, outData, 0);
  // equivalent to returning col(output) in MATLAB with IRT
  return vectorise(outData);
}
