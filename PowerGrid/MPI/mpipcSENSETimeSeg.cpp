/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [mpipcSENSETimeSeg.cpp]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/
#include "mpipcSENSETimeSeg.h"

template <typename T1> mpipcSENSETimeSeg<T1>::~mpipcSENSETimeSeg() {

  for (uword jj = 0; jj < (*taskList)[pWorld->rank()].size(); jj++) {
    delete GObj[jj];
    delete AObj[jj];
  }
  std::cout << "Deleted Objects in AObj array" << std::endl;
  delete[] AObj;
  delete[] GObj;
  std::cout << "Deleted Objects AObj array" << std::endl;
  delete taskList;
  std::cout << "Deleted taskList array" << std::endl;
}

template <typename T1>
mpipcSENSETimeSeg<T1>::mpipcSENSETimeSeg(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx,
                           uword ny, uword nz, uword nc, Col<T1> t, uword l, uword type_interp,
                           Col<complex<T1>> SENSEmap, Col<T1> FieldMap,
                           Col<T1> ShotPhaseMap, bmpi::environment &en,
                           bmpi::communicator &wor) {
  pEnv = &en;
  pWorld = &wor;
  Ni = nx * ny * nz;
  Nc = nc;
  Ns = ShotPhaseMap.n_elem / Ni;
  Nd = kx.n_elem / Ns;
  L = l;
  type = type_interp; 
  if (pWorld->rank() == 0) {
    std::cout << "Nd = " << Nd << std::endl;
    std::cout << "Ns = " << Ns << std::endl;
    std::cout << "Nc = " << Nc << std::endl;
    std::cout << "Ni = " << Ni << std::endl;
  }
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
  uword ShotRank = pWorld->rank();
  uword numShots = Ns / pWorld->size();

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
  taskList = new std::vector<std::vector<uword>>(pWorld->size(), init);
  //uword remaining = Ns * Nc;
  uword process = 0;
  for (uword task = 0; task < Ns * Nc; task++) {
    (*taskList)[process].push_back(task);
    process++;
    if (process == pWorld->size()) {
      process = 0;
    }
  }

  GObj = new Gnufft<T1> *[(*taskList)[pWorld->rank()].size()];
  AObj = new TimeSegmentation<T1, Gnufft<T1>> *[(*taskList)[pWorld->rank()].size()];

  uword taskIndex;
  //uword coilIndex;
  uword shotIndex;
  // Initialize the field correction and G objects we need for this
  // reconstruction
  for (uword jj = 0; jj < (*taskList)[pWorld->rank()].size(); jj++) {
    taskIndex = (*taskList)[pWorld->rank()].at(jj);
    shotIndex = shotList(taskIndex);
    
    GObj[jj] = new Gnufft<T1>(Nd, 2.0, Nx, Ny, Nz, Kx.col(shotIndex), Ky.col(shotIndex), Kz.col(shotIndex), Ix, Iy, Iz);
    AObj[jj] = new TimeSegmentation <T1, Gnufft<T1>>(*GObj[jj], vectorise(FMap),
                                                      vectorise(Tvec.col(shotIndex)), (uword) Nd, 
                                                      (uword)(Nx * Ny * Nz), (uword) L, type, Ns);
  }

}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1>
Col<complex<T1>> mpipcSENSETimeSeg<T1>::operator*(const Col<complex<T1>> &d) const {
  // uword ShotRank = pWorld->rank();
  Mat<complex<T1>> outData = zeros<Mat<complex<T1>>>(Nd, Ns * Nc);
  Mat<complex<T1>> tempOutData =
      zeros<Mat<complex<T1>>>(Nd, (*taskList)[pWorld->rank()].size());
  Mat<complex<T1>> tempOutData2 =
      zeros<Mat<complex<T1>>>(Nd, (*taskList)[pWorld->rank()].size());
  uword taskIndex;
  uword coilIndex;
  uword shotIndex;
  // Shot loop. Each shot has it's own kspace trajectory
  // for (unsigned int jj = 0; jj < Ns; jj++) {
  std::cout << "Task Number = " << (*taskList)[pWorld->rank()].size()
            << std::endl;
  // Coil loop. Each coil exists for each shot, so we need to work with these.
  for (uword jj = 0; jj < (*taskList)[pWorld->rank()].size(); jj++) {
    taskIndex = (*taskList)[pWorld->rank()].at(jj);
    shotIndex = shotList(taskIndex);
    coilIndex = coilList(taskIndex);
    tempOutData.col(jj) =
        (*AObj[jj]) *
        (d % (SMap.col(coilIndex) % exp(-i * (PMap.col(shotIndex)))));
  }

  //}
  // Now let's do some MPI stuff here.
  if (pWorld->rank() == 0) {
    std::vector<Mat<complex<T1>>> OutDataGather(pWorld->size());
    // Collect all the data into OutDataGather an std::vector collective
    // std::cout << "Rank #: " << pWorld->rank() << " reached foward xform
    // gather" << std::endl;
    bmpi::gather<Mat<complex<T1>>>(*pWorld, tempOutData, OutDataGather, 0);
    // std::cout << "Rank #: " << pWorld->rank() << " passed forward xform
    // gather" << std::endl;
    for (uword jj = 0; jj < pWorld->size(); jj++) {
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
    //std::cout << "Rank #: " << pWorld->rank() << " reached foward xform gather"
    //          << std::endl;
    bmpi::gather<Mat<complex<T1>>>(*pWorld, tempOutData, 0);
  }
  //std::cout << "Rank #: " << pWorld->rank() << " reached foward xform broadcast"
  //          << std::endl;
  bmpi::broadcast(*pWorld, outData, 0);
  //std::cout << "Rank #: " << pWorld->rank() << " passed foward xform broadcast"
  //          << std::endl;

  // equivalent to returning col(output) in MATLAB with IRT
  return vectorise(outData);
}

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
template <typename T1>
Col<complex<T1>> mpipcSENSETimeSeg<T1>::operator/(const Col<complex<T1>> &d) const {
  // uword ShotRank = pWorld->rank();
  Mat<complex<T1>> inData = reshape(d, Nd, Ns * Nc);
  uword taskIndex;
  uword coilIndex;
  uword shotIndex;
  // Col <complex<T1>> tempOutData = zeros < Col < complex<T1> >> (Ni);
  Col<complex<T1>> outData = zeros<Col<complex<T1>>>(Ni);
  // Shot Loop. Each shot has it's own k-space trajectory
  // for (unsigned int jj = 0; jj < Ns; jj++) {
  // Coil Loop - for each shot we have a full set of coil data.
  for (uword jj = 0; jj < (*taskList)[pWorld->rank()].size(); jj++) {
    taskIndex = (*taskList)[pWorld->rank()].at(jj);
    shotIndex = shotList(taskIndex);
    coilIndex = coilList(taskIndex);
    outData += conj(SMap.col(coilIndex) % exp(-i * (PMap.col(shotIndex)))) %
               ((*AObj[jj]) / inData.col(shotIndex + coilIndex * Ns));
    // std::cout << "Processed shot # " << jj << " coil # " << ii <<
    // std::endl;
  }

  if (pWorld->rank() == 0) {
    std::vector<Col<complex<T1>>> OutDataGather;
    // Collect all the data into OutDataGather an std::vector collective
    bmpi::gather<Col<complex<T1>>>(*pWorld, outData, OutDataGather, 0);
    outData.zeros(Ni);
    for (uword jj = 0; jj < pWorld->size(); jj++) {
      // Skip the zeroth rank because that is the outData we started with on
      // this processor.
      // Now we will manually reduce the data by adding across the elements.
      outData += OutDataGather.at(jj);
    }
  } else {
    bmpi::gather<Col<complex<T1>>>(*pWorld, outData, 0);
  }
  // Broadcast will send the data to all machines if the rank==0 and receive
  // the broadcasted data from rank=0 to
  // all other machines in the communicator, overwriting the outData value
  // already there.
  bmpi::broadcast<Col<complex<T1>>>(*pWorld, outData, 0);
  // equivalent to returning col(output) in MATLAB with IRT
  return vectorise(outData);
}

// Explict Instantiation
template class mpipcSENSETimeSeg<float>;
template class mpipcSENSETimeSeg<double>;

#ifdef PowerGridMPI
#pragma message("building reconSolve MPI verison")
template  Col<complex<float>> reconSolve(Col<complex<float>>, mpipcSENSETimeSeg<float>&,
                               QuadPenalty<float>, Col<float>, Col<float>,
                               Col<float>, uword, uword, uword, Col<float>,
                               uword);
#endif
