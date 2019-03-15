/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [processNIFTI.hpp]

    Synopsis    [Helper functions for working with NIFTI files.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [2018/02/14]

 *****************************************************************************/

#ifndef POWERGRID_PROCESSNIFTI_HPP
#define POWERGRID_PROCESSNIFTI_HPP

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

#include "PowerGrid.h"

using namespace arma;

#include <iostream>
#include <string>
#include <cstring>
#include "../Support/nifti1.h"

template<typename T1>
void writeNiftiRealImage(std::string filename, Col<T1> imageData, uword Nx, uword Ny, uword Nz) {

  nifti_1_header hdr;
  nifti1_extender pad = {0,0,0};

  // Initialize header
  hdr.sizeof_hdr = MIN_HEADER_SIZE;
  hdr.dim[0] = 3;
  hdr.dim[1] = Nx;
  hdr.dim[2] = Ny;
  hdr.dim[3] = Nz;
  hdr.dim[4] = 1;
  hdr.datatype = NIFTI_TYPE_FLOAT32;
  hdr.bitpix = 32;
  hdr.pixdim[0] = 1.0;
  hdr.pixdim[1] = 1.0;
  hdr.pixdim[2] = 1.0;
  hdr.pixdim[3] = 1.0; //STUPID!!

  hdr.quatern_b = 0;
  hdr.quatern_c = 0;
  hdr.quatern_d = 0;
  hdr.qoffset_x = 0;
  hdr.qoffset_y = 0;
  hdr.qoffset_z = 0;

  for (int ii = 0; ii < 4; ii++) {
    hdr.srow_x[ii] = 0.0f;
    hdr.srow_y[ii] = 0.0f;
    hdr.srow_z[ii] = 0.0f;
  }



  hdr.vox_offset = (float) NII_HEADER_SIZE;
  hdr.scl_slope = 100.0;
  hdr.xyzt_units = NIFTI_UNITS_MM;

  // set magic number for file format
  std::strncpy(hdr.magic, "n+1\0",4);

  //Writing the file based on nifti1_read_write.c from NIH.
  // Open file
  std::ofstream ofsNii(filename);

  // Write header to disk
  ofsNii.write((const char *) &hdr, MIN_HEADER_SIZE);

  // Need to write a pad to disk
  ofsNii.write((const char *) &pad,4);
  float temp = 0;
  // Write image data to disk
  for(uword ii = 0; ii < Nx*Ny*Nz; ii++) {
    temp = static_cast<float>(imageData(ii));
    ofsNii.write(reinterpret_cast<char *>(&temp), sizeof(float));
  }

  ofsNii.close();
}

  // Start with Magnitude

  template<typename T1>
  void writeNiftiMagPhsImage(std::string filename, Col<std::complex<T1>> imageData, uword Nx, uword Ny, uword Nz) {
    // Separate data into magnitude and phase to be written.
    Col<T1> magImage = abs(imageData);
    Col<T1> phsImage = arg(imageData);

    std::string filenameMag = filename + "_mag" + ".nii";
    std::string filenamePhs = filename + "_phs" + ".nii";

    //Write out separate magnitude and phase images
    writeNiftiRealImage<T1>(filenameMag, magImage,Nx,Ny,Nz);
    writeNiftiRealImage<T1>(filenamePhs, phsImage,Nx,Ny,Nz);

  }




#endif //POWERGRID_PROCESSNIFTI_HPP
