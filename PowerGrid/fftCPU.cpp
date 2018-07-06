/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [fftCPU.cpp]

    Synopsis    [Wrappers to the FFTW library supporting single and double
                                precision]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]
                [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/13/2016]

 *****************************************************************************/

#include "fftCPU.h"

// Like Armadillo, we're using SFINAE here to choose between float and double.
// (Maybe FP16 some day in the future)
// We need enable_if to choose which version to run based on the type of the
// template parameter.

template <typename T1,
          typename std::enable_if<std::is_same<T1, float>::value, uword>::type>
void fft1dCPU(T1 *d_data, uword nx) {

  fftwf_plan plan;
  plan =
      fftwf_plan_dft_1d(nx, (fftwf_complex *)d_data,
                        (fftwf_complex *)d_data, FFTW_FORWARD, FFTW_ESTIMATE);

  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, float>::value, uword>::type>
void ifft2dCPU(T1 *d_data, uword nx, uword ny) {
  cout << "Running backward xform 2d" << endl;
  fftwf_plan plan;
  plan =
      fftwf_plan_dft_2d(ny, nx, (fftwf_complex *)d_data,
                        (fftwf_complex *)d_data, FFTW_BACKWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place.

  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, float>::value, uword>::type>
void fft2dCPU(T1 *d_data, uword nx, uword ny) {
  cout << "Running forward xform 2d" << endl;
  fftwf_plan plan;
  plan =
      fftwf_plan_dft_2d(ny, nx, (fftwf_complex *)d_data,
                        (fftwf_complex *)d_data, FFTW_FORWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place.
  // fftwf_pruword_plan(plan);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, float>::value, uword>::type>
void ifft3dCPU(T1 *d_data, uword nx, uword ny, uword nz) {
  cout << "Running backward xform 3d" << endl;
  fftwf_plan plan;
  plan =
      fftwf_plan_dft_3d(nz, ny, nx, (fftwf_complex *)d_data,
                        (fftwf_complex *)d_data, FFTW_BACKWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place.
  // fftwf_pruword_plan(plan);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, float>::value, uword>::type>
void fft3dCPU(T1 *d_data, uword nx, uword ny, uword nz) {
  cout << "Running forward xform 3d" << endl;
  fftwf_plan plan;
  plan =
      fftwf_plan_dft_3d(nz, ny, nx, (fftwf_complex *)d_data,
                        (fftwf_complex *)d_data, FFTW_FORWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place.
  // fftwf_pruword_plan(plan);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, double>::value, uword>::type>
void fft1dCPU(T1 *d_data, uword nx) {
  cout << "Running forward xform 2d" << endl;

  fftw_plan plan;
  plan = fftw_plan_dft_1d(nx, (fftw_complex *)d_data,
                          (fftw_complex *)d_data, FFTW_FORWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place.
  fftw_execute(plan);
  fftw_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, double>::value, uword>::type>
void ifft2dCPU(T1 *d_data, uword nx, uword ny) {
  cout << "Running backward xform 2d" << endl;

  fftw_plan plan;
  plan = fftw_plan_dft_2d(ny, nx, (fftw_complex *)d_data,
                          (fftw_complex *)d_data, FFTW_BACKWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place
  fftw_execute(plan);
  fftw_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, double>::value, uword>::type>
void fft2dCPU(T1 *d_data, uword nx, uword ny) {
  cout << "Running forward xform 2d" << endl;

  fftw_plan plan;
  plan = fftw_plan_dft_2d(ny, nx, (fftw_complex *)d_data,
                          (fftw_complex *)d_data, FFTW_FORWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place.
  fftw_execute(plan);
  fftw_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, double>::value, uword>::type>
void ifft3dCPU(T1 *d_data, uword nx, uword ny, uword nz) {
  cout << "Running backward xform 3d" << endl;

  fftw_plan plan;
  plan = fftw_plan_dft_3d(nz, ny, nx, (fftw_complex *)d_data,
                          (fftw_complex *)d_data, FFTW_BACKWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place.
  fftw_execute(plan);
  fftw_destroy_plan(plan);
}

template <typename T1,
          typename std::enable_if<std::is_same<T1, double>::value, uword>::type>
void fft3dCPU(T1 *d_data, uword nx, uword ny, uword nz) {
  cout << "Running forward xform 3d" << endl;

  fftw_plan plan;
  plan = fftw_plan_dft_3d(nz, ny, nx, (fftw_complex *)d_data,
                          (fftw_complex *)d_data, FFTW_FORWARD, FFTW_ESTIMATE);

  // Inverse transform 'gridData_d' in place.
  fftw_execute(plan);
  fftw_destroy_plan(plan);
}

// Explicit Instantiations
template void fft1dCPU<float>(float *, uword);
template void fft1dCPU<double>(double *, uword);

template void ifft2dCPU<float>(float *, uword, uword);
template void ifft2dCPU<double>(double *, uword, uword);

template void fft2dCPU<float>(float *, uword, uword);
template void fft2dCPU<double>(double *, uword, uword);

template void ifft3dCPU<float>(float *, uword, uword, uword);
template void ifft3dCPU<double>(double *, uword, uword, uword);

template void fft3dCPU<float>(float *, uword, uword, uword);
template void fft3dCPU<double>(double *, uword, uword, uword);
