/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [TimeSegmentation.h]

    Synopsis    [Wrappers to the cuFFT library supporting single and double
                                precision for GPU accelerated FFTs]

    Description []

    Revision    [0.1.0; Giang-Chau Ngo, BIOE UIUC
                 0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/2/2016]

 *****************************************************************************/
// We are using two template types at the moment. One for the type of data to be
// processed (ie Col<cx_double>) and one for the type of G object (ie
// Gfft<Col<cx_double>>
// T1 is the data type for complex, T2 is the data type for real data

#ifndef PowerGrid_TimeSegmentation_h
#define PowerGrid_TimeSegmentation_h

#include "Gdft.h"
#include "Gnufft.h"
#include "PGIncludes.h"
using namespace std;

template <typename T1, typename Tobj> class TimeSegmentation {
  typedef complex<T1> CxT1;

public:
  TimeSegmentation();

  // Class variables go here
  uword n1;     // Data size
  uword n2;     // Image size
  int L;        // number of time segments
  uword type;   // type of time segmentation
  uword Nshots; // Number of shots, used to reduce complexity of calculating
                // interpolator
  T1 tau;       // time segment length
  T1 T_min;     // minimum time in the time vector (i.e. TE for spiral out)
  Tobj *obj;
  Col<T1> fieldMap; // Field map (in radians per second)
  Col<T1> timeVec;  // timing vector of when each data point was collected
                    // relative to the echo time (in seconds)
  Mat<CxT1> AA;     // interpolator coefficients for the different time segments
  CxT1 i = CxT1(0., 1.);
  Mat<CxT1> Wo;
  Mat<CxT1> WoH;
  Col<T1> RowOnes;
  mutable Mat<complex<T1>> outData;
  mutable Mat<complex<T1>> outImg;
  mutable Mat<complex<T1>> tempD;
  mutable Mat<complex<T1>> tempAD;

  // Class constructor
  TimeSegmentation(Tobj &G, Col<T1> map_in, Col<T1> timeVec_in, uword a,
                   uword b, uword c, uword interptype = 1, uword shots = 1);

  // Overloaded operators go here
  Col<CxT1> operator*(const Col<CxT1> &d) const;
  Col<CxT1> operator/(const Col<CxT1> &d) const;

  protected:

  // Handy utility function for fftw via armadillo
  Col<CxT1> calcFFT1D(const Col<CxT1> &d, uword KK) const;
};

// Now we insert the explicit instantiations we need
extern template class TimeSegmentation<float, Gnufft<float>>;
extern template class TimeSegmentation<float, Gdft<float>>;
extern template class TimeSegmentation<double, Gnufft<double>>;
extern template class TimeSegmentation<double, Gdft<double>>;
#endif // PowerGrid_TimeSegmentation_h
