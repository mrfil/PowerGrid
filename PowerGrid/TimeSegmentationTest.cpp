/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [TimeSegmentation.cpp]

    Synopsis    [Wrappers to the cuFFT library supporting single and double
        precision for GPU accelerated FFTs]

    Description []

    Revision    [0.1.0; Giang-Chau Ngo, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/
#include "TimeSegmentation.h"
#include <chrono>  // for high_resolution_clock
// This using field correction by time segmentation
// The data is corrected to time 0 with reference to the time vector passed
//
//

// We are using two template types at the moment. One for the type of data to be
// processed (ie Col<cx_double>) and one for the type of G object (ie
// Gfft<Col<cx_double>>
// T1 is the data type for complex, T2 is the data type for real data
template <typename T1, typename Tobj>
TimeSegmentation<T1, Tobj>::TimeSegmentation(Tobj &G, Col<T1> map_in,
                                             Col<T1> timeVec_in, uword a,
                                             uword b, uword c, uword interptype,
                                             uword shots) {

  cout << "Entering Class constructor" << endl;
  n1 = a;            // Data size
  n2 = b;            // Image size

  // L == 0 or 1 is ambiguous. In matlab we would often specific L = 0 to mean no time segmentation.
  if ((L == 0) || (L == 1)) {
    L = 1;
  } else {  
    L = c;  
  }           // number of time segments
  type = interptype; // type of time segmentation performed
  Nshots = shots;    // number of shots
  obj = &G;
  fieldMap = map_in;
  outData.set_size(n1,L);
  outImg.set_size(n2,L);
  tempD.set_size(n1,L);
  tempAA.set_size(n1,L);
  RowOnes.ones(L,1);
  
  cout << "N1 = " << n1 << endl;
  cout << "N2 = " << n2 << endl;
  cout << "L = " << L << endl;

  AA.set_size(n1, L); // time segments weights
  timeVec = timeVec_in;
  T_min = timeVec.min();
  T1 rangt = timeVec.max() - T_min;
  tau = (rangt + datum::eps) / (L - 1); // it was L-1 before
  timeVec = timeVec - T_min;

  uword NOneShot = n1 / Nshots;
  if ((L == 1)) {
    tau = 0;
    AA.ones();
    Wo.ones(n2, 1);
    WoH.ones(n2, 1);
    
  } else {
    Mat<complex<T1>> tempAA(NOneShot, L);
    if (type == 1) { // Hanning interpolator
      cout << "Hanning interpolation" << endl;
      for (unsigned int ii = 0; ii < L; ii++) {
        for (unsigned int jj = 0; jj < NOneShot; jj++) {
          if ((abs(timeVec(jj) - ((ii)*tau))) <= tau) {
            tempAA(jj, ii) =
                0.5 +
                0.5 * std::cos((datum::pi) * (timeVec(jj) - ((ii)*tau)) / tau);
          } else {
            tempAA(jj, ii) = 0.0;
          }
        }
      }
      cout << "Precalculating interpolators" << endl;

      AA = repmat(tempAA, Nshots, 1);
      Wo.set_size(n2, L);
      WoH.set_size(n2, L);
      #pragma omp parallel for shared(Wo, WoH)
      for (unsigned int ii = 0; ii < L; ii++) {
        Wo.col(ii) =
            exp(-i * (this->fieldMap) * ((ii) * this->tau + this->T_min));
        WoH.col(ii) =
            exp(i * (this->fieldMap) * ((ii) * this->tau + this->T_min));
      }
      //Wo.save("Wo.dat",raw_ascii);
      //WoH.save("WoH.dat",raw_ascii);
    } else if (type == 2) { // Min-max interpolator: Exact LS interpolator

      cout << "Min Max time segmentation" << endl;

      Mat<complex<T1>> Ltp;
      Ltp.ones(1, L);
      Col<complex<T1>> ggtp;
      ggtp.ones(n2, 1);
      Mat<complex<T1>> gg;
      gg = exp(i * fieldMap * tau) * Ltp;
      Mat<complex<T1>> iGTGGT;
      iGTGGT.set_size(L + 1, n2);
      Mat<complex<T1>> gl;
      gl.zeros(n2, L);

      for (unsigned int ii = 0; ii < L; ii++) {
        for (unsigned int jj = 0; jj < n2; jj++) {
          gl(jj, ii) = pow(gg(jj, ii), (T1)(ii + 1));
        }
      }

      Mat<complex<T1>> G;
      G.set_size(n2, L);

      for (unsigned int jj = 0; jj < L; jj++) {
        if (jj == 0) {
          G.col(jj) = ggtp;
        } else {
          G.col(jj) = gl.col(jj - 1);
        }
      }

      Col<complex<T1>> glsum;
      Mat<complex<T1>> GTG;
      GTG.zeros(L, L);
      GTG.diag(0) += n2;
      glsum = sum(gl.t(), 1);
      Mat<complex<T1>> GTGtp(L, L);
      for (unsigned int ii = 0; ii < (L - 1); ii++) {
        GTGtp.zeros();
        GTGtp.diag(-(T1)(ii + 1)) += glsum(ii);
        GTGtp.diag((T1)(ii + 1)) += std::conj(glsum(ii));
        GTG = GTG + GTGtp;
      }

      T1 rcn = 1 / cond(GTG);
      if (rcn > 10 * 2e-16) { // condition number of GTG
        iGTGGT = inv(GTG) * G.t();

      } else {
        iGTGGT = pinv(GTG) * G.t(); // pseudo inverse
      }

      Mat<complex<T1>> iGTGGTtp;
      Mat<complex<T1>> ftp;
      Col<complex<T1>> res, temp;

      for (unsigned int ii = 0; ii < NOneShot; ii++) {
        ftp = exp(i * fieldMap * timeVec(ii));
        res = iGTGGT * ftp;
        tempAA.row(ii) = res.t();
      }
      AA = repmat(tempAA, Nshots, 1);
    }
  }
  // savemat("aamat.mat", "AA", vectorise(AA));
  cout << "Exiting class constructor." << endl;
}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1, typename Tobj>
inline Col<complex<T1>> TimeSegmentation<T1, Tobj>::
operator*(const Col<complex<T1>> &d) const {
  RANGE()
  auto start = std::chrono::high_resolution_clock::now();

  Tobj *G = this->obj;
  // output is the size of the kspace data
  //Col<complex<T1>> outData = zeros<Col<complex<T1>>>(this->n1);
  //outData.zeros();
  
  // cout << "OutData size = " << this->n1 << endl;
  //Col<complex<T1>> Wo;
  
  //Col<complex<T1>> temp;
  //tempD = this->Wo.each_col() % d;
  //uvec dataMaskTrimmed;
  // loop through time segments
  
  #pragma omp parallel for shared(tempD)
  for (unsigned int ii = 0; ii < this->L; ii++){
    tempD.unsafe_col(ii) = Wo.unsafe_col(ii) % d;
  }

  for (unsigned int ii = 0; ii < this->L; ii++) {
    

    // apply a phase to each time segment
    //Wo = exp(-i * (this->fieldMap) * ((ii) * this->tau + this->T_min));

    // perform multiplication by the object and sum up the time segments
    //temp = (this->Wo.col(ii)) % d;

    //outData.unsafe_col(ii) =  (*G * tempD.unsafe_col(ii));
    outData.col(ii) =  (*G * tempD.col(ii));
    //outData +=  AA.unsafe_col(ii) % (*G * tempD.unsafe_col(ii));
 
  }

  #pragma omp parallel for shared(outData)
  for (unsigned int ii = 0; ii < this->L; ii++){
    outData.col(ii) %= AA.col(ii);
  }
  //outData %= AA;


auto finish = std::chrono::high_resolution_clock::now();

  std::chrono::duration<float> elapsed = finish - start;
  std::cout << "Time Segmentation For Elapsed time: " << elapsed.count() << " s\n";
  
  return sum(outData,1);
}
template <typename T1, typename Tobj>
inline Col<complex<T1>> TimeSegmentation<T1, Tobj>::
operator/(const Col<complex<T1>> &d) const {
  RANGE()
  auto start = std::chrono::high_resolution_clock::now();

  Tobj *G = this->obj;
  //outImg.zeros();
  /*
  #pragma omp parallel for shared(tempAA)
  for (unsigned int ii = 0; ii < this->L; ii++){
    tempAA.col(ii) = AA.col(ii) % d;
  }
  */
  // loop through the time segments
  for (unsigned int ii = 0; ii < this->L; ii++) {

    // perform adjoint operation by the object and sum up the time segments
    outImg.col(ii) = WoH.col(ii) % ((*G) / (AA.col(ii) % d));
  }
  /*
  #pragma omp parallel for shared(outImg)
  for (unsigned int ii = 0; ii < this->L; ii++) {
    outImg.col(ii) %= WoH.col(ii);
  }
  */
  auto finish = std::chrono::high_resolution_clock::now();

  std::chrono::duration<float> elapsed = finish - start;
  std::cout << "Time Segmentation Adj Elapsed time: " << elapsed.count() << " s\n";
  
  return sum(outImg,1);
}

// Explicit Instantiations
template class TimeSegmentation<float, Gnufft<float>>;
template class TimeSegmentation<double, Gnufft<double>>;
