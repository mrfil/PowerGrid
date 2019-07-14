/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [solve_pwls_pcg.hpp]

    Synopsis    [Object implementing a preconditioned weighted least squares
                    solver using preconditioned conjugate gradient.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef POWERGRID_SOLVE_PWLS_PCG_HPP_
#define POWERGRID_SOLVE_PWLS_PCG_HPP_

#include <cstdlib>
#include <chrono>
using namespace arma;

template <typename T1>
inline complex<T1> dot_double(const Col<complex<T1>>& A,
                              const Col<complex<T1>>& B) {
  complex<T1> sumReturn = accu(A % B);
  return sumReturn;
}

template <typename T1>
inline T1 norm_grad(const Col<complex<T1>> &g, const Col<complex<T1>> &yi,
                    const Col<T1> &W) {
  T1 normGrad = conv_to<T1>::from(norm(g) / real(trans(yi) * (W % yi)));
  return normGrad;
}

template <typename T1, typename Tobj, typename Robj>
Col<complex<T1>> solve_pwls_pcg(const Col<complex<T1>> &xInitial, Tobj const &A,
                                Col<T1> const &W, Col<complex<T1>> const &yi,
                                Robj const &R, uword niter) {
  typedef complex<T1> CxT1;
  // Initialize projection
  cout << "Entering solve_pwls_pcg" << endl;
  Col<CxT1> Ax = A * xInitial;
  if (Ax.has_nan())
    cout << "Warning: Ax has NaN in solve_pwls_pcg" << endl;
  // cout << "Ax length = " << Ax.n_rows << endl;
  Col<CxT1> x = xInitial;
  CxT1 oldinprod = 0;
  CxT1 gamma = 0.0;
  Col<CxT1> ddir;
  Col<CxT1> Adir;
  CxT1 dAWAd;
  CxT1 dAWr;
  CxT1 pdenom;
  CxT1 denom;

  Col<CxT1> ngrad;
  Col<CxT1> pgrad;
  CxT1 pdot;
  Col<CxT1> WAdir;
  //CxT1 proj;
  Col<CxT1> stepIntermediate;
  CxT1 step;
  //CxT1 rdenom;
  CxT1 newinprod;

  cout << "Entering solve_pwls_pcg iteration loop" << endl;
  for (unsigned int ii = 0; ii < niter; ii++) {
    // Compute negative gradient

    ngrad = A / (W % (yi - Ax));
    if(ngrad.has_nan())
      cout << "Warning: ngrad has NaN in solve_pwls_pcg" << endl;


    if (norm_grad<T1>(ngrad, yi, W) < 1e-10) {
      cout << "Terminating early due to zero gradient." << endl;
      return x;
    }
    ngrad -= R.Gradient(x);

    // Direction
    newinprod = real(cdot(ngrad,ngrad));
    if (ii == 0) {
      ddir = ngrad;

    } else {
      if (std::abs(oldinprod) < 1e-10) {
        gamma = 0.0;
      } else {
        gamma = newinprod / oldinprod;
      }

      ddir = ngrad + gamma * ddir;
    }

    Col<CxT1> oldgrad = ngrad;
    oldinprod = newinprod;

    // Check if descent direction
    if (real(cdot(ddir, ngrad)) < 0) {
      cout << " Warning descent direction not negative" << endl;
      return x;
    }

    // Step size in search direction
    Adir = A * ddir;
    if (Adir.has_nan())
      cout << "Warning: NaN found in Adir in solve_pwls_pcg" << endl;

    WAdir = W % Adir;
    dAWAd = as_scalar(real(cdot(Adir, WAdir)));
    dAWr  = as_scalar(real(Adir.t() * (W % (yi - Ax))));

    step = 0.0;

    for (unsigned int j = 0; j < 3; j++) {
      // pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir); Original MATLAB code
      // from pwls_pcg1.m
      pdenom = R.Denom(ddir, x + step * ddir);
      denom = dAWAd + pdenom;
      if( denom != denom)
        cout << "Warning: denom has NaN in solve_pwls_pcg" << endl;

      if (std::abs(denom) < 1e-20 || std::abs(denom) > 1e25) {
        if (norm(ngrad, 2) == 0) {
          cout << " Found exact solution" << endl;
          return x;
        } else {
          cout << "inf denom" << endl;
          return x;
        }
      }

      pgrad = R.Gradient(x + step * ddir);
      pdot = real(cdot(ddir, pgrad));

      stepIntermediate = (-dAWr + step * dAWAd + pdot) / denom;
      step -= as_scalar(stepIntermediate);
    }
    //cout << "denom = " << std::abs(denom) << endl;
    // Check downhill direction
    if (as_scalar(real(step)) < 0) {
      cout << "Warning downhill?" << endl;
    }

    // Update
    Ax += step * Adir;
    x += (step * ddir);
    cout << "Iteration Error Norm = " << norm(yi - Ax, 2) << endl;
    cout << "Iteration Complete = " << ii << endl;

  }
  return x;
}

#endif /* POWERGRID_SOLVE_PWLS_PCG_HPP_ */
