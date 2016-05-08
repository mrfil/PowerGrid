/*
 This file is distributed under a BSD 3-Clause license.
Copyright (c) 2014, Ceemple Software Ltd.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
 following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
 disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
 disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef CEEMPLEARMADILLO_H
#define CEEMPLEARMADILLO_H

#include "CeempleComplex.h"
#include "armadillo"

template <typename T>
inline arma::Mat<T> join_cols(const arma::Mat<T> &A, const double b) {
  arma::Mat<T> B;
  B = (T)b;

  return join_cols(A, B);
}

template <typename T>
inline arma::Mat<T> join_cols(const double a, const arma::Mat<T> &B) {
  arma::Mat<T> A;
  A = (T)a;

  return join_cols(A, B);
}

template <typename T>
inline arma::Mat<T> join_rows(const arma::Mat<T> &A, const double b) {
  arma::Mat<T> B;
  B = (T)b;

  return join_rows(A, B);
}

template <typename T>
inline arma::Mat<T> join_rows(const double a, const arma::Mat<T> &B) {
  arma::Mat<T> A;
  A = (T)a;

  return join_rows(A, B);
}

template <typename T> inline arma::uword size(const arma::Mat<T> &A, double d) {
  if (d == 1)
    return A.n_rows;
  else if (d == 2)
    return A.n_cols;
  else
    return 0;
}

#endif // CEEMPLEARMADILLO_H
