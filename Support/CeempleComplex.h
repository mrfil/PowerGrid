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

#ifndef CEEMPLECOMPLEX_H
#define CEEMPLECOMPLEX_H

#include <complex>
// complex, int
template <class T> std::complex<T> operator+(const std::complex<T> &x, int y) {
  return x + double(y);
}
template <class T> std::complex<T> operator+(int x, const std::complex<T> &y) {
  return double(x) + y;
}
template <class T> std::complex<T> operator-(const std::complex<T> &x, int y) {
  return x - double(y);
}
template <class T> std::complex<T> operator-(int x, const std::complex<T> &y) {
  return double(x) - y;
}
template <class T> std::complex<T> operator*(const std::complex<T> &x, int y) {
  return x * double(y);
}
template <class T> std::complex<T> operator*(int x, const std::complex<T> &y) {
  return double(x) * y;
}
template <class T> std::complex<T> operator/(const std::complex<T> &x, int y) {
  return x / double(y);
}
template <class T> std::complex<T> operator/(int x, const std::complex<T> &y) {
  return double(x) / y;
}

template <class T> std::complex<T> operator+(const std::complex<T> &x, float y) {
  return x + double(y);
}
template <class T> std::complex<T> operator+(float x, const std::complex<T> &y) {
  return double(x) + y;
}
template <class T> std::complex<T> operator-(const std::complex<T> &x, float y) {
  return x - double(y);
}
template <class T> std::complex<T> operator-(float x, const std::complex<T> &y) {
  return double(x) - y;
}
template <class T> std::complex<T> operator*(const std::complex<T> &x, float y) {
  return x * double(y);
}
template <class T> std::complex<T> operator*(float x, const std::complex<T> &y) {
  return double(x) * y;
}
template <class T> std::complex<T> operator/(const std::complex<T> &x, float y) {
  return x / double(y);
}
template <class T> std::complex<T> operator/(float x, const std::complex<T> &y) {
  return double(x) / y;
}

#endif // CEEMPLECOMPLEX_H
