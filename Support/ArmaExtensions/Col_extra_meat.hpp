// Code based on mlpack serializaton code for Armadillo Mat objects.
// Modified by Alex Cerjanic on 1/18/2016
//mlpack is provided without any warranty of fitness for any purpose.  You
//        can redistribute the library and/or modify it under the terms of the 3-clause
//        BSD license.  The text of the 3-clause BSD license is contained below.
//
//----
//Copyright (c) 2007-2014, mlpack contributors (see COPYRIGHT.txt)
//All rights reserved.
//
//Redistribution and use of mlpack in source and binary forms, with or without
//        modification, are permitted provided that the following conditions are met:
//
//1. Redistributions of source code must retain the above copyright notice, this
//list of conditions and the following disclaimer.
//
//2. Redistributions in binary form must reproduce the above copyright notice,
//this list of conditions and the following disclaimer in the documentation and/or
//other materials provided with the distribution.
//
//3. Neither the name of the copyright holder nor the names of its contributors
//        may be used to endorse or promote products derived from this software without
//specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
//ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
//ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POWERGRID_COL_EXTRA_MEAT_HPP
#define POWERGRID_COL_EXTRA_MEAT_HPP

#pragma message "Including Col_extra_meat.hpp"

template<typename eT>
template<typename Archive>
void Col<eT>::serialize(Archive &ar, const unsigned int /* version */) {
    using boost::serialization::make_nvp;
    using boost::serialization::make_array;

    const uword old_n_elem = Mat<eT>::n_rows;

    // This is accurate from Armadillo 3.6.0 onwards.
    // We can't use BOOST_SERIALIZATION_NVP() because of the access::rw() call.
    ar & make_nvp("n_rows", access::rw(Mat<eT>::n_rows));
    ar & make_nvp("n_cols", access::rw(Mat<eT>::n_cols));
    ar & make_nvp("n_elem", access::rw(Mat<eT>::n_elem));
    ar & make_nvp("vec_state", access::rw(Mat<eT>::vec_state));

    // mem_state will always be 0 on load, so we don't need to save it.
    if (Archive::is_loading::value) {
        // Don't free if local memory is being used.
        if (Mat<eT>::mem_state == 0 && Mat<eT>::mem != NULL && old_n_elem > arma_config::mat_prealloc) {
            memory::release(access::rw(Mat<eT>::mem));
        }

        access::rw(Mat<eT>::mem_state) = 0;

        // We also need to allocate the memory we're using.
        Mat<eT>::init_cold();
    }

    ar & make_array(access::rwp(Mat<eT>::mem), Mat<eT>::n_elem);
}

#endif //POWERGRID_COL_EXTRA_MEAT_HPP
