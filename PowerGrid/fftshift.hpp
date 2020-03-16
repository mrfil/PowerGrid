//
//  fftshift.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/2/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__fftshift__hpp
#define __PowerGrid__fftshift__hpp
#include "armadillo"


using namespace arma;
using namespace std;

template<typename T1>
arma_inline
T1 fftshift(T1& X, uword dim);

template<typename T1>
arma_inline
T1 fftshift(T1 X);

template<typename T1>
arma_inline
T1 fftshift(
            T1& X,
            uword dim
            ) {
    uword size = 0;
    if (dim == 0) {
        size = X.n_rows;
    } else if(dim == 1) {
        size = X.n_cols;
    } else if(dim > 1) {
        arma_extra_debug_sigprint();
        arma_debug_print( "fftshift(): trying to shift dimension greater than 2");
    }

    T1 out = circshift(X,dim,std::floor(size/2));
    return out;
}

template<typename T1>
arma_inline
T1 fftshift(
            T1 X
            ) {

    T1 out = circshift(circshift(X,0,std::floor(X.n_rows/2)),1,std::floor(X.n_cols/2));
    return out;
}



#endif /* defined(__PowerGrid__fftshift__hpp) */
