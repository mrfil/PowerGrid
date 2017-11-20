//
//  fftshift.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/12/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_fn_circshift_hpp
#define PowerGrid_fn_circshift_hpp




using namespace std;

//  Create an armadillo compatible circshift function similar to MATLAB's

template<typename T1>
arma_inline
const Op<T1, op_circshift>
circshift
(
  const T1& X,
  const uword n_dim,
  const uword shift
)
{
    arma_extra_debug_sigprint();
    arma_debug_check ( (n_dim > 3), "circshfit(): trying to shift dimension greater than 3");
    
    return Op<T1, op_circshift>(X, n_dim, shift);
    
}




#endif
