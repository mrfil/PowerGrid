//
//  op_circshift_meat.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/1/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_op_circshift_meat_hpp
#define PowerGrid_op_circshift_meat_hpp

template<typename T1>
inline
void
op_circshift::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_circshift>& in)
{
    arma_extra_debug_sigprint();
    
    typedef typename T1::elem_type eT;
    
    const unwrap<T1>   tmp(in.m);
    const Mat<eT>& X = tmp.M;
    
    if(X.is_empty()) { out.copy_size(X); return; }
    
    const uword dim = in.aux_uword_a;
    
    arma_debug_check( (dim > 1), "circshift(): dim must be 0 or 1" );
    
    const uword shift = in.aux_uword_b;
    
    const bool is_alias = (&out == &X);
    
    uword xshift = 0;
    uword yshift = 0;
    uword xdim = X.n_cols;
    uword ydim = X.n_rows;

    if(is_alias == 0) {
        
        out.copy_size(X);

        if(dim == 0) {
            xshift = shift;
            for (uword i = 0; i < xdim; i++) {
                uword ii = positive_modulo(i + xshift, xdim);
                for (uword j = 0; j < ydim; j++) {
                    uword jj = positive_modulo(j + yshift, ydim);
                    out(jj,ii) = X(j,i);
                }
            }
        } else {
            yshift = shift;
            for (uword i = 0; i < xdim; i++) {
                uword ii = positive_modulo(i + xshift, xdim);
                for (uword j = 0; j < ydim; j++) {
                    uword jj = positive_modulo(j + yshift, ydim);
                    out(jj,ii) = X(j,i);
                }
            }
        }
    } else { //X is an alias of out (same memory address)
        Mat<eT> temp;
        temp.copy_size(X);
        temp = X;
        
        if(dim == 0) {
            xshift = shift;
            for (uword i = 0; i < xdim; i++) {
                uword ii = positive_modulo(i + xshift, xdim);
                for (uword j = 0; j < ydim; j++) {
                    uword jj = positive_modulo(j + yshift, ydim);
                    out(jj,ii) = temp(j,i);
                }
            }
        } else {
            yshift = shift;
            for (uword i = 0; i < xdim; i++) {
                uword ii = positive_modulo(i + xshift, xdim);
                for (uword j = 0; j < ydim; j++) {
                    uword jj = positive_modulo(j + yshift, ydim);
                    out(jj,ii) = temp(j,i);
                }
            }
        }
    }
    
}

template<typename T1>
inline
T1
op_circshift::positive_modulo(T1 i, T1 n)
{
    return (i % n + n) % n;
}

#endif
