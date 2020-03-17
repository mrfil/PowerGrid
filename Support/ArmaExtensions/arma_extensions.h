//
// Created by Alex Cerjanic on 1/18/16.
//

#ifndef POWERGRID_ARMA_EXTENSIONS_H
#define POWERGRID_ARMA_EXTENSIONS_H

#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/serialization.hpp>

// The armadillo classes define these macros to include files at the end of the
// class so you can modify them
// without touching the header files, must come prior to <armadillo>

#define ARMA_EXTRA_MAT_PROTO Mat_extra_bones.hpp
#define ARMA_EXTRA_MAT_MEAT Mat_extra_meat.hpp

// The armadillo classes define these macros to include files at the end of the
// class so you can modify them
// without touching the header files, must come prior to <armadillo>
#define ARMA_EXTRA_COL_PROTO Col_extra_bones.hpp
#define ARMA_EXTRA_COL_MEAT Col_extra_meat.hpp

#include <armadillo>

// These headers add a circshift routine to armadillo in the arma namespace
namespace arma {

#include "op_circshift_bones.hpp"
#include "fn_circshift.hpp"
#include "op_circshift_meat.hpp"
//#include "permute.hpp"
}

#endif // POWERGRID_ARMA_EXTENSIONS_H
