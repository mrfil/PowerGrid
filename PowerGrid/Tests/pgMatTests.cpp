#include "catch.hpp"

#ifdef _OPENACC
#include "accel.h"
#endif

#include "../PGIncludes.h"
#include "../pgMat.hpp"
#include <cmath>
#include <iostream>


TEST_CASE("pgMat<float>: operators", "[pgMat<float>") {

    // Setup Prerequsites for test
    arma::uword lengthA = 256;
    arma::uword lengthB = 1000;
    
    pgMat<float> pgA(lengthA,lengthA);
    pgMat<float> pgB(lengthB,lengthB);
    pgMat<float> pgC;
    pgMat<float> pgD;
    
    pgC.set_size(lengthA,lengthA);
    pgD.set_size(lengthB,lengthB);

    pgMat<float> pgCZeros(lengthA,lengthA);
    pgMat<float> pgDZeros(lengthB,lengthB);
    
    SECTION( "Test constructor and n_elems " ) {

        REQUIRE(pgA.n_elem == lengthA * lengthA);
        REQUIRE(pgB.n_elem == lengthB * lengthB);
    } 

    SECTION( "Test constructor and n_rows " ) {

        REQUIRE(pgA.n_rows == lengthA);
        REQUIRE(pgB.n_rows == lengthB);
    } 

    SECTION( "Test constructor and n_cols " ) {

        REQUIRE(pgA.n_cols == lengthA);
        REQUIRE(pgB.n_cols == lengthB);
    } 

    SECTION( "Test set_size() and n_elem" ) {

        REQUIRE(pgC.n_elem == lengthA * lengthA);
        REQUIRE(pgD.n_elem == lengthB * lengthB);
    }

    SECTION( "Test set_size() and n_rows " ) {

        REQUIRE(pgC.n_rows == lengthA);
        REQUIRE(pgD.n_rows == lengthB);
    } 

    SECTION( "Test set_size() and n_cols " ) {

        REQUIRE(pgC.n_cols == lengthA);
        REQUIRE(pgD.n_cols == lengthB);
    } 


    SECTION( "Test zeros and sum" ) {
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum(sum(pgC)) == 0);
        REQUIRE(sum(sum(pgD)) == 0);
    }

    SECTION( "Test ones and sum()" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(sum(pgC)) == (float)(lengthA * lengthA));
        REQUIRE(sum(sum(pgD)) == (float)(lengthB * lengthB));
    }

    SECTION( "Test Operator+" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(sum( pgC + pgC )) == 2 * lengthA * lengthA);
        REQUIRE(sum(sum( pgD + pgD )) == 2 * lengthB * lengthB);
        
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum(sum( pgC + pgC )) == 0);
        REQUIRE(sum(sum( pgD + pgD )) == 0);
    }

    SECTION( "Test Operator-" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(sum( pgC - pgC )) == 0);
        REQUIRE(sum(sum( pgD - pgD )) == 0);
    }

    SECTION( "Test Operator%" ) {
        pgC.ones();
        pgD.ones();

        pgCZeros.zeros();
        pgDZeros.zeros();
        REQUIRE(sum(sum( pgC % pgC )) == lengthA * lengthA);
        REQUIRE(sum(sum( pgD % pgD )) == lengthB * lengthB);
        REQUIRE(sum(sum( pgC % pgCZeros)) == 0);
        REQUIRE(sum(sum( pgD % pgDZeros)) == 0);
        
    }

}
