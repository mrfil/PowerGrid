#include "catch.hpp"

#include "../PGIncludes.h"
#include "../pgCol.hpp"
#include <cmath>
#include <iostream>
#include "accel.h"

TEST_CASE("pgCol<float>: operators", "[pgCol<float>]") {

    // Setup Prerequsites for test
    arma::uword lengthA = 100000;
    arma::uword lengthB = 256 * 256;
    
    pgCol<float> pgA(lengthA);
    pgCol<float> pgB(lengthB);
    pgCol<float> pgC;
    pgCol<float> pgD;
    
    pgC.set_size(lengthA);
    pgD.set_size(lengthB);

    pgCol<float> pgCZeros(lengthA);
    pgCol<float> pgDZeros(lengthB);
    
    SECTION( "Test constructor and n_elems " ) {

        REQUIRE(pgA.n_elem == lengthA);
        REQUIRE(pgB.n_elem == lengthB);
    } 



    SECTION( "Test set_size()" ) {

        REQUIRE(pgC.n_elem == lengthA);
        REQUIRE(pgD.n_elem == lengthB);
    }

    SECTION( "Test zeros and sum" ) {
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum(pgC) == 0);
        REQUIRE(sum(pgD) == 0);
    }

    SECTION( "Test ones and sum()" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(pgC) == (float)lengthA);
        REQUIRE(sum(pgD) == (float)lengthB);
    }

    SECTION( "Test Operator+" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum( pgC + pgC ) == 2 * lengthA);
        REQUIRE(sum( pgD + pgD ) == 2 * lengthB);
        
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum( pgC + pgC ) == 0);
        REQUIRE(sum( pgD + pgD ) == 0);
    }

    SECTION( "Test Operator-" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum( pgC - pgC ) == 0);
        REQUIRE(sum( pgD - pgD ) == 0);
    }

    SECTION( "Test Operator%" ) {
        pgC.ones();
        pgD.ones();

        pgCZeros.zeros();
        pgDZeros.zeros();
        REQUIRE(sum( pgC % pgC ) == lengthA);
        REQUIRE(sum( pgD % pgD ) == lengthB);
        REQUIRE(sum( pgC % pgCZeros) == 0);
        REQUIRE(sum( pgD % pgDZeros) == 0);
        
    }

}
