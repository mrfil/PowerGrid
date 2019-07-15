#include "catch.hpp"

#include "../pgComplex.hpp"
#include <complex>
#include <cmath>
#include <iostream>


TEST_CASE("pgComplex<float>: operators", "[pgComplex]") {

    // Setup Prerequsites for test
    float re_ = 1.2f;
    float im_ = -3.4f;
    
    pgComplex<float> pgTest(re_,im_); // 1.2f - j3.4
    pgComplex<float> pgZero = {};        // 0 + j0
    pgComplex<float> pgJ(0.0f, 1.0f); // 0 + j1

    SECTION( "Test constructor .real and .imag()" ) {
        REQUIRE( pgZero.real() == 0.0f);
        REQUIRE( pgZero.imag() == 0.0f);
        REQUIRE( pgTest.real() == re_);
        REQUIRE( pgTest.imag() == im_);
    }

    SECTION( "Test operator==" ) {
        REQUIRE( pgTest == pgComplex<float>(re_, im_));
    }

    SECTION( "Test operator+" ) {
        REQUIRE( pgZero + pgTest == pgTest);
        REQUIRE( pgTest + pgTest == pgComplex<float>(2 * re_, 2 * im_));
        REQUIRE( pgTest + 1.0f == pgComplex<float>(re_ + 1.0f, im_));
        REQUIRE( pgTest + pgJ == pgComplex<float>(re_, 1.0f + im_));
    }

    SECTION( "Test operator- " ) {
        REQUIRE( pgTest - pgZero == pgTest);
        REQUIRE( pgTest - pgTest == pgZero);
        REQUIRE( pgTest - 1.0f == pgComplex<float>(re_ - 1.0f, im_));
        REQUIRE( pgTest - pgJ == pgComplex<float>(re_, im_ - 1.0f));
    }

    SECTION( "Test operator* " ) {
        REQUIRE( pgTest * pgZero == pgZero);
        REQUIRE( pgTest * pgComplex<float>(1.0f, 0.0f) == pgTest);
        REQUIRE( pgTest * 1.0f == pgTest);
    }

    SECTION( "Test operator/ " ) {
        REQUIRE( pgZero  / pgTest == pgZero);
        REQUIRE( pgTest / pgComplex<float>(1.0f, 0.0f) == pgTest);
        REQUIRE( pgTest / 1.0f == pgTest);
    }
}
TEST_CASE("pgComplex<float>: operations against std::complex", "[pgComplex]") {
    // Setup Prerequsites for test
    float re_ = 1.2f;
    float im_ = -3.4f;
    
    pgComplex<float> pgTest(re_,im_); // 1.2f - j3.4
    pgComplex<float> pgZero = {};        // 0 + j0
    pgComplex<float> pgJ(0.0f, 1.0f); // 0 + j1

    std::complex<float> cplxTest(re_,im_); // 1.2f - j3.4
    std::complex<float> cplxZero = {};     // 0 + j0
    std::complex<float> cplxJ(0.0f,1.0f);  // 0 + j1

    SECTION( "Test exponential function " ) {
        REQUIRE(real(exp(pgZero)) == Approx(real(std::exp(cplxZero))));
        REQUIRE(real(exp(pgJ)) == Approx(real(std::exp(cplxJ))));
        REQUIRE(real(exp(pgTest)) == Approx(real(std::exp(cplxTest))));
        REQUIRE(imag(exp(pgZero)) == Approx(imag(std::exp(cplxZero))));
        REQUIRE(imag(exp(pgJ)) == Approx(imag(std::exp(cplxJ))));
        REQUIRE(imag(exp(pgTest)) == Approx(imag(std::exp(cplxTest))));
    }

    SECTION( "Test abs function " ) {
        REQUIRE(abs(pgZero) == Approx(abs(cplxZero)));
        REQUIRE(abs(pgTest) == Approx(abs(cplxTest)));
        REQUIRE(abs(pgJ) == Approx(abs(cplxJ)));
    }

    SECTION( "Test norm function " ) {
        REQUIRE(norm(pgZero) == Approx(norm(cplxZero)));
        REQUIRE(norm(pgTest) == Approx(norm(cplxTest)));
        REQUIRE(norm(pgJ) == Approx(norm(cplxJ)));
    }

    SECTION( "Test arg function " ) {
        REQUIRE(arg(pgZero) == Approx(arg(cplxZero)));
        REQUIRE(arg(pgTest) == Approx(arg(cplxTest)));
        REQUIRE(arg(pgJ) == Approx(arg(cplxJ)));
    }

}
