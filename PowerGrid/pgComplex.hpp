// pgComplex.hpp

#ifndef POWER_GRID_pgComplex_h
#define POWER_GRID_pgComplex_h

#include <cmath>
#include <sstream>
#include <complex>

template<typename T>
class pgComplex {

private:
    T __re_;
    T __im_;

public:
    pgComplex(const T& real = T(), const T& imag = T() )
        : __re_(real),
          __im_(imag) 
          {}

    template<typename X>            
    pgComplex(const pgComplex<X> &pgC)
        : __re_(pgC.real()),
          __im_(pgC.imag())
          {}

    template<typename X>   
    pgComplex(const std::complex<X> &sCplx) 
        : __re_(sCplx.real()),
          __im_(sCplx.imag())
         {}

    inline
    T real() const {
        return __re_;
    }

    inline
    T imag() const {
        return __im_;
    }

    inline
    void real(T _re) {
        __re_ = _re;
    }
    
    inline
    void imag(T _im) {
        __im_ = _im;
    }

    ~pgComplex() {};

    pgComplex& operator=(const T &_re) {
        __re_ = _re;
        __im_ = T();
        return *this;
    }  
      
    // In-place operations with real operator of same type
    inline
    pgComplex& operator+=(const T &_re) {
        __re_ += _re;
        return *this;
    }
    inline
    pgComplex& operator-=(const T &_re) {
        __re_ -= _re;
        return *this;
    }
    inline
    pgComplex& operator*=(const T &_re) {
        __re_ *= _re;
        __im_ *= _re;
        return *this;
    }
    inline
    pgComplex& operator/=(const T &_re) {
        __re_ /= _re;
        __im_ /= _re;
        return *this;
    }

    // In-place operations on pgComplex with different types

    template<typename X> 
    inline
    pgComplex& operator=(const pgComplex<X> pgC) {
        __re_ = pgC.real();
        __im_ = pgC.imag();
        return *this;
    }

    template<typename X> 
    inline
    pgComplex& operator-=(const pgComplex<X> pgC) {
        __re_ -= pgC.real();
        __im_ -= pgC.imag();
        return *this;
    }

    template<typename X> 
    pgComplex& operator+=(const pgComplex<X> pgC) {
        __re_ += pgC.real();
        __im_ += pgC.imag();
        return *this;
    }

    template<typename X> 
    inline
    pgComplex& operator*=(const pgComplex<X> pgC) {
        *this = *this * pgComplex<T>(pgC.real(), pgC.imag());
        return *this;
    }

    template<typename X> 
    inline
    pgComplex& operator/=(const pgComplex<X> pgC) {
        *this = *this / pgComplex<T>(pgC.real(), pgC.imag());
        return *this;
    }

};
// Non-member function templates
// Binary +,-,*,/ operators
// Defined in terms of in place operators
template<typename T>
inline
pgComplex<T> operator+(const pgComplex<T> &A, const pgComplex<T> &B) {
    pgComplex<T> pgC(A);
    pgC += B;
    return pgC;
}
template<typename T>
inline
pgComplex<T> operator+(const pgComplex<T> &A, const T &B) {
    pgComplex<T> pgC(A);
    pgC += B;
    return pgC;
}
template<typename T>
inline
pgComplex<T> operator+(const T &A, const pgComplex<T> &B) {
    pgComplex<T> pgC(B);
    pgC += A;
    return pgC;
}
template<typename T>
inline
pgComplex<T> operator-(const pgComplex<T> &A, const pgComplex<T> &B) {
    pgComplex<T> pgC(A);
    pgC -= B;
    return pgC;
}
template<typename T>
inline
pgComplex<T> operator-(const pgComplex<T> &A, const T &B) {
    pgComplex<T> pgC(A);
    pgC -= B;
    return pgC;
}
template<typename T>
inline
pgComplex<T> operator-(const T &A, const pgComplex<T> &B) {
    pgComplex<T> pgC(B);
    pgC -= A;
    return pgC;
}

template<typename T>
inline
pgComplex<T> operator*(const pgComplex<T> &pgA, const pgComplex<T> &pgB) {
    T re = pgA.real() * pgB.real() - pgA.imag() * pgB.imag();
    T im = pgA.real() * pgB.imag() + pgA.imag() * pgB.real();
    return pgComplex<T>(re, im);
}

template<typename T>
inline
pgComplex<T> operator*(const pgComplex<T> &pgA, const T &B) {
    pgComplex<T> pgC(pgA);
    pgC *= B;
    return pgC;
}

template<typename T>
inline
pgComplex<T> operator*(const T &A, const pgComplex<T> &pgB) {
    pgComplex<T> pgC(A);
    pgC *= pgB;
    return pgC;
}

template<typename T>
inline
pgComplex<T> operator/(const pgComplex<T> &pgA, const pgComplex<T> &pgB) {
    T denom = (pgB.real() * pgB.real() + pgB.imag() * pgB.imag());
    T re = (pgA.real() * pgB.real() + pgA.imag() * pgB.imag()) / denom; 
    T im = (pgA.imag() * pgB.real() - pgA.real() * pgB.imag()) / denom;
    return pgComplex<T>(re, im);
}

template<typename T>
inline
pgComplex<T> operator/(const pgComplex<T> &pgA, const T &B) {

    return pgComplex<T>(pgA.real() / B, pgA.imag() / B);
}

template<typename T>
inline
pgComplex<T> operator/(const T &A, const pgComplex<T> &pgB) {
    pgComplex<T> pgC(A);
    pgC /= pgB;
    return pgC;
}

// Unitary + and - 

template<typename T>
inline
pgComplex<T> operator+(const pgComplex<T> &pgA) {
    return pgA;
}

template<typename T>
inline
pgComplex<T> operator-(const pgComplex<T> &pgA) {
    return pgComplex<T>(-pgA.real(), -pgA.imag());
}

// Equals and Not Equals operators
template<typename T>
inline
bool operator==(const pgComplex<T> &pgA, const pgComplex<T> &pgB) {
    return (pgA.real() == pgB.real()) && (pgA.imag() == pgB.imag());
}

template<typename T>
inline
bool operator==(const pgComplex<T> &pgA, const T &B) {
    return (pgA.real() == B) && (pgA.imag() == 0);
}

template<typename T>
inline
bool operator==(const T &A, const pgComplex<T> &pgB) {
    return (A == pgB.real()) && (pgB.imag() == 0);
}

template<typename T>
inline
bool operator!=(const pgComplex<T> &pgA, const pgComplex<T> &pgB) {
    return !(pgA == pgB);
}

template<typename T>
inline
bool operator!=(const pgComplex<T> &pgA, const T &B) {
    return !(pgA == B);
}

template<typename T>
inline
bool operator!=(const T &A, const pgComplex<T> &pgB) {
    return !(A == pgB);
}

// Non-member functions real() and imag(): 
//        to support real(pgC) and real(T) syntax
template<typename T>
inline
T real(const pgComplex<T> &pgA){
    return pgA.real();
}

template<typename T>
inline
T real(const T &A){
    return T;
}

template<typename T>
inline
T imag(const pgComplex<T> &pgA){
    return pgA.imag();
}

template<typename T>
inline
T imag(const T &A){
    return 0;
}


// Non-member math functions
template<typename T>
inline
T abs(const pgComplex<T> &pgA) {
    return hypot(pgA.real(), pgA.imag());
}

template<typename T>
inline
T arg(const pgComplex<T> &pgA) {
    return atan2(pgA.imag(), pgA.real());
}

template<typename T>
inline
T arg(const T& A) {
    return atan2(T(), A);
}

template<typename T>
inline
T norm(const pgComplex<T> &pgA) {
    return pgA.real() * pgA.real() + pgA.imag() * pgA.imag();
}

template<typename T>
inline
pgComplex<T> conj(const pgComplex<T> &pgC) {
    return pgComplex<T>(pgC.real(), -pgC.imag());
}

template<typename T>
inline
pgComplex<T> polar(const T& rho, const T& theta = T()) {
    T x = rho * cos(theta);
    T y = rho * sin(theta);
    return pgComplex<T>(x ,y);
}

template<typename T>
inline
pgComplex<T> log(const pgComplex<T> &pgA) {
    return pgComplex<T>(log(abs(pgA)), arg(pgA));
}

template<typename T>
inline
pgComplex<T> log10(const pgComplex<T> &pgA) {
    return log(pgA)/log(T(10));
}

template<typename T>
inline
pgComplex<T> sqrt(const pgComplex<T> &pgA) {
    return polar( sqrt(abs(pgA)), arg(pgA) / T(2));
}

template<typename T>
inline
pgComplex<T> exp(const pgComplex<T> pgA) {
    T i = pgA.imag();
    T e = exp(pgA.real());
    return pgComplex<T>(e * cos(i), e * sin(i));
}

template<typename T>
inline
pgComplex<T> pow(const pgComplex<T> pgA, const pgComplex<T> pgB) {
    return exp(pgB * log(pgA));
}

// stream manipulators
template<typename T, typename charT>
std::basic_ostream<charT>& operator<<(std::basic_ostream<charT>& os, const pgComplex<T> pgA) {
    std::basic_ostringstream<charT> s;
    s.flags(os.flags());
    s.imbue(os.getloc());
    s.precision(os.precision());
    s << '(' << pgA.real() << ',' << pgA.imag() << ')';
    return os << s.str();
}

#endif //POWER_GRID_pgComplex_h
