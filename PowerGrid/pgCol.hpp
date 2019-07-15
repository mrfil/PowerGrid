// pgCol.hpp

#ifndef POWER_GRID_pgCol_hpp
#define POWER_GRID_pgCol_hpp

#include "PGIncludes.h"
#include "pgComplex.hpp"
#include <type_traits>

#ifdef _OPENACC
#include "openacc.h"
#include "accel.h"
#endif


template<typename T>
class pgCol {

private:

T *mem; //Pointer to raw data
bool isInitialized; 
bool isOnGPU; 
bool isCopy;
// Number of elements in array
//arma::uword n_elem;

public:

const arma::uword n_elem;

// Constructors

pgCol<T>() :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {  
        #pragma acc enter data copyin(this) 
    }


pgCol<T>(arma::uword length) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    #pragma acc enter data create(this) 
    set_size(length);
}


template <typename U = T>
pgCol<U>(arma::Col<T> &cSCplx, typename std::enable_if<std::is_integral<U>::value,
                                               U>::type dummy = 0) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    std::cout << "Entering pgCol<T> copy from Armadillo constructor" << std::endl;

    #pragma acc enter data create(this) 
    set_size(cSCplx.n_elem);
    memcpy(this->mem, cSCplx.memptr(), sizeof(T) * cSCplx.n_elem);
    
    #pragma acc update device(mem[0:n_elem])

}

template <typename U = T>
pgCol<U>(arma::Col<std::complex<float>> &cSCplx, typename std::enable_if<std::is_same<T, pgComplex<float>::value,
                                               int>::type dummy = 0) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    std::cout << "Entering pgCol<T> copy from Armadillo cplx constructor" << std::endl;

    #pragma acc enter data create(this) 
    set_size(cSCplx.n_elem);
    memcpy(this->mem, reinterpret_cast<float *>(cSCplx.memptr()), sizeof(T) * cSCplx.n_elem);

    #pragma acc update device(mem[0:n_elem])

}                                               

template <typename U = T>
pgCol<U>(arma::Col<std::complex<double>> &cSCplx, typename std::enable_if<std::is_same<T, pgComplex<double>::value,
                                               int>::type dummy = 0) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    std::cout << "Entering pgCol<T> copy from Armadillo cplx constructor" << std::endl;

    #pragma acc enter data create(this) 
    set_size(cSCplx.n_elem);
    memcpy(this->mem, reinterpret_cast<double *>(cSCplx.memptr()), sizeof(T) * cSCplx.n_elem);

    #pragma acc update device(mem[0:n_elem])

}                                               



// Copy Constructor
pgCol<T>(const pgCol<T>& pgA) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(true),
    n_elem(0) {
    #pragma acc enter data create(this) 
    set_size(pgA.n_elem);
    size_t bytes = sizeof(T) * pgA.n_elem;
    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    #pragma acc update device(this) 
    memcpy(mem, pgA.memptr(), bytes);
    #ifdef _OPENACC
        acc_memcpy(acc_deviceptr(mem), acc_deviceptr(pgA.memptr()), bytes); 
    #endif
    
}

// Move Constructor
pgCol<T>(pgCol<T>&& pgA) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    #pragma acc enter data create(this) 
    access::rw(n_elem) = pgA.n_elem;
    mem = pgA.memptr();
    isInitialized = true;
    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    
    #pragma acc update device(this) 
    #ifdef _OPENACC
        acc_attach((void **) &mem);
    #endif
    

    pgA.reset_mem();

}



// Default Destructor
~pgCol<T>() {

    #ifdef _OPENACC
        if( acc_deviceptr(mem) != NULL) {
            acc_delete((void *)mem, sizeof(T) * n_elem);
        }
    #endif
    if (mem != NULL) {
        delete[] mem;
    }

    #pragma acc exit data delete(this)

}

T* memptr() const {
    return mem;
}

void reset_mem() {
    #ifdef _OPENACC
        acc_detach((void **) &mem);
    #endif
    mem = NULL;
    isInitialized = false;
    access::rw(n_elem) = 0;
    isOnGPU = false;
    #pragma acc update device(this)

}

void set_size(arma::uword length) {
    if (isInitialized) {
        if (isOnGPU) {
            #pragma acc exit data finalize detach(mem) delete(mem[0:n_elem])
            isOnGPU  = false; 
        }
        delete[] mem;
        mem = NULL;
    }
    
    isInitialized = true;
    arma::access::rw(n_elem) = length;
    //mem = (T*)malloc(sizeof(T) * n_elem); 
    mem = new T[n_elem];
    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    #pragma acc update device(this)
    #pragma acc enter data create(mem[0:n_elem])

} 

void zeros() {
    #pragma acc parallel loop present(mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] = T(0.0);
    }
}

void ones() {
    #pragma acc parallel loop present(mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] = T(1.0);
    }
}

// Conversion from pgCol to arma::Col
template <typename U = T>
arma::Col<T> getArma(typename std::enable_if<std::is_integral<U>::value,
                                               U>::type dummy = 0) {
    #pragma acc update host(mem[0:n_elem])
    arma::Col<T> armaT(mem, n_elem, true, false);

    return armaT;
}
template <typename U = T>
arma::Col<std::complex<double>> getArma(typename std::enable_if<std::is_same<U, pgComplex<double>::value,
                                               int>::type dummy = 0) {
    #pragma acc update host(mem[0:n_elem])
    arma::Col<T> armaT(reinterpret_cast<float *>(mem), n_elem, true, false);

    return armaT;
}

template <typename U = T>
arma::Col<std::complex<float>> getArma(typename std::enable_if<std::is_same<T, pgComplex<float>::value,
                                               int>::type dummy = 0) {
    #pragma acc update host(mem[0:n_elem])
    arma::Col<T> armaT(reinterpret_cast<double *>(mem), n_elem, true, false);

    return armaT;
}
// Operators for element manipulation
// We'll assume .at() is for fast, GPU manipulation
inline
const T at(const arma::uword d) const {
    return mem[d];
}

inline
T& at(const arma::uword d) {
    return mem[d];
}

inline
T& operator()(const arma::uword d) {
    //#pragma acc update_host(mem)
    return mem[d];
}

inline
const T operator()(const arma::uword d) const {
    //#pragma acc update_host(mem)
    return mem[d];
}

pgCol<T>& operator=(const pgCol<T>& d) {
    #pragma acc parallel loop present(mem[0:n_elem],d)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] = d.at(ii);
    }
    return *this;
}

pgCol<T>& operator=(pgCol<T>&& d) {
    //size_t bytes = sizeof(T) * d.n_elem;
    if (isInitialized) {
        #pragma acc exit data delete(mem[0:n_elem])
        delete[] mem;
        access::rw(n_elem) = 0;
        isInitialized = false;
        isOnGPU = false;
    }
    access::rw(n_elem) = d.n_elem;
    isInitialized = true;
    isOnGPU = true;
    mem = d.memptr();
    #pragma acc update device(this)   
    #ifdef _OPENACC
        acc_attach((void **) &mem);
    #endif


    d.reset_mem();

    return *this;
}

pgCol<T>& operator+=(const T& A) {

    #pragma acc parallel loop present(mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] += A;
    }
    return *this;
}

pgCol<T>& operator-=(const T& A) {
    #pragma acc parallel loop present(this, mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] -= A;
    }
    return *this;
}

template<typename X>
pgCol<T>& operator%=(const T& A) {
    #pragma acc parallel loop present(this, mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] *= A;
    }
    return *this;
}

template<typename X>
pgCol<T>& operator/=(const T& A) {
    #pragma acc parallel loop present(this, mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] /= A;
    }
    return *this;
}

template<typename X>
pgCol<T>& operator+=(const pgCol<X> &pgA) {

    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] += pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgCol<T>& operator-=(const pgCol<X> &pgA) {
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] -= pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgCol<T>& operator%=(const pgCol<X> &pgA) {
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] *= pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgCol<T>& operator/=(const pgCol<X> &pgA) {
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] /= pgA.at(ii);
    }
    return *this;
}

// Operators(pgCol, scalar)
template<typename T>
const pgCol<T> operator+(const T& B) const {
    pgCol<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] + B;
    
    }
    return std::move(pgC);
}

template<typename T>
const pgCol<T> operator-( const T& B) const {
    pgCol<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] - B;
    }
    return std::move(pgC);
}

template<typename T>
const pgCol<T> operator%(const T& B) const {
    pgCol<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], SpgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] * B;
    }
    return std::move(pgC);
}

template<typename T>
const pgCol<T> operator/(const T& B) const {
    pgCol<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] / B;
    }
    return std::move(pgC);
}

// Operators(pgCol, pgCol)
template<typename X>
const pgCol<T> operator+(const pgCol<X>& pgB) const {
    pgCol<T> pgC(n_elem);

    #pragma acc parallel loop present( mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] + pgB.at(ii);
    }
    return std::move(pgC);
}

template<typename X>
const pgCol<T> operator-(const pgCol<X>& pgB) const {
    pgCol<T> pgC(n_elem);

    #pragma acc parallel loop present( mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] - pgB.at(ii);
    }
    return std::move(pgC);
}
template<typename X>
const pgCol<T> operator%(const pgCol<X>& pgB) const {
    pgCol<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] * pgB.at(ii);
    }

    return std::move(pgC);
}

template<typename X>
const pgCol<T> operator/(const pgCol<X>& pgB) const {
    pgCol<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] / pgB.at(ii);
    }
    return std::move(pgC);
}

};

template<typename T>
const pgComplex<T> sum(const pgCol<pgComplex<T>> &pgA) {
    T sumReal = {};
    T sumImag = {};

    #pragma acc parallel loop present(pgA) reduction(+:sumReal)
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumReal += real(pgA.at(ii));
    }

    #pragma acc parallel loop present(pgA) reduction(+:sumImag)
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumImag += imag(pgA.at(ii));
    }

    return pgComplex<T>(sumReal, sumImag);

}

template<typename T>
const T sum(const pgCol<T> &pgA) {
    T sumA = {};

    #pragma acc parallel loop present(pgA) reduction(+:sumA)
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumA += pgA.at(ii);
    }
    return sumA;
}

#endif //POWER_GRID_pgCol_hpp