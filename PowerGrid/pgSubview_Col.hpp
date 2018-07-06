// pgCol.hpp

#include "PGIncludes.h"

template<typename T>
class pgCol<T> {

private:

T *mem; //Pointer to raw data
bool isInitialized; 
bool isOnGPU; 

// Number of elements in array
unsigned int n_elem;

public:

void set_size() {
    if (isInitialized {
        if (isOnGPU) {
            #pragma acc exit data delete(mem) 
        }
        delete[](mem);
        mem = malloc(); 
        #pragma acc enter data create(mem[0:n_elem])
        #ifdef _OPENACC
        acc_attach();
        #endif
    }
} 

unsigned int n_elem() const{
    return n_elem;
}

void zeros() {
    #pragma acc parallel loop present(mem[0:n_elem])
    for(unsigned int ii = 0; ii < n_elem; ii++) {
        mem[ii] = T();
    }
}

void ones() {
    #pragma acc parallel loop present(mem[0:n_elem])
    for(unsigned int ii = 0; ii < n_elem; ii++) {
        mem[ii] = T(1.0);
    }
}

operator=(const pgCol<T>& d) {
    if()
}

operator()(const unsigned int d) {
    return
}

const operator()(const unsigned int d) const {
    return
}

};


