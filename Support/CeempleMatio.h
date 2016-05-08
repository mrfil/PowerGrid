/*
 This file is distributed under a BSD 3-Clause license.
Copyright (c) 2014, Ceemple Software Ltd.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
 following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
 disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
 disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef CEEMPLEMATIO_H
#define CEEMPLEMATIO_H

#include "CeempleArmadillo.h"
#include "matio.h"
#include <iterator>
#include <algorithm>

using namespace arma;

template <class T>struct MatioTypes {
    static const matio_types DataType = MAT_T_UNKNOWN;
};
template <> struct MatioTypes<std::complex<double>> {
    static const matio_types DataType = MAT_T_DOUBLE;
    static const matio_flags FlagType = MAT_F_COMPLEX;
};
template <> struct MatioTypes<std::complex<float>> {
    static const matio_types DataType = MAT_T_SINGLE;
    static const matio_flags FlagType = MAT_F_COMPLEX;
};
template <> struct MatioTypes<uint8_t> {
    static const matio_types DataType = MAT_T_UINT8;
};
template <> struct MatioTypes<float> {
    static const matio_types DataType = MAT_T_SINGLE;
};
template <> struct MatioTypes<double> {
    static const matio_types DataType = MAT_T_DOUBLE;
};



// Return true on success.
template <typename T>
bool
savemat(const std::string &FileName, const std::string &VarName,
        const T &Var,
        const typename arma_not_cx<typename T::elem_type>::result* junk = 0) {
	//cout << "Calling real savemat" << endl;
    arma_ignore(junk);
    mat_t *matfp = Mat_CreateVer(FileName.c_str(), NULL, MAT_FT_MAT5);
    if (!matfp) {
        fprintf(stderr, "saveMat: could not create the file '%s'.\n",
                FileName.c_str());
        return false;
    }
	typedef typename T::elem_type Telem;
	const Telem* Vardata = Var.memptr();
	double* temp = new double[Var.n_rows*Var.n_cols];
	for (int ii = 0; ii<Var.n_rows*Var.n_cols; ii++) {
		temp[ii] = (double) Vardata[ii];
	}
    size_t dims[2] = {Var.n_rows, Var.n_cols};
    matvar_t *matvar =
            Mat_VarCreate(VarName.c_str(), MAT_C_DOUBLE, MatioTypes<typename T::elem_type>::DataType, 2,
		            dims, (void*) temp, 0);
    if (!matvar) {
        fprintf(stderr, "saveMat: error creating variable '%s'.\n.",
                VarName.c_str());
        Mat_Close(matfp);
        return false;
    }
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
    Mat_Close(matfp);
	delete[] temp;
    return true;
}

// Return true on success.

template<typename T1, typename std::enable_if<
		std::is_same<typename T1::elem_type, std::complex<double>> ::value, int>::type = 0>
inline
bool
savemat(const std::string &FileName, const std::string &VarName,
        const T1 &Var,
        const typename arma_cx_only<typename T1::elem_type>::result* junk = 0) {
    arma_extra_debug_sigprint();
    arma_ignore(junk);
	//cout << "Calling complex savemat" << endl;

    mat_t *matfp = Mat_CreateVer(FileName.c_str(), NULL, MAT_FT_MAT5);
    if (!matfp) {
        fprintf(stderr, "saveMat: could not create the file '%s'.\n",
                FileName.c_str());
        return false;
    }
    arma::Mat<cx_double> temp = conv_to<arma::Mat<cx_double>>::from(Var);
    size_t dims[2] = {temp.n_rows, temp.n_cols};
    arma::Mat<double> reTemp = arma::real(temp);
    arma::Mat<double> imTemp = arma::imag(temp);
    struct mat_complex_split_t z = {reTemp.memptr(),imTemp.memptr()};

	//std::cout << "MatioTypes<T1>::DataType = " << MatioTypes<typename T1::elem_type>::DataType << std::endl;

    matvar_t *matvar =
            Mat_VarCreate(VarName.c_str(), MAT_C_DOUBLE, MatioTypes<typename T1::elem_type>::DataType, 2,
		            dims, &z, MAT_F_COMPLEX);
    if (!matvar) {
        fprintf(stderr, "saveMat: error creating variable '%s'.\n.",
                VarName.c_str());
        Mat_Close(matfp);
        return false;
    }
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return true;
}

template<typename T1, typename std::enable_if<
		std::is_same<typename T1::elem_type, std::complex<float>> ::value, int>::type = 0>

inline
bool
savemat(const std::string& FileName, const std::string& VarName,
		const T1& Var,
		const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
{
	arma_extra_debug_sigprint();
	arma_ignore(junk);
	//cout << "Calling complex savemat" << endl;

	mat_t* matfp = Mat_CreateVer(FileName.c_str(), NULL, MAT_FT_MAT5);
	if (!matfp) {
		fprintf(stderr, "saveMat: could not create the file '%s'.\n",
				FileName.c_str());
		return false;
	}
	arma::Mat <std::complex<float>> temp = conv_to<arma::Mat<std::complex<float >>>::from(Var);
	size_t dims[2] = {temp.n_rows, temp.n_cols};
	arma::Mat<float> reTemp = arma::real(temp);
	arma::Mat<float> imTemp = arma::imag(temp);
	float* reData = (float*) reTemp.memptr();
	float* imData = (float*) imTemp.memptr();
	double* reDataOut = (double*) malloc(temp.n_rows*temp.n_cols*sizeof(double));
	double* imDataOut = (double*) malloc(temp.n_rows*temp.n_cols*sizeof(double));

	//Cast data from double to desired type
	for (int ii = 0; ii<temp.n_rows*temp.n_cols; ii++) {
		reDataOut[ii] = (double) (reData[ii]);
		imDataOut[ii] = (double) (imData[ii]);
	}

	struct mat_complex_split_t z = {reDataOut, imDataOut};

	//std::cout << "MatioTypes<T1>::DataType = " << MatioTypes<typename T1::elem_type>::DataType << std::endl;

	matvar_t* matvar =
	Mat_VarCreate(VarName.c_str(), MAT_C_DOUBLE, MatioTypes<std::complex<double>>
	::DataType, 2,
			dims, &z, MAT_F_COMPLEX);
	if (!matvar) {
		fprintf(stderr, "saveMat: error creating variable '%s'.\n.",
				VarName.c_str());
		Mat_Close(matfp);
		return false;
	}
	Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
	Mat_VarFree(matvar);
	Mat_Close(matfp);
	free(reDataOut);
	free(imDataOut);
	return true;
}

template<typename T1>
bool
loadmat(const std::string &FileName, const std::string &VarName,
		T1* Var,
		const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
{
    arma_ignore(junk);
    cout << "Calling real loadmat" << endl;
	typedef typename T1::elem_type Telems;
    mat_t *matfp = Mat_Open(FileName.c_str(), MAT_ACC_RDONLY);
    if (!matfp) {
        fprintf(stderr, "loadMat: could not open the file '%s'.\n",
                FileName.c_str());
        return false;
    }
    matvar_t *matvar = Mat_VarRead(matfp, VarName.c_str());
    if (!matvar) {
        fprintf(stderr, "loadMat: variable '%s' not found in file '%s'.\n",
                VarName.c_str(), FileName.c_str());
        Mat_Close(matfp);
        return false;
    }
    bool matvar2D = (matvar->rank == 2) && (matvar->class_type == MAT_C_DOUBLE);
    if (matvar2D) {
        unsigned Rows = matvar->dims[0], Cols = matvar->dims[1];

	    Telems* dataOut = (Telems*) malloc(Rows*Cols*sizeof(Telems));
	    double* Tdata = (double*) matvar->data;
	    //Cast data from double to desired type
	    for (int ii = 0; ii<Rows*Cols; ii++) {
		    dataOut[ii] = (Telems) Tdata[ii];
	    }

	    arma::Mat <Telems> out(dataOut, Rows, Cols, true, true);
	    (*Var) = conv_to<T1>::from(out);
	    //memcpy(Var->memptr(), matvar->data, Rows * Cols * sizeof(double));
	    free(dataOut);
    } else {
        fprintf(stderr,
                "loadMat: Variable '%s' is not 2-dimensional double matrix.\n"
                        "rank =** %d class_type = %d\n",
                VarName.c_str(), matvar->rank, matvar->class_type);
    }
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return matvar2D;
}

template<typename T1, typename std::enable_if<
		std::is_same<typename T1::elem_type, std::complex<double>> ::value, int>::type = 0>
bool
loadmat(const std::string &FileName, const std::string &VarName,
		T1* Var, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
{
    arma_ignore(junk);
    cout << "Calling complex loadmat" << endl;

    mat_t *matfp = Mat_Open(FileName.c_str(), MAT_ACC_RDONLY);
    if (!matfp) {
        fprintf(stderr, "loadMat: could not open the file '%s'.\n",
                FileName.c_str());
        return false;
    }
    matvar_t *matvar = Mat_VarRead(matfp, VarName.c_str());
    if (!matvar) {
        fprintf(stderr, "loadMat: variable '%s' not found in file '%s'.\n",
                VarName.c_str(), FileName.c_str());
        Mat_Close(matfp);
        return false;
    }
    bool matvar2D = (matvar->rank == 2) && (matvar->class_type == MAT_C_DOUBLE);
    if (matvar2D) {
        unsigned Rows = matvar->dims[0], Cols = matvar->dims[1];
        //(*Var) = arma::Mat<T1>(Rows, Cols);

        //If complex data, we get a pointer to a mat_complex_split_t where we can get real and imaginary pointers
        mat_complex_split_t* z = (mat_complex_split_t*)matvar->data;

        arma::Mat<double> re((double*)z->Re,Rows,Cols,false,true);
        arma::Mat<double> im((double*)z->Im,Rows,Cols,false,true);
        arma::Mat<arma::cx_double> out(re,im);
        (*Var) = conv_to<T1>::from(out);
        //memcpy(Var->memptr(), matvar->data, Rows * Cols * sizeof(T));
    }
    else {
	    fprintf(stderr,
			    "loadMat: Variable '%s' is not 2-dimensional double matrix.\n"
					    "rank = %d class_type = %d\n",
			    VarName.c_str(), matvar->rank, matvar->class_type);
    }
	Mat_VarFree(matvar);
	Mat_Close(matfp);
	return matvar2D;
}

template<typename T1, typename std::enable_if<
		std::is_same<typename T1::elem_type, std::complex<float>> ::value, int>::type = 0>

bool
loadmat(const std::string& FileName, const std::string& VarName,
		T1* Var, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
{
	arma_ignore(junk);
	cout << "Calling complex loadmat" << endl;

	mat_t* matfp = Mat_Open(FileName.c_str(), MAT_ACC_RDONLY);
	if (!matfp) {
		fprintf(stderr, "loadMat: could not open the file '%s'.\n",
				FileName.c_str());
		return false;
	}
	matvar_t* matvar = Mat_VarRead(matfp, VarName.c_str());
	if (!matvar) {
		fprintf(stderr, "loadMat: variable '%s' not found in file '%s'.\n",
				VarName.c_str(), FileName.c_str());
		Mat_Close(matfp);
		return false;
	}
	bool matvar2D = (matvar->rank==2) && (matvar->class_type==MAT_C_DOUBLE);
	if (matvar2D) {
		unsigned Rows = matvar->dims[0], Cols = matvar->dims[1];
		//(*Var) = arma::Mat<T1>(Rows, Cols);

		//If complex data, we get a pointer to a mat_complex_split_t where we can get real and imaginary pointers
		mat_complex_split_t* z = (mat_complex_split_t*) matvar->data;
		double* reData = (double*) z->Re;
		double* imData = (double*) z->Im;

		float* reDataOut = (float*) malloc(Rows*Cols*sizeof(float));
		float* imDataOut = (float*) malloc(Rows*Cols*sizeof(float));

		//Cast data from double to desired type
		for (int ii = 0; ii<Rows*Cols; ii++) {
			reDataOut[ii] = (float) (reData[ii]);
			imDataOut[ii] = (float) (imData[ii]);
		}

		arma::Mat<float> re(reDataOut, Rows, Cols, true, true);
		arma::Mat<float> im(imDataOut, Rows, Cols, true, true);
		arma::Mat <std::complex<float>> out(re, im);
		(*Var) = conv_to<T1>::from(out);
		//memcpy(Var->memptr(), matvar->data, Rows * Cols * sizeof(T));
		free(reDataOut);
		free(imDataOut);
    } else {
        fprintf(stderr,
                "loadMat: Variable '%s' is not 2-dimensional double matrix.\n"
                        "rank = %d class_type = %d\n",
                VarName.c_str(), matvar->rank, matvar->class_type);
    }
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return matvar2D;
}

#endif // CEEMPLEMATIO_H//
