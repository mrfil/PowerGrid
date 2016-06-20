---
layout: project
title: Creating Objects
category: docs
order_page: 3
---
## Creating Objects
{: .content-subhead }

Like the Image Reconstruction Toolbox, PowerGrid uses a similar approach to construct iterative reconstruction routines. The majority of the work is accomplished via objects that implement forward and adjoint operations. Performing a Fourier transform is achieved by using these objects, such as the Gfft or Gdft object.

Other forms of transforms, such as SENSE, sensitivity encoding, are also implemented. Furthermore, other transforms can be created by creating a new C++ class and implementing a forward and adjoint operator function.

### Creating the SENSE object

Let's look how the SENSE (Sensitivity Encoding) object for parallel imaging is implemented.

First we have boilerplate code called the include guard. This makes sure your code isn't included multiple times by accident and is standard practice in C/C++. Don't overthink this one.
```C++
#ifndef PowerGrid_SENSE_hpp
#define PowerGrid_SENSE_hpp

//Code goes here

#endif //PowerGrid_SENSE_hpp
```

We are going to use a lot of code from [Armadillo](http://arma.sourceforge.net) so let's include the namespace to shorten variable names

```C++
// Before: arma::mat<complex<float>>
using namespace arma;
// After: mat<complex<float>>
```
In PowerGrid, we use a very mild implementation of template metaprogramming to generalize classes and objects without resorting to complex C++ code. For example, the SENSE object needs to perform many forward and adjoint Fourier transforms. These can be field corrected or not, Gdft or Gnufft, etc. Templates let us wait until the code is compiled to worry about which type of variables or transform objects we are using. As long as all of the operators and functions called can be found at compile time, we can write the code in terms of a general variable type, called a templated or generic type.

```C++
template<typename T1, typename Tobj>
```

Our object class will use these types in place of the actual objects, keeping in mind that at compile time T1 and Tobj will be replaced by actual variable and object types.

Typically, T1 is a floating point variable type, float or double.

Tobj is a PowerGrid object, such as Gdft or Gnufft or TimeSegmentation.

From here we can define a typical class.

```C++
class SENSE {
    // define a helper type for a complex
    typedef complex <T1> CxT1;
public:
    SENSE();

    //Class variables go here
    uword n1 = 0; //Data size
    uword n2 = 0; //Image size
    uword nc = 0; //number of coils
    Tobj* G_obj;  //The transform object used to implement the recon.
    Mat <CxT1> SMap; //Coil sensitivity matrix with dimensions Image size b (n1 by number of coils (nc)

    //Class constructor
    SENSE(Tobj& G, Col <CxT1> SENSEmap, uword a, uword b, uword c)
    {
	    n1 = a;
	    n2 = b;
	    nc = c;
	    G_obj = &G;
	    SMap = reshape(SENSEmap, n2, nc);
    }
```

The most important part goes here. All PowerGrid objects, like their IRT/MATLAB counterparts, must implement a forward and adjoint transform operation. Due to some of the differences between operator precedence and overloading in C++ and MATLAB, we use operator*() to implement the forward transform and operator/() to implement the adjoint.

```C++
    //Overloaded operators go here

    //Forward transformation is *
    // d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col <CxT1> operator*(const Col <CxT1>& d) const
    {

	    Mat <CxT1> outData = zeros<Mat<CxT1 >> (this->n1, this->nc);
	    //Col<CxT1> temp;
	    //In SENSE we store coil data using the columns of the data matrix, and we weight the data by the coil sensitivities from the SENSE map

	    for (unsigned int ii = 0; ii<this->nc; ii++) {

		    outData.col(ii) = (*this->G_obj)*(d%(this->SMap.col(ii)));

	    }
	    Col <CxT1> out = vectorise(outData);
	    //equivalent to returning col(output) in MATLAB with IRT
	    return out;

    }

    //For the adjoint operation, we have to weight the adjoint transform of the coil data by the SENSE map.

    Col <CxT1> operator/(const Col <CxT1>& d) const
    {

	    Mat <CxT1> inData = reshape(d, this->n1, this->nc);

	    Col <CxT1> outData = zeros<Col<CxT1 >> (this->n2);
		  //Mat <CxT1> coilImages(n2,nc);

	    for (unsigned int ii = 0; ii<this->nc; ii++) {
			//coilImages.col(ii) = (*this->G_obj)/inData.col(ii);
		    outData += conj(this->SMap.col(ii))%((*this->G_obj)/inData.col(ii));
	    }
		  //outData = sum(conj(SMap)%coilImages,2);
	    //equivalent to returning col(output) in MATLAB with IRT
	    return vectorise(outData);


    }

```
Now we close the class and close the include guard.

```C++
};
#endif //PowerGrid_SENSE_hpp
```
