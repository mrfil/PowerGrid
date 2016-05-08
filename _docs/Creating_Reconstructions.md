---
layout: project
title: Creating Reconstructions
category: docs
order_page: 3
---
## Performing an Image Reconstruction

Like the Image Reconstruction Toolbox, PowerGrid uses a similar approach to construct iterative reconstruction routines. A reconstruction routine will leverage the use of objects to perform matrix multiplications. This section will present an example of how the objects can be combined in order to create a reconstruction routine.  The example will use several models and one solver but could be modified and extended to use different solvers and reconstruction models.  Please see the documentation on objects and solver routines to see more provided options.

### Example 1: Single Coil Image Reconstruction

An example of code that can be used to perform a reconstruction using the a discrete Fourier transform with data from a single coil is demonstrated.

The first step is to declare the main variables that the reconstruction routine will use

```C++
//Setup image space coordinates/trajectory

//Pixel locations in image space
Col<double> ix;
Col<double> iy;
Col<double> iz;
//location of data in kspace
Col<double> kx;
Col<double> ky;
Col<double> kz;

//Field correction parameters
//(can use a vector of zeros if field correction is not desired)
Col<double> FM;  //field map in radians/second
Col<double> tvec; //Time each data point was collected in seconds

//define the size of the reconstruction in image and kspace
uword Nk; //number of kspace data points
uword Nx; //Size of image in the x direction
uword Ny; //Size of image in the y direction
uword Nz; //Size of image in the z direction

//declare data which is complex valued
Col<complex<double>> data;
```

The next step to is to fill the declare variables with the correct information.  This can be achieved through loading data from an external file or header.   Code that creates the correct information can also be used.

An example of code that defines image space coordinates on a cartesain grid is given below.


```C++
for(uword ii = 0; ii < Ny; ii++) { //y
    for (uword jj = 0; jj < Nx; jj++) { //x
        for (uword kk = 0; kk < Nz; kk++) { //z
            ix(ii+ jj*Ny+ kk*Ny*Nx) = ((double) jj - (T1) Nx / 2.0) / ((double) Nx);
            iy(ii+ jj*Ny+ kk*Ny*Nx) = ((double) ii - (T1) Ny / 2.0) / ((double) Ny);
            iz(ii+ jj*Ny+ kk*Ny*Nx) = ((double) kk - (T1) Nz / 2.0) / ((double) Nz);
        }
    }
}
```
The object model that is desired to be used in the reconstruction can now be created

```c++
Gdft<double> Gd(Nk,Nx*Ny*Nz,kx,ky,kz,ix,iy,iz,FM,tvec);
```
Most currently written solvers support the use of regularization penalties.

Here is an example of how a quadratic penalty object could be created.

```C++
double beta;
beta = 0.1;
QuadPenalty < double > R(Nx, Ny, Nz, beta);
```

All of the created objects along with the kspace data can now be passed to a solver routine to a solver routine to perform the image reconstruction.

```C++
Col<complex<double>> xinit(Nx*Ny*Nz); // initial estimate of x
xinit.zeros();
Col < double > W; //is a data weighting term
W.ones(Nk);
uword niter;
niter = 30; //number of CG iterations to peform

Col<complex<double>> out_image;
out_image = solve_pwls_pcg<double,Gdft<double>,QuadPenalty<double>>(xinit, Gd, W, data, R, niter);
```
### Example 2: Incorporating SENSE

The code above can be modified to use coil sensitiivity encoding

A coil sensitivity map is now needed

```C++
uword Nc;
Nc = 4; //number of coils
Col<complex<double>> SMap;
//fill the SENSE map with the coil sensitivity information

Col<complex<double>> data;
//data must now be define as the data from all the coils
```
The Gdft object can then be wrapped by a  SENSE object

```C++
 SENSE<double, Gdft<double>> Sd(Gd,SMap,Nk,Nx*Ny*Nz,Nc);
```
The modified new object can now be used by the solver

```C++
out_image = solve_pwls_pcg<double,SENSE<double, Gdft<double>>,QuadPenalty<double>>(xinit, Sd, W, data, R, niter);
```
