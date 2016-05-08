---
layout: project
title: Objects in PowerGrid
category: docs
order_page: 1
---

## Object Reference
{: .content-subhead }


### Transform objects

These objects implement forward and adjoint transforms. All can be used via the syntax below. The forward operations and adjoint operations are implemented via operator overloading in C++.

```C++

arma::Col<complex<double>> data;

//G is a transform object

ForwardTransform = G * data;
AdjointTransform = G / data;
```

#### Gfft

Implements the FFT (Fast Fourier Transform) sampled uniformly on the grid in both dimensions. This transform supports both GPU and CPU computation, relying on cufft for GPU computation and FFTW for GPU computation.

Gfft(nx, ny, nz)

| Variables        | Description           |
| ------------- |:-------------|
| nx | image size in the x direction|
| ny | image size in the y direction|
| nz | number of slices   |


#### Gdft

Implements a field corrected discrete Fourier transform (DFT). This object supports both CPU and GPU computation, although the CPU implementation not recommended for production work. The CPU version is not multithreaded and is O(n^2) in complexity, scaling very poorly.

Gdft(sizeKspace, sizeImg, kx, ky, kz, ix, iy, iz, FieldMap, timeVec)

| Variables        | Description           |
| ------------- |:-------------|
| dataLength | k-space data size|
| sizeImg | image size |
| kx | k-space coordinates in the x direction    |
| ky | k-space coordinates in the y direction    |
| kz | k-space coordinates in the z direction    |
| ix | image coordinates in the x direction    |
| iy | image coordinates in the y direction    |
| iz | image coordinates in the z direction   |
| FieldMap | field map   |
| timeVec |  timing vector for a single shot (all shots assumed to have same timing vector)   |


#### Gnufft

Implements a non-Uniform Fast Fourier Transform (NUFFT). Field correction can be achieved with this object by combining with the TimeSegmentation object. This transform runs on both CPU and GPU, achieving good performance on both.

Gnufft(dataLength, gridOS, nx, ny, nz, k1, k2, k3, i1, i2, i3)

| Variables        | Description           |
| ------------- |:-------------|
| dataLength | length of k-space trajectory|
| gridOS | grid oversampling parameter |
| nx | image size in the x direction   |
| ny | image size in the x direction    |
| nz |number of slices   |
| kx | k-space coordinates in the x direction    |
| ky | k-space coordinates in the y direction    |
| kz | k-space coordinates in the z direction    |
| ix | image coordinates in the x direction    |
| iy | image coordinates in the y direction    |
| iz | image coordinates in the z direction   |


### Model objects

These objects allow to perform different types of image reconstructions. Forward and adjoint operations are also defined via operator overloading in C++.

#### TimeSegmentation

Implements corrections for image distortions due to magnetic field susceptibility using a time segmentation approach. Hanning interpolator and min-max formulation are implemented.

TimeSegmentation(G, FieldMap, timeVec, dataLength, sizeImg, L, interpType, shots)

| Variables        | Description   |
| ------------- |:-------------|
| G | Transform object such as TimeSegmentation|
| FieldMap | field map|
| dataLength | length of k-space trajectory  |
| sizeImg | image size |
| L | number of time segments|
| interpType | type of time interpolator (1 = Hanning interpolation ; 2= min-max interpolation). Default is 1|
| shots | number of shots. Default is 1|

#### SENSE

Implements a sensitivity encoding (SENSE) operator.

SENSE( G, SENSEmap, dataLength, sizeImg, nc)

| Variables        | Description           |
| ------------- |:-------------|
| G | Transform object such as Ggrid|
| SENSEmap | coil sensitivity map: dimensions are the image size (nx*ny*nz) by number of coils (nc)|
| dataLength | length of k-space trajectory   |
| sizeImg | image size |
| nc | number of coils|



#### pcSENSE

Implements SENSE operator with phase correction for motion induced phase errors in DTI.   

pcSENSE( kx, ky, kz, nx, ny, nz, nc, timeVec, SENSEmap, FieldMap, ShotPhaseMap)

| Variables        | Description           |
| ------------- |:-------------|
|    kx   | kspace coordinates in the x direction |
|    ky   | kspace coordinates in the y direction      |
|    kz   | kspace coordinates in the z direction      |
|    nx   | image size in the x direction |
|    ny   | image size in the y direction     |
|    nz    | number of slices      |
|    nc   | number of coils   |
|    timeVec   | timing vector for a single shot (all shots assumed to have same timing vector)     |
|    SENSEmap   | coil sensitivity map: dimensions are the image size (nx*ny*nz) by number of coils (nc) |
|    FieldMap   | field map      |
|    ShotPhaseMap   | shot phase map in radians: dimensions are image size (nx*ny*nz) by number of coils (nc)       |


#### mpipcSENSE

Implements the same operator than pcSENSE with MPI.

mpipcSENSE(kx, ky, kz, nx, ny, nz, nc, timeVec, SENSEmap, FieldMap, ShotPhaseMap, en, wor)

Same variables as pcSENSE. Two additional variables en , wor are boost::mpi variables for the MPI environment and MPI world communicator.

### Penalty objects

These objects correspond to different penalty functions that can be used in the solver.

#### Robject


#### TVPenalty

Implements a total variation penalty. It inherits from the Robject.

TVPenalty(nx, ny, nz, beta, delta)

| Variables        | Description           |
| ------------- |:-------------|
|    nx   | image size in the x direction |
|    ny   | image size in the y direction      |
|    nz   | number of slice      |
|    beta   | regularization parameter      |
|    delta   | regularization parameter      |

#### QuadPenalty

Implements a quadratic penalty. It inherits from the Robject class.

QuadPenalty(nx, ny, nz, beta);

| Variables        | Description           |
| ------------- |:-------------|
|    nx   | image size in the x direction |
|    ny   | image size in the y direction      |
|    nz   | number of slices      |
|    beta   | regularization parameter      |
