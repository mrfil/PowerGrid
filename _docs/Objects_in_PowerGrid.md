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

#### Gdft

Implements a field corrected discrete Fourier transform (DFT). This object supports both CPU and GPU computation, although the CPU implementation not recommended for production work. The CPU version is not multithreaded and is O(n^2) in complexity, scaling very poorly.

#### Ggrid

Implements a non-Uniform Fast Fourier Transform (NUFFT). Field correction can be achieved with this object by combining with the TimeSegmentation object. This transform runs on both CPU and GPU, achieving good performance on both.

### Model objects

These objects allow to perform different types of image reconstructions. Forward and adjoint operations are also defined via operator overloading in C++.

#### TimeSegmentation

Implements corrections for image distortions due to magnetic field susceptibility using a time segmentation approach. Hanning interpolator and min-max formulation are implemented.

#### SENSE

Implements a sensitivity encoding (SENSE) operator.

#### pcSENSE

#### mpipcSENSE

### Penalty objects

These objects correspond to different penalty functions that can be used in the solver.

#### Robject


#### TVPenalty

Implements a total variation penalty. It inherits from the Robject.

#### QuadPenalty

Implements a quadratic penalty. It inherits from the Robject.
