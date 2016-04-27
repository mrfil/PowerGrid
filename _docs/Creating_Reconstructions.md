---
layout: project
title: Creating Reconstructions
category: docs
order_page: 3
---
## Creating Objects
{: .content-subhead }

Like the Image Reconstruction Toolbox, PowerGrid uses a similar approach to construct iterative reconstruction routines. The majority of the work is accomplished via objects that implement forward and adjoint operations. Performing a Fourier transform is achieved by using these objects, such as the Gfft or Gdft object.

Other forms of transforms, such as SENSE, sensitivity encoding, are also implemented. Furthermore, other transforms can be created by creating a new C++ class and implementing a forward and adjoint operator function.


### Creating the SENSE object

UNDER CONSTRUCTION
