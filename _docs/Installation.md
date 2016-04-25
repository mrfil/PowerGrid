---
layout: project
title: Installation
category: install
order_page: 1
---

## Building PowerGrid
{: .content-subhead }

Building PowerGrid requires several dependencies, depending on the features you want to use or develop. For testing or algorithm development without GPU acceleration, fewer dependencies are required. Testing with GPU support, distributed memory support (MPI), or ISMRMRD support requires additional dependencies.

### Required Dependencies

Supported OSes:

  * Linux (CPU, GPU, MPI)
  * Mac OS X (CPU, MPI)

#### Dependencies for CPU execution
 * Cmake
 * libarmadillo
 * FFTW
 * Xerces-C++

#### Dependencies for ISMRMRD Support
 * libismrmrd

#### Dependencies for GPU execution
 * cufft library
 * C++ compiler supporting [OpenACC](http://www.openacc.org)

We recommend PGC++ 15.7 from [NVIDIA/The Portland Group](http://www.pgroup.com) as the version we have used most extensively. There is a free license available as part of the [OpenACC Toolkit](https://developer.nvidia.com/openacc-toolkit)
 for academic users.

## Download PowerGrid
{: .content-subhead }

You can get a copy of PowerGrid either by cloning the Git repository or by downloading a copy of the repository.

### Clone Git Repository from GitHub

  git clone git@github.com:mrfil/PowerGrid.git

### Download archive from GitHub
To download PowerGrid, go to [PowerGrid GitHub page](https://github.com/mrfil/PowerGrid) and click "Download Zip" on the right side of the page.

If you would like to use git to host your website on GitHub Pages you can also click "Clone on Desktop" given that you've got git installed.

![Download Zip or Clone]( {{site.baseurl}}/assets/img/docs/github-download-clone.png)
{: .pure-image }

If you want to contribute code back to the project, we recommend forking the project.
