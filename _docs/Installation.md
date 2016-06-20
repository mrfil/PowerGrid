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

  * Linux (CPU, GPU, MPI) - Ubuntu 14.04 and 16.04 have been tested by the developers
  * Mac OS X (CPU, MPI) - El Capitan has been tested by the developers

#### Dependencies for CPU execution
 * Cmake - Version 2.8 or higher required.
 * libarmadillo[http://arma.sourceforge.net] - Version 6.x required, 6.700.7 recommended, Version 7.x not supported yet.
 * FFTW - We only support FFTW3 and need a reasonably recent version.
 * Xerces-C++ - Version 3.1.0 or greater.

#### Dependencies for Experimental ISMRMRD Support
 * libismrmrd - We are using the master branch in the repository but not version 2.0 yet.

#### Dependencies for GPU execution
 * cufft library - The version accompanying CUDA 7.0 has been tested
 * C++ compiler supporting [OpenACC](http://www.openacc.org)

We have experience with PGC++ 15.7 from [NVIDIA/The Portland Group](http://www.pgroup.com) as the version we have used most extensively. There is a free license available as part of the [OpenACC Toolkit](https://developer.nvidia.com/openacc-toolkit) for academic users.

GCC 6.1 has OpenACC support but has not yet been tested by the developers, we welcome reports of anyone trying to compile with it. We hope to support it alongside PGI compilers in the near future.

For those lucky enough to have access to Cray supercomputers, the Cray compiler does support OpenACC, but we have not tried to build with it. Because the Cray compilers are not available on desktops, workstations, or non-Cray branded clusters, we cannot dedicate resources to testing PowerGrid on it.

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
