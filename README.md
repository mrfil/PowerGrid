# PowerGrid

Software for CPU and GPU accelerated iterative magnetic resonance imaging reconstruction. Quickly transate code from MATLAB/Image Reconstruction Toolbox to C++ code. Implements GPU accelerated non-uniform Fast Fourier Transforms and field correction via two algorithms. Also supports distributed memory computations via MPI for pcSENSE for correction of motion induced phase error.

## Depedenencies 
*   Xerces-C++ (https://xerces.apache.org/xerces-c/) - XML Parser used with the CodeSynthesis Generated files for parsing our config.xml files.

*   Armadillo (http://arma.sourceforge.net) - Templated Linear Algebra library. Gives us MATLAB like syntax in C++

*   MATIO  (http://matio.sourceforge.net) - Library for reading and writing MATLAB .mat files in C/C++

*   FFTW (http://www.fftw.org) - Fastest Fourier Transform in the West - Used for CPU implementations of the FFTs used in Gridding.

*   ISMRMRD (http://ismrmrd.github.io) ISMRM Raw Data Format - HDF5 Based data format for Magnetic Resonance Imaging data

