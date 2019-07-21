# PowerGrid

Software for CPU and GPU accelerated iterative magnetic resonance imaging reconstruction. Quickly transate code from MATLAB/Image Reconstruction Toolbox to C++0 code. Implements GPU accelerated non-uniform Fast Fourier Transforms and field correction via two algorithms. Also supports distributed memory computations via MPI for pcSENSE for correction of motion induced phase error.

## Depedenencies 
*   Armadillo (http://arma.sourceforge.net) - Templated Linear Algebra library. Gives us MATLAB like syntax in C++

*   ISMRMRD (http://ismrmrd.github.io) ISMRM Raw Data Format - HDF5 Based data format for Magnetic Resonance Imaging data

*   FFTW (http://www.fftw.org) - Fastest Fourier Transform in the West - Used for CPU implementations of the FFTs used in Gridding.

## Installing PowerGrid

### Using Docker (Recommended)

* 	Install Docker (https://docs.docker.com/install/linux/docker-ce/ubuntu)
* 	Install Nvidia-docker (https://github.com/NVIDIA/nvidia-docker)

Useful options include `-v /Source/On/Host:/Target/On/Container` to mount a directory of data and/or code onto the scanner.

```shell
docker pull mrfil/powergrid-dev
docker run --runtime=nvidia -it mrfil/powergrid-dev
```

### Installing dependencies on Ubuntu 16.04 (Not recommended - for advanced users only)

### Install PGI Community Edition compilers

You can find instructions for installing the PGI compilers with a free (as in beer) license at (https://www.pgroup.com/products/community.htm).

#### Add needed packages
```shell
sudo apt-get -y update && apt-get install -y g++ gcc curl wget libopenblas-dev \
						libarpack2-dev nano wget git libboost-all-dev \
						libhdf5-dev libfftw3-dev
```
#### Install Latest CMake
```shell
curl -O -J -L http://cmake.org/files/v3.14/cmake-3.14.0-Linux-x86_64.tar.gz
sudo apt-get -y purge cmake
tar -xvf ./cmake-3.14.0-Linux-x86_64.tar.gz
cd ./cmake-3.14.0-Linux-x86_64
cp -r bin /usr/
cp -r doc /usr/share/
cp -r man /usr/share/
cp -r share /usr/
```

#### Install SuperLU5
```shell
curl -O -J -L crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.2.1.tar.gz
tar xvf superlu_5.2.1.tar.gz
mkdir /root/SuperLU_5.2.1/build
cd /root/SuperLU_5.2.1/build
cmake ../ -DBUILD_SHARED_LIBS=ON
make
sudo make install
```

#### Install Armadillo
```shell
git clone https://gitlab.com/conradsnicta/armadillo-code
cd armadillo-code
git checkout 9.200.x
mkdir build
cd build
echo $PATH
cmake ../ -DCMAKE_INSTALL_PREFIX=/opt
make
sudo make install
sudo ldconfig
```


#### Install ISMRMRD 
```shell
git clone https://github.com/ismrmrd/ismrmrd.git
cd ismrmrd
git fetch
git checkout tags/v1.4.0
mkdir ./build
cd build
cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_INSTALL_PREFIX=/opt
make
sudo make install
```

#### Set up environment variables
```shell
export LD_LIBRARY_PATH="/opt/lib:/opt/lib64:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="/opt/PowerGrid/lib:${LD_LIBRARY_PATH}"
export PATH="${PATH}:/opt/PowerGrid/bin"
```

