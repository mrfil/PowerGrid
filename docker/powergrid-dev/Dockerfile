#FROM nvidia/cuda:7.5-devel-ubuntu16.04
# Base NVIDIA Image
#FROM nvidia/cuda:9.0-devel-ubuntu16.04
FROM nvcr.io/hpc/pgi-compilers:ce as builder
#FROM cuda-ssh:latest
# Offers X access via ssh for NVVP and PGProf
# Run everything in /root
WORKDIR /root

# Patch error in /etc/sshd_conf
# RUN  echo 'AddressFamily inet' >> /etc/ssh/sshd_config

ENV LD_LIBRARY_PATH="/opt/pgi/linux86-64/2018/cuda/9.1/lib64:${LD_LIBRARY_PATH}"

# Add needed packages
RUN apt-get -y update && apt-get install -y g++ gcc curl wget libopenblas-dev \
						libarpack2-dev nano wget cmake git libboost-all-dev \
						libhdf5-dev libfftw3-dev && apt-get clean all;

ENV LD_LIBRARY_PATH="/opt/pgi/linux86-64/2018/cuda/9.1/lib64:${LD_LIBRARY_PATH}"

# Install SuperLU5
RUN curl -O -J -L crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.2.1.tar.gz; \
	tar xvf superlu_5.2.1.tar.gz; \
	mkdir /root/SuperLU_5.2.1/build; \
	cd /root/SuperLU_5.2.1/build; \
	cmake ../ -DBUILD_SHARED_LIBS=ON; \
	make; \
	make install;

# Install Armadillo

RUN git clone https://gitlab.com/conradsnicta/armadillo-code; \
	cd armadillo-code; \
	git checkout 9.200.x; \
	mkdir build; \
	cd build; \
	echo $PATH; \
	cmake ../; \
	make; \
	make install; \
	ldconfig

ENV LD_LIBRARY_PATH="/usr/local/lib:${PATH}"


# Install Latest CMake

# Install Latest CMake
RUN curl -O -J -L http://cmake.org/files/v3.14/cmake-3.14.0-Linux-x86_64.tar.gz; \
	apt-get -y purge cmake; \
	tar -xvf ./cmake-3.14.0-Linux-x86_64.tar.gz; \
	cd ./cmake-3.14.0-Linux-x86_64; \
	cp -r bin /usr/; \
	cp -r doc /usr/share/; \
	cp -r man /usr/share/; \
	cp -r share /usr/


RUN git clone https://github.com/ismrmrd/ismrmrd.git; \
    cd ismrmrd; \
    git fetch; \
    git checkout master; \
    mkdir ./build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc; \
    make;  \
    make install
ENV LD_LIBRARY_PATH="/opt/PowerGrid/lib:${LD_LIBRARY_PATH}"
ENV PATH="${PATH}:/opt/PowerGrid/bin"