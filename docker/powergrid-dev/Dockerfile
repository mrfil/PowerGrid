FROM nvcr.io/hpc/pgi-compilers:ce as builder

# Run everything in /root
WORKDIR /root

# Add needed packages
RUN apt-get -y update && apt-get install -y g++ gcc curl wget libopenblas-dev \
						libarpack2-dev nano wget git libboost-all-dev \
						libhdf5-dev libfftw3-dev && apt-get clean all;

ENV LD_LIBRARY_PATH="/opt/pgi/linux86-64/2018/cuda/9.1/lib64:${LD_LIBRARY_PATH}"

# Install Latest CMake
RUN curl -O -J -L http://cmake.org/files/v3.14/cmake-3.14.0-Linux-x86_64.tar.gz; \
	apt-get -y purge cmake; \
	tar -xvf ./cmake-3.14.0-Linux-x86_64.tar.gz; \
	cd ./cmake-3.14.0-Linux-x86_64; \
	cp -r bin /usr/; \
	cp -r doc /usr/share/; \
	cp -r man /usr/share/; \
	cp -r share /usr/

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
	cmake ../ -DCMAKE_INSTALL_PREFIX=/opt; \
	make; \
	make install; \
	ldconfig

ENV LD_LIBRARY_PATH="/opt/lib:/opt/lib64:${LD_LIBRARY_PATH}"

RUN git clone https://github.com/ismrmrd/ismrmrd.git; \
    cd ismrmrd; \
    git fetch; \
    git checkout tags/v1.4.0; \
    mkdir ./build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_INSTALL_PREFIX=/opt; \
    make;  \
    make install

ENV LD_LIBRARY_PATH="/opt/PowerGrid/lib:${LD_LIBRARY_PATH}"
ENV PATH="${PATH}:/opt/PowerGrid/bin"

# For singularity compatibility
RUN ldconfig
