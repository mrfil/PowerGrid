# Base NVIDIA Image
#FROM nvidia/cuda:9.0-devel-ubuntu16.04 as builder
FROM nvcr.io/hpc/pgi-compilers:ce as builder
#FROM cuda-ssh:latest
# Offers X access via ssh for NVVP and PGProf
# Run everything in /root
WORKDIR /root

# Add needed packages
RUN  apt-get -y update && apt-get install -y g++ gcc curl wget libopenblas-dev \
												libarpack2-dev nano wget cmake git \
												libboost-all-dev libhdf5-dev libfftw3-dev \
												&& apt-get clean all;



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
	cmake ../ -DCMAKE_INSTALL_PREFIX=/opt; \
	make; \
	make install; \
	ldconfig

ENV LD_LIBRARY_PATH="/usr/local/lib:${PATH}"

# Install Latest CMake
RUN curl -O -J -L http://cmake.org/files/v3.9/cmake-3.9.4-Linux-x86_64.tar.gz; \
	apt-get -y purge cmake; \
	tar -xvf ./cmake-3.9.4-Linux-x86_64.tar.gz; \
	cd ./cmake-3.9.4-Linux-x86_64; \
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

ARG GITHUBTOKEN=983877d3d082851c3560ddbee074b124bdf46150

RUN git clone https://acerjanic:$GITHUBTOKEN@github.com/acerjanic/PrivatePG.git; \
    cd PrivatePG; \
    git fetch; \
    git checkout SharedLib; \
    mkdir build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=pgc++ -DCMAKE_INSTALL_PREFIX=/opt/PowerGrid; \
    make -j4; make install

ENV LD_LIBRARY_PATH="/opt/PowerGrid/lib:${LD_LIBRARY_PATH}"
ENV PATH="${PATH}:/opt/PowerGrid/bin"
	
#Stage 2 of build to create a runtime only version

#FROM nvidia/cuda:9.0-runtime-ubuntu16.04
#ARG PGI_VERSION=18.7
#ARG PGI_INSTALL_DIR=/usr/local/pgi
#COPY --from=builder /opt /opt
#COPY --from=builder /opt/lib/* /opt/lib/
#COPY --from=builder /usr/local/lib/* /opt/lib/
#
#RUN apt-get -y update && apt-get install -y libopenblas-base liblapack3 libarpack2 \
#											libhdf5-10 libhdf5-cpp-11 libfftw3-3 libatomic1 \
#											 && apt-get clean all;
ENV LD_LIBRARY_PATH="/opt/lib:${LD_LIBRARY_PATH}"





