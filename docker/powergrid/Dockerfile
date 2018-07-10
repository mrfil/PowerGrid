# Base NVIDIA Image
FROM nvidia/cuda:9.0-devel-ubuntu16.04 as builder
#FROM cuda-ssh:latest
# Offers X access via ssh for NVVP and PGProf
# Run everything in /root
WORKDIR /root

# Add needed packages
RUN apt-get -y update && apt-get install -y g++ gcc curl wget libopenblas-dev \
												libarpack2-dev nano wget cmake git \
												xsdcxx libxerces-c-dev libmatio-dev \
												libhdf5-dev libfftw3-dev && apt-get clean all;

# Set MPI library to MPICH2. Needed for Shifter to transparently patch MPI for Cray libraries
#RUN update-alternatives --set mpi /usr/include/mpich;

# Instal MPICH 3.2 using example from NERSC
RUN cd /root && wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz && tar xvzf mpich-3.2.tar.gz && cd /root/mpich-3.2 \
	&& ./configure && make -j4 && make install && make clean && rm /root/mpich-3.2.tar.gz

# Install custom compile of boost for mpi.
RUN wget -O ~/boost_1_66_0.tar.bz2 -c 'http://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.bz2' --max-redirect=100\
 && tar xjf ~/boost_1_66_0.tar.bz2 \
 && cd ~/boost_1_66_0/ \
 && ./bootstrap.sh --prefix=/usr/local --with-libraries=atomic,chrono,container,context,coroutine,date_time,exception,filesystem,graph,iostreams,locale,log,math,mpi,program_options,random,regex,serialization,signals,system,test,thread,timer,wave \
 && cp ~/boost_1_66_0/tools/build/example/user-config.jam ~/boost_1_66_0/tools/build/src/user-config.jam \
 && echo "using mpi ; " >> ~/boost_1_66_0/tools/build/src/user-config.jam

RUN cd ~/boost_1_66_0 && ./b2 -j 4 toolset=gcc install || true

RUN mkdir /opt/lib
ENV LD_LIBRARY_PATH="/opt/lib:${LD_LIBRARY_PATH}"
# Install SuperLU5
RUN curl -O -J -L crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.2.1.tar.gz; \
	tar xvf superlu_5.2.1.tar.gz; \
	mkdir /root/SuperLU_5.2.1/build; \
	cd /root/SuperLU_5.2.1/build; \
	cmake ../ -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=/opt; \
	make; \
	make install;

# Install PGI Compiler Community Edition
RUN git clone https://github.com/Hopobcn/pgi-travis.git && \
    cd pgi-travis && \
    ./install-pgi.sh --nvidia --managed --mpi --mpi-gpu && \
    cd - && rm -rf pgi-travis

ARG PGI_VERSION=18.4
ARG PGI_INSTALL_DIR=/usr/local/pgi

ENV PGI_VERSION ${PGI_VERSION}
ENV PGI_INSTALL_DIR ${PGI_INSTALL_DIR}
ENV PGI_HOME ${PGI_INSTALL_DIR}/linux86-64/${PGI_VERSION}
ENV PGI_BIN_DIR ${PGI_HOME}/bin
ENV PGI_LIB_DIR ${PGI_HOME}/lib
ENV PGI_MAN_DIR ${PGI_HOME}/man

ENV PATH="${PGI_BIN_DIR}:${PATH}"
ENV LD_LIBRARY_PATH="${PGI_LIB_DIR}:${LD_LIBRARY_PATH}"
ENV MANPATH="${PGI_MAN_DIR}:${MANPATH}"

RUN cp -r --dereference ${PGI_INSTALL_DIR}/linux86-64/${PGI_VERSION}/lib/* /opt/lib/

# Install Armadillo
RUN git clone https://github.com/conradsnicta/armadillo-code; \
	cd armadillo-code; \
	git checkout 8.500.x; \
	mkdir build; \
	cd build; \
	echo $PATH; \
	cmake ../ -DCMAKE_INSTALL_PREFIX=/opt; \
	make; \
	make install; \
	ldconfig

# Install Latest CMake
RUN curl -O -J -L http://cmake.org/files/v3.9/cmake-3.9.4-Linux-x86_64.tar.gz; \
	apt-get -y purge cmake; \
	tar -xvf ./cmake-3.9.4-Linux-x86_64.tar.gz; \
	cd ./cmake-3.9.4-Linux-x86_64; \
	cp -r bin /usr/; \
	cp -r doc /usr/share/; \
	cp -r man /usr/share/; \
	cp -r share /usr/

RUN git clone https://github.com/acerjanic/ismrmrd.git; \
    cd ismrmrd; \
    git fetch; \
    git checkout matlab; \
    mkdir ./build; \
    cd build; \
    cmake ../ -DCMAKE_INSTALL_PREFIX=/opt; \
    make;  \
    make install

ARG GITHUBTOKEN=983877d3d082851c3560ddbee074b124bdf46150

RUN git clone https://acerjanic:$GITHUBTOKEN@github.com/acerjanic/PrivatePG.git; \
    cd PrivatePG; \
    git fetch; \
    git checkout master; \
    mkdir build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=pgc++ -DCMAKE_INSTALL_PREFIX=/opt/PowerGrid; \
    make;  \
    make install
	
#Stage 2 of build to create a runtime only version

FROM nvidia/cuda:9.0-runtime-ubuntu16.04
ARG PGI_VERSION=18.4
ARG PGI_INSTALL_DIR=/usr/local/pgi
COPY --from=builder /opt /opt
COPY --from=builder /opt/lib/* /opt/lib/
COPY --from=builder /usr/local/lib/* /opt/lib/

RUN apt-get -y update && apt-get install -y libopenblas-base liblapack3 libarpack2 \
											libhdf5-10 libhdf5-cpp-11 libfftw3-3 libatomic1 \
											 && apt-get clean all;
ENV LD_LIBRARY_PATH="/opt/lib:${LD_LIBRARY_PATH}"





