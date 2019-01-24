# Base NVIDIA Image
#FROM nvidia/cuda:9.0-devel-ubuntu16.04 as builder
FROM nvidia/cuda:9.0-devel-ubuntu16.04 AS pgi-ce

# PGI compiler version 18.10
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        libnuma1 \
        perl \
        wget && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -O /var/tmp/pgi-community-linux-x64-latest.tar.gz --referer https://www.pgroup.com/products/community.htm?utm_source=hpccm\&utm_medium=wgt\&utm_campaign=CE\&nvid=nv-int-14-39155 -P /var/tmp https://www.pgroup.com/support/downloader.php?file=pgi-community-linux-x64 && \
    mkdir -p /var/tmp/pgi && tar -x -f /var/tmp/pgi-community-linux-x64-latest.tar.gz -C /var/tmp/pgi -z && \
    cd /var/tmp/pgi && PGI_ACCEPT_EULA=accept PGI_INSTALL_DIR=/opt/pgi PGI_INSTALL_MPI=false PGI_INSTALL_NVIDIA=true PGI_MPI_GPU_SUPPORT=false PGI_SILENT=true ./install && \
    echo "variable LIBRARY_PATH is environment(LIBRARY_PATH);" >> /opt/pgi/linux86-64/18.10/bin/siterc && \
    echo "variable library_path is default(\$if(\$LIBRARY_PATH,\$foreach(ll,\$replace(\$LIBRARY_PATH,":",), -L\$ll)));" >> /opt/pgi/linux86-64/18.10/bin/siterc && \
    echo "append LDLIBARGS=\$library_path;" >> /opt/pgi/linux86-64/18.10/bin/siterc && \
    rm -rf /var/tmp/pgi-community-linux-x64-latest.tar.gz /var/tmp/pgi
ENV LD_LIBRARY_PATH=/opt/pgi/linux86-64/18.10/lib:$LD_LIBRARY_PATH \
    PATH=/opt/pgi/linux86-64/18.10/bin:$PATH

FROM pgi-ce as develop
WORKDIR /root

# Add needed packages
RUN  apt-get -y update && apt-get install -y curl libopenblas-dev \
										libarpack2-dev nano wget cmake git \
										libhdf5-dev libfftw3-dev && apt-get clean all;




# Instal MPICH 3.2 using example from NERSC
#RUN unset F90 && unset F90FLAGS && unset F77 && unset F77FLAGS && CPP=cpp && cd /root && wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz && tar xvzf mpich-3.2.tar.gz && cd /root/mpich-3.2 \
#	&& ./configure && make -j6 && make install && rm /root/mpich-3.2.tar.gz

# Install Open-MPI >= 2.1.0 for transparent MPI support via Singularity.
#RUN cd /root && wget https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.0.tar.gz && tar xzvf openmpi-2.1.0.tar.gz && cd /root/openmpi-2.1.0 \
#	&& env CC=pgcc CFLAGS=-O2 CPP=cpp CXX=pgc++ CXXFLAGS=-O2 -tp=bulldozer,sandybridge,px FC=pgfortran FCFLAGS=-O2 ./configure --enable-shared --enable-static --without-tm --enable-mpi-cxx && make -j4 && make install && rm /root/openmpi-2.1.0.tar.gz

RUN cd /root && wget https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.2.tar.gz && tar xzvf openmpi-2.1.2.tar.gz && cd /root/openmpi-2.1.2 \
	&& ./configure --enable-shared --enable-static --without-tm --enable-mpi-cxx && make -j4 && make install && rm /root/openmpi-2.1.2.tar.gz

ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

RUN cd openmpi-2.1.2/examples && make && cp ring_c /usr/local/bin
# Install custom compile of boost for mpi.
RUN wget -O ~/boost_1_66_0.tar.bz2 -c 'http://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.bz2' --max-redirect=100\
 && tar xjf ~/boost_1_66_0.tar.bz2 \
 && cd ~/boost_1_66_0/ \
 && ./bootstrap.sh --prefix=/usr/local --with-libraries=atomic,chrono,container,context,coroutine,date_time,exception,filesystem,graph,iostreams,locale,log,math,mpi,program_options,random,regex,serialization,signals,system,test,thread,timer,wave \
 && cp ~/boost_1_66_0/tools/build/example/user-config.jam ~/boost_1_66_0/tools/build/src/user-config.jam \
 && echo "using mpi ; " >> ~/boost_1_66_0/tools/build/src/user-config.jam

RUN cd ~/boost_1_66_0 && ./b2 -j 4 toolset=gcc install || true

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
	
ENV LD_LIBRARY_PATH="/opt/lib:${LD_LIBRARY_PATH}"

#RUN wget http://mvapich.cse.ohio-state.edu/download/mvapich/osu-micro-benchmarks-5.3.2.tar.gz; \
#    tar xf osu-micro-benchmarks-5.3.2.tar.gz; \
#    cd osu-micro-benchmarks-5.3.2; \
#    ./configure CC=$(which mpicc) CXX=$(which mpicxx); \
#    make; \
#    make install