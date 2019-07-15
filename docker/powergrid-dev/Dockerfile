FROM nvcr.io/hpc/pgi-compilers:ce as builder

# Run everything in /root
WORKDIR /root

# Add needed packages
RUN apt-get -y update && apt-get install -y g++ gcc curl wget libopenblas-dev \
						libarpack2-dev nano wget git \
						libhdf5-dev libfftw3-dev && apt-get clean all;

#ENV LD_LIBRARY_PATH="/opt/pgi/linux86-64/2019/cuda/9.1/lib64:${LD_LIBRARY_PATH}"

# Install Latest CMake
RUN curl -O -J -L http://cmake.org/files/v3.14/cmake-3.14.0-Linux-x86_64.tar.gz; \
	apt-get -y purge cmake; \
	tar -xvf ./cmake-3.14.0-Linux-x86_64.tar.gz; \
	cd ./cmake-3.14.0-Linux-x86_64; \
	cp -r bin /usr/; \
	cp -r doc /usr/share/; \
	cp -r man /usr/share/; \
	cp -r share /usr/

	# Install custom compile of boost for mpi.
RUN wget -O ~/boost_1_70_0.tar.bz2 -c 'http://dl.bintray.com/boostorg/release/1.70.0/source/boost_1_70_0.tar.bz2' --max-redirect=100\
 && tar xjf ~/boost_1_70_0.tar.bz2 \
 && cd ~/boost_1_70_0/ \
 && ./bootstrap.sh --prefix=/usr/local --with-toolset=pgi --with-libraries=mpi,program_options,serialization \
 && cp ~/boost_1_70_0/tools/build/example/user-config.jam ~/boost_1_70_0/tools/build/src/user-config.jam \
 && echo "using mpi : /opt/pgi/linux86-64/2019/mpi/openmpi-3.1.3/bin/mpic++ ; " >> ~/boost_1_70_0/tools/build/src/user-config.jam

RUN cd ~/boost_1_70_0 && ./b2 -j 4 toolset=pgi install || true

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

RUN git clone https://github.com/acerjanic/ismrmrd.git; \
    cd ismrmrd; \
    git fetch; \
    git checkout matlab; \
    mkdir ./build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DBoost_FOUND=OFF; \
    make;  \
    make install

ENV LD_LIBRARY_PATH="/opt/PowerGrid/lib:${LD_LIBRARY_PATH}"
ENV PATH="${PATH}:/opt/PowerGrid/bin"

# For singularity compatibility
RUN ldconfig
