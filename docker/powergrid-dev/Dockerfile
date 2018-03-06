#FROM nvidia/cuda:7.5-devel-ubuntu16.04
# Base NVIDIA Image
FROM nvidia/cuda:9.0-devel-ubuntu16.04
#FROM cuda-ssh:latest
# Offers X access via ssh for NVVP and PGProf
# Run everything in /root
WORKDIR /root

# Patch error in /etc/sshd_conf
# RUN  echo 'AddressFamily inet' >> /etc/ssh/sshd_config

# Add needed packages
RUN apt-get -y update && apt-get install -y g++ gcc curl wget libopenblas-dev \
																					libarpack2-dev nano wget cmake git \
																					xsdcxx libxerces-c-dev libmatio-dev \
																					libhdf5-dev \
																					libfftw3-dev && apt-get clean all;

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


# Install SuperLU5
RUN curl -O -J -L crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.2.1.tar.gz; \
	tar xvf superlu_5.2.1.tar.gz; \
	mkdir /root/SuperLU_5.2.1/build; \
	cd /root/SuperLU_5.2.1/build; \
	cmake ../ -DBUILD_SHARED_LIBS=ON; \
	make; \
	make install;

# Install PGI Compiler Community Edition (17.4)
# Download a copy of PGI 17.4 Community Edition and place in the same directory as the Dockerfile
#COPY pgilinux-2017-174-x86_64.tar.gz /root
#RUN mkdir /root/pgi; \
#    tar xvf /root/pgilinux-2017-174-x86_64.tar.gz -C /root/pgi;
#ENV PGI_SILENT true
#ENV PGI_ACCEPT_EULA accept
#ENV PGI_INSTALL_DIR /opt/pgi
#ENV PGI_INSTALL_NVIDIA true
#ENV PGI_INSTALL_AMD true
#ENV PGI_INSTALL_JAVA true
#ENV PGI_INSTALL_MPI true
#ENV PGI_MPI_GPU_SUPPORT true
#ENV PGI_INSTALL_MANAGED true

#RUN /root/pgi/install

#RUN echo "# PGI Compiler Stuff\n" >> /root/.bashrc; \
#	echo "PGI=/opt/pgi; export PGI\n" >> /root/.bashrc; \
# 	echo "MANPATH=\$MANPATH:\$PGI/linux86-64/17.4/man; export MANPATH\n" >> /root/.bashrc; \
#	echo "PATH=\$PATH:\$PGI/linux86-64/17.4/bin\n" >> /root/.bashrc

#ENV PGI /opt/pgi
#ENV PATH $PATH:$PGI/linux86-64/17.4/bin

RUN git clone https://github.com/Hopobcn/pgi-travis.git && \
    cd pgi-travis && \
    ./install-pgi.sh --nvidia --managed --mpi --mpi-gpu && \
    cd - && rm -rf pgi-travis

ENV PGI_VERSION 17.10
ENV PGI_INSTALL_DIR "/usr/local/pgi"
ENV PGI_HOME "${PGI_INSTALL_DIR}"/linux86-64/"${PGI_VERSION}"
ENV PGI_BIN_DIR "${PGI_HOME}"/bin
ENV PGI_LIB_DIR "${PGI_HOME}"/lib
ENV PGI_MAN_DIR "${PGI_HOME}"/man

ENV PATH            "${PGI_BIN_DIR}:${PATH}"
ENV LD_LIBRARY_PATH "${PGI_LIB_DIR}:${LD_LIBRARY_PATH}"
ENV MANPATH         "${PGI_MAN_DIR}:${MANPATH}"


# Install Armadillo

RUN git clone https://github.com/conradsnicta/armadillo-code; \
	cd armadillo-code; \
	git checkout 8.200.x; \
	mkdir build; \
	cd build; \
	echo $PATH; \
	cmake -DCMAKE_CXX_COMPILER=pgc++ ../; \
	make; \
	make install; \
	ldconfig

# Install Latest CMake

RUN curl -O -J -L http://cmake.org/files/v3.8/cmake-3.8.2-Linux-x86_64.tar.gz; \
	apt-get -y purge cmake; \
	tar -xvf ./cmake-3.8.2-Linux-x86_64.tar.gz; \
	cd ./cmake-3.8.2-Linux-x86_64; \
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
    cmake ../; \
    make;  \
    make install
