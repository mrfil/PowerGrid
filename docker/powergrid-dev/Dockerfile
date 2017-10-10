#FROM nvidia/cuda:7.5-devel-ubuntu16.04
# Base NVIDIA Image
FROM cuda-ssh:latest
# Offers X access via ssh for NVVP and PGProf
# Run everything in /root
WORKDIR /root

# Patch error in /etc/sshd_conf
RUN  echo 'AddressFamily inet' >> /etc/ssh/sshd_config

# Add needed packages
RUN apt-get -y update && apt-get install -y g++ gcc curl libopenblas-dev libarpack2-dev nano cmake git xsdcxx libxerces-c-dev libmatio-dev libboost-all-dev libhdf5-dev libfftw3-dev

# Install SuperLU5
RUN curl -O -J -L crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.2.1.tar.gz; \
	tar xvf superlu_5.2.1.tar.gz; \
	mkdir /root/SuperLU_5.2.1/build; \
	cd /root/SuperLU_5.2.1/build; \
	cmake ../ -DBUILD_SHARED_LIBS=ON; \
	make; \
	make install;

# Install PGI Compiler Community Edition (17.4)

COPY pgilinux-2017-174-x86_64.tar.gz /root
RUN mkdir /root/pgi
RUN tar xvf /root/pgilinux-2017-174-x86_64.tar.gz -C /root/pgi
ENV PGI_SILENT true
ENV PGI_ACCEPT_EULA accept
ENV PGI_INSTALL_DIR /opt/pgi
ENV PGI_INSTALL_NVIDIA true
ENV PGI_INSTALL_AMD true
ENV PGI_INSTALL_JAVA true
ENV PGI_INSTALL_MPI true
ENV PGI_MPI_GPU_SUPPORT true
ENV PGI_INSTALL_MANAGED true

RUN /root/pgi/install

RUN echo "# PGI Compiler Stuff\n" >> /root/.bashrc; \
	echo "PGI=/opt/pgi; export PGI\n" >> /root/.bashrc; \
 	echo "MANPATH=\$MANPATH:\$PGI/linux86-64/17.4/man; export MANPATH\n" >> /root/.bashrc; \
	echo "PATH=\$PATH:\$PGI/linux86-64/17.4/bin\n" >> /root/.bashrc

ENV PGI /opt/pgi
ENV PATH $PATH:$PGI/linux86-64/17.4/bin

# Install Armadillo

RUN curl -O -J -L  https://sourceforge.net/projects/arma/files/armadillo-7.800.2.tar.xz/download; \
	tar xvf armadillo-7.800.2.tar.xz; \
	cd armadillo-7.800.2; \
	mkdir build; \
	cd build; \
	echo $PATH; \
	cmake -DCMAKE_CXX_COMPILER=pgc++ ../; \
	make; \
	make install;

	# Install Latest CMake

	RUN curl -O -J -L http://cmake.org/files/v3.8/cmake-3.8.2-Linux-x86_64.tar.gz; \
		apt-get -y purge cmake; \
		tar -xvf ./cmake-3.8.2-Linux-x86_64.tar.gz; \
		cd ./cmake-3.8.2-Linux-x86_64; \
		cp -r bin /usr/; \
		cp -r doc /usr/share/; \
		cp -r man /usr/share/; \
		cp -r share /usr/;
