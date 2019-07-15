FROM mrfil/powergrid-dev:latest

# Run everything in /root
WORKDIR /root

# Looks like we lost the PGI environmental variable
ENV PGI=/opt/pgi
RUN echo $PGI; 

RUN git clone https://github.com/mrfil/PowerGrid.git; \
    cd PowerGrid; \
    git fetch; \
    git checkout master; \
    mkdir build; \
    cd build; \
    cmake ../ -DCMAKE_CXX_COMPILER=pgc++ -DCMAKE_INSTALL_PREFIX=/opt/PowerGrid; \
    make -j4; make install

# For singularity compatibility
RUN ldconfig





