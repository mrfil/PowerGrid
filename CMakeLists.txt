cmake_minimum_required(VERSION 3.2) #CMake version 3.2 required to properly support CMAKE_CXX_STANDARD

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 11 REQUIRED)
#set(CMAKE_CXX_STANDARD_REQUIRED 1)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -flto")
endif()

project("PowerGrid" CXX)

option(OPENACC "Enable OpenACC support" OFF)
option(ISMRMRD "Enable ISMRMRD Suppport" OFF)
option(MPI "Enable MPI Support" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Find Matio
find_package(MATIO REQUIRED)
include_directories(${MATIO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MATIO_LIBRARIES})

# Find Boost
if(MPI)
find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization mpi)
else(MPI)
find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization)
set(LIBS ${LIBS} ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})
endif()

if(OPENACC)
	find_path(cuFFT_INCLUDE_DIR cufft.h HINTS /usr/local/cuda/include
$ENV{CUDATOOLKIT_HOME}/include )
	include_directories(${cuFFT_INCLUDE_DIR})
	find_library(cuFFT_LIBRARIES cufft HINTS /usr/local/cuda/lib64
	$ENV{CUDATOOLKIT_HOME}/lib64)
	set(LIBS ${LIBS} ${cuFFT_LIBRARIES})
	#find_library(cudaRT_LIBRARIES cudart HINTS /usr/local/cuda/lib64)
	#set(LIBS ${LIBS} ${cudaRT_LIBRARIES})
else()
	#find_package(FFTW REQUIRED)
	find_path(FFTW_INCLUDE_DIR fftw3.h)
        find_library(FFTW_LIBRARIES NAMES fftw3 HINTS /usr/lib/x86_64-linux-gnu/)
        find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
endif()

if(USE_ACML)
	find_package(ACML)
	set(LIBS ${LIBS} ${ACML_LIBRARIES})
endif()

# Find Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
message(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})

# Find Xerces
find_package(Xerces)
include_directories(${XERCESC_INCLUDE_DIRS})
set(LIBS ${LIBS} ${XERCESC_LIBRARY})

if(ISMRMRD)
    find_package(HDF5)
    include_directories(${HDF5_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
		set(LIBS ${LIBS} ${HDF5_CXX_LIBRARIES})
endif()



set(POWERGRID_TEST_GDFT
    ./PowerGrid/PowerGridGdft.cpp
    ./PowerGrid/config.cxx
)

set(POWERGRID_TEST_Gnufft
	./PowerGrid/PowerGridGnufft.cpp
	./PowerGrid/config.cxx
)

set(POWERGRID_TEST_DWI3DDFT
		./PowerGrid/PowerGridDWI3DDFT.cpp
		./PowerGrid/config.cxx
		)

set(POWERGRID_STEAM
		./PowerGrid/PowerGridSTEAM.cpp
		./PowerGrid/config.cxx
		)

add_executable(PowerGridGdft ${POWERGRID_TEST_GDFT})
add_executable(PowerGridGnufft ${POWERGRID_TEST_Gnufft})
add_executable(PowerGridDWI3DDFT ${POWERGRID_TEST_DWI3DDFT})
add_executable(PowerGridSTEAM ${POWERGRID_STEAM})
target_link_libraries(PowerGridGdft ${LIBS})
target_link_libraries(PowerGridGnufft ${LIBS})
target_link_libraries(PowerGridDWI3DDFT ${LIBS})
target_link_libraries(PowerGridSTEAM ${LIBS})

if (ISMRMRD)
	set(POWERGRID_ISMRMRD
			./PowerGrid/PowerGridIsmrmrd.cpp
			./PowerGrid/config.cxx
	)
	find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)
	set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})
	add_executable(PowerGridIsmrmrd ${POWERGRID_ISMRMRD})
	target_link_libraries(PowerGridIsmrmrd ${LIBS})
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include")
if(APPLE)
else()
	if(CMAKE_COMPILER_IS_GNUCXX)

	else()
		if(OPENACC)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -O4 -gopt -acc  -Minfo=accel,mp -Msmartalloc -ta=tesla:cc35 -pgf77libs -fpermissive")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gopt -fast -O4 -fpermissive -Msmartalloc -pgf77libs ")
		endif()
	endif()
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lfftw3")

#add_definitions(-D__STRICT_ANSI__)
add_definitions(-DARMA_NO_DEBUG)
#add_definitions(-DARMA_DONT_USE_CXX11)

if (MPI)
    add_subdirectory(PowerGrid/MPI)
endif ()
