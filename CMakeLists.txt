cmake_minimum_required(VERSION 2.8)
project("PowerGrid" CXX)


option(OPENACC "Enable OpenACC support" OFF)
option(ISMRMRD "Enable ISMRMRD Suppport" OFF)
option(MPI "Enable MPI Support" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED 1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MATIO)
include_directories(${MATIO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MATIO_LIBRARIES})
set(BOOST_LIBRARYDIR /usr/lib/x86_64-linux-gnu/)
set(BOOST_ROOT /u/sciteam/cerjanic/PGProj/libraries)
find_package(Boost 1.43 COMPONENTS program_options serialization mpi)
set(LIBS ${LIBS} ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})


if(OPENACC)
	find_path(cuFFT_INCLUDE_DIR cufft.h HINTS /usr/local/cuda/include
$ENV{CUDATOOLKIT_HOME}/include )
	include_directories(${cuFFT_INCLUDE_DIR})
	find_library(cuFFT_LIBRARIES cufft HINTS /usr/local/cuda/lib64
	$ENV{CUDATOOLKIT_HOME}/lib64)
	set(LIBS ${LIBS} ${cuFFT_LIBRARIES})
	#find_library(cudaRT_LIBRARIES cudart HINTS /usr/local/cuda/lib64)
	#set(LIBS ${LIBS} ${cudaRT_LIBRARIES})
else(OPENACC)
	#find_package(FFTW REQUIRED)
	find_path(FFTW_INCLUDE_DIR fftw3.h)
find_library(FFTW_LIBRARIES NAMES fftw3 HINTS /usr/lib/x86_64-linux-gnu/)
find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
endif(OPENACC)
if(${USE_ACML})
	find_package(ACML)
	set(LIBS ${LIBS} ${ACML_LIBRARIES})
endif()


find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
message(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})
find_package(Xerces)
include_directories(${XERCESC_INCLUDE_DIRS})
set(LIBS ${LIBS} ${XERCESC_LIBRARY})
if(ISMRMRD)
    find_package(HDF5)
    include_directories(${HDF5_INCLUDE_DIRS})
endif(ISMRMRD)

set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_CXX_LIBRARIES})

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

endif (ISMRMRD)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include")
if(APPLE)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lgfortran")
else()
	if(CMAKE_COMPILER_IS_GNUCXX)
		#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ")
	else()
		if(OPENACC)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -gopt -Minfo=accel -ta=tesla:cc35 -fpermissive")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gopt -O -fpermissive -pgf77libs ")
		endif()
	endif()
endif()
#Needed for CMake 2.8 or so to not add the -rdynamic flag to compilers like PGI
#PGI chokes on link because it does not understand the -rdynamic flag
#Why CMake insists on passing this flag, I don't know...
#set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
#set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} -lfftw3")

#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX11)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#else()
#        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
#endif()
if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -flto")
endif(APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lfftw3")
#add_definitions(-D__STRICT_ANSI__)
#remove_definitions(-DARMA_NO_DEBUG)
add_definitions(-DARMA_DONT_USE_CXX11)


if (MPI)
    add_subdirectory(PowerGrid/MPI)
endif (MPI)
