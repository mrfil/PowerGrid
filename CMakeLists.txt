cmake_minimum_required(VERSION 3.2) #CMake version 3.2 required to properly support CMAKE_CXX_STANDARD

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED 1)
set(CMAKE_CXX_EXTENSIONS OFF)

project("PowerGrid" CXX)

option(OPENACC "Enable OpenACC support" ON)
option(ISMRMRD "Enable ISMRMRD Suppport" ON)
option(MPISupport "Enable MPI Support" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Find Matio
find_package(MATIO REQUIRED)
include_directories(${MATIO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MATIO_LIBRARIES})

# Find Boost
set(BOOST_LIBRARYDIR /usr/lib/x86_64-linux-gnu)

if(MPISupport)
    find_package(MPI)
    set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
    include_directories(${MPI_CXX_INCLUDE_PATH})
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization mpi)
else(MPISupport)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization)
endif(MPISupport)

set(Boost_USE_STATIC_LIBS OFF)
set(LIBS ${LIBS} ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})


if(OPENACC)
	find_path(cuFFT_INCLUDE_DIR cufft.h HINTS /usr/local/cuda/include
$ENV{CUDATOOLKIT_HOME}/include )
	include_directories(${cuFFT_INCLUDE_DIR})
	find_library(cuFFT_LIBRARIES cufft HINTS /usr/local/cuda/lib64
	$ENV{CUDATOOLKIT_HOME}/lib64)
	set(LIBS ${LIBS} ${cuFFT_LIBRARIES})
	find_library(cudaRT_LIBRARIES cudart HINTS /usr/local/cuda/lib64)
	set(LIBS ${LIBS} ${cudaRT_LIBRARIES})
	find_path(FFTW_INCLUDE_DIR fftw3.h)
				find_library(FFTW_LIBRARIES NAMES fftw3 HINTS /usr/lib/x86_64-linux-gnu/)
				find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
else()
	#find_package(FFTW REQUIRED)
	find_path(FFTW_INCLUDE_DIR fftw3.h)
        find_library(FFTW_LIBRARIES NAMES fftw3 HINTS /usr/lib/x86_64-linux-gnu/)
        find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
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
    find_package(HDF5 COMPONENTS CXX)
    include_directories(${HDF5_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
		set(LIBS ${LIBS} ${HDF5_CXX_LIBRARIES})
endif()

include_directories(Support/ArmaExtensions)

#Define files to be compiled
set(PG_FILES
	${PROJECT_SOURCE_DIR}/PowerGrid/Robject.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/QuadPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/TVPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/fftCPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/fftGPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/ftCpu.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gdft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gfft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/TimeSegmentation.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/griddingSupport.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/gridding.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gnufft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/SENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/pcSENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/reconSolve.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/GdftRecon.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/GnufftRecon.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/DWIRecon.cpp

)

set(POWERGRID_TEST_GDFT
    ./PowerGrid/PowerGridGdft.cpp
    ./PowerGrid/config.cxx
)

set(POWERGRID_TEST_Gnufft
	./PowerGrid/PowerGridGnufft.cpp
	./PowerGrid/config.cxx
)

add_executable(PowerGridGdft ${PG_FILES} ${POWERGRID_TEST_GDFT})
add_executable(PowerGridGnufft ${PG_FILES} ${POWERGRID_TEST_Gnufft})

target_link_libraries(PowerGridGdft ${LIBS})
target_link_libraries(PowerGridGnufft ${LIBS})

if(MPISupport)
 include_directories(MPI)
endif(MPISupport)

if (ISMRMRD)
	set(POWERGRID_ISMRMRD
			${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.cpp
			${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.h
			./PowerGrid/PowerGridIsmrmrd.cpp
			./PowerGrid/config.cxx
	)
	find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)
	set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})
	add_executable(PowerGridIsmrmrd ${POWERGRID_ISMRMRD} ${PG_FILES})
	target_link_libraries(PowerGridIsmrmrd ${LIBS})

	set(POWERGRID_PCSENSE
			${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.cpp
			${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.h
			./PowerGrid/PowerGridPcSense.cpp
			)
	add_executable(PowerGridPcSense ${POWERGRID_PCSENSE} ${PG_FILES})
	target_link_libraries(PowerGridPcSense ${LIBS})
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include")

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -flto")
else()
	if(CMAKE_COMPILER_IS_GNUCXX)

	else()
	add_definitions(-DBOOST_NO_ALIGNMENT)
		if(OPENACC)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -gopt -fast -O4 -Minfo=accel,ccff -ta=tesla:cc60,cuda9.0,fastmath -pgf77libs -fpermissive")
			#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -acc -g -Minfo=accel,ccff -ta=multicore -pgf77libs -fpermissive")
			#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Minfo=ccff -fpermissive -pgf77libs")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --c++11")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -Minfo=ccff -fpermissive -pgf77libs")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --c++11")

        endif()
    endif()
endif()


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --c++11")

#add_definitions(-D__STRICT_ANSI__)
add_definitions(-DARMA_NO_DEBUG)
#add_definitions(-DARMA_DONT_USE_CXX11)
#add_definitions(-DARMA_USE_OPENMP)

if (MPISupport)
    add_subdirectory(PowerGrid/MPI)
endif (MPISupport)
