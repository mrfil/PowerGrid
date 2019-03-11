cmake_minimum_required(VERSION 3.2) #CMake version 3.2 required to properly support CMAKE_CXX_STANDARD

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_EXTENSIONS ON)

project("PowerGrid" C CXX)

option(OPENACC "Enable OpenACC support" ON)
option(MPISupport "Enable MPI Support" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Find Boost
set(BOOST_LIBRARYDIR /usr/lib/x86_64-linux-gnu)
#set(BOOST_ROOT /usr/local/lib)

if(MPISupport)
    find_package(MPI)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization mpi)
else(MPISupport)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization)
endif(MPISupport)

set(Boost_USE_STATIC_LIBS OFF)
set(LIBS ${LIBS} ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})


if(OPENACC)
	find_path(cuFFT_INCLUDE_DIR cufft.h HINTS $ENV{PGI}/linux86-64/2018/cuda/9.1/include /usr/local/pgi/linux86-64/2018/cuda/9.1/include /usr/local/cuda/include
		$ENV{CUDATOOLKIT_HOME}/include )
	include_directories(${cuFFT_INCLUDE_DIR})
	find_library(cuFFT_LIBRARIES cufft HINTS $ENV{PGI}/linux86-64/2018/cuda/9.1/lib64 /usr/local/pgi/linux86-64/2018/cuda/9.1/lib64 /usr/local/cuda/lib64
	$ENV{CUDATOOLKIT_HOME}/lib64)
	set(LIBS ${LIBS} ${cuFFT_LIBRARIES})
	find_library(cudaRT_LIBRARIES cudart HINTS $ENV{CRAY_CUDATOOLKIT_DIR}/lib64 $ENV{PGI}/linux86-64/2018/cuda/9.1/lib64 /usr/local/pgi/linux86-64/2018/cuda/9.1/lib64 /usr/local/cuda/lib64)
	set(LIBS ${LIBS} ${cudaRT_LIBRARIES})
	find_library(nvToolsExt_LIBRARIES nvToolsExt HINTS /software/cuda-toolkit-9.1/lib64 $ENV{CRAY_CUDATOOLKIT_DIR}/lib64 $ENV{PGI}/linux86-64/2018/cuda/9.1/lib64 /usr/local/pgi/linux86-64/2018/cuda/9.1/lib64 /usr/local/cuda/lib64)
	set(LIBS ${LIBS} ${nvToolsExt_LIBRARIES})
	find_path(FFTW_INCLUDE_DIR fftw3.h HINTS $ENV{FFTW_INC})
	find_library(FFTW_LIBRARIES NAMES fftw3 HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
	find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
else()
	#find_package(FFTW REQUIRED)
	find_path(FFTW_INCLUDE_DIR fftw3.h HINTS $ENV{FFTW_INC})
    find_library(FFTW_LIBRARIES NAMES fftw3 HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
    find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
endif()

# Find Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
message(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})
#set(HDF5_ROOT /usr/local)
#find_package(HDF5)
find_package(HDF5 COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_CXX_LIBRARIES})
#find_library(GPERFTOOLS_LIBRARIES NAMES profiler libprofiler.so.0 HINTS /usr/)
#set(LIBS ${LIBS} ${GPERFTOOLS_LIBRARIES})
include_directories(Support/ArmaExtensions)

#Define files to be compiled
set(PG_FILES
	${PROJECT_SOURCE_DIR}/PowerGrid/Robject.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/QuadPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/TVPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/fftCPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/fftGPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/ftCpu.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/ftCpuWithGrads.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gdft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/GdftR2.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gfft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/TimeSegmentation.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/griddingSupport.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/gridding.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gnufft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/SENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/pcSENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/pcSenseTimeSeg.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/reconSolve.cpp

)

if(MPISupport)
	include_directories(MPI)
endif(MPISupport)


set(POWERGRID_ISMRMRD
	${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.h
	./PowerGrid/PowerGridIsmrmrd.cpp
)
find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/ ~/CIVIC/lib)
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

set(POWERGRID_PCSENSE_TIMESEG
	${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.h
	./PowerGrid/PowerGridPcSenseTimeSeg.cpp
)
add_executable(PowerGridPcSenseTimeSeg ${POWERGRID_PCSENSE_TIMESEG} ${PG_FILES})
target_link_libraries(PowerGridPcSenseTimeSeg ${LIBS})

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -flto")
else()
	if(CMAKE_COMPILER_IS_GNUCXX)

	else()
	add_definitions(-DBOOST_NO_ALIGNMENT)
		if(OPENACC)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -mp -gopt -fast -O3 -tp=px -ta=tesla:cc35,cc60,cuda9.1,flushz -pgf77libs")
			#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -acc -g -Minfo=accel,ccff -ta=multicore -pgf77libs -fpermissive")
			#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Minfo=ccff -fpermissive -pgf77libs")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -Minfo=ccff -fpermissive ")
		endif()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --c++11")
    endif()
endif()

#add_definitions(-D__STRICT_ANSI__)
add_definitions(-DARMA_NO_DEBUG)
add_definitions(-DARMA_DONT_USE_HDF5)
add_definitions(-DARMA_DONT_USE_OPENMP)
#add_definitions(-DARMA_SIMPLE_LOOPS)
#add_definitions(-DUSE_NVTX)
if (MPISupport)
    add_subdirectory(PowerGrid/MPI)
endif (MPISupport)

install(TARGETS PowerGridPcSense DESTINATION bin)
install(TARGETS PowerGridIsmrmrd DESTINATION bin)