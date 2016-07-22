# File:   FindFFT.cmake
# Brief:  Find the FFT includes and library

# FFTW_INCLUDE_DIR = fftw3.h
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found

set(Libfftw fftw3)
if(QMC_BUILD_STATIC)
  set(Libfftw libfftw3.a)
endif()

if(FFTW_INCLUDE_DIRS)
  find_path(FFTW_INCLUDE_DIR fftw3.h  ${FFTW_INCLUDE_DIRS})
  find_library(FFTW_LIBRARY ${Libfftw} ${FFTW_LIBRARY_DIRS})
else()
  find_path(FFTW_INCLUDE_DIR fftw3.h ${FFTW_HOME}/include $ENV{FFTW_HOME}/include $ENV{FFTW_INC} )
  find_library(FFTW_LIBRARIES ${Libfftw} /usr/lib/x86_64-linux-gnu ${FFTW_HOME}/lib $ENV{FFTW_HOME}/lib $ENV{FFTW_DIR}) 
endif()

set(FFTW_FOUND FALSE)
if(FFTW_INCLUDE_DIR AND FFTW_LIBRARIES)
  message(STATUS "FFTW_INCLUDE_DIR=${FFTW_INCLUDE_DIR}")
  message(STATUS "FFTW_LIBRARIES=${FFTW_LIBRARIES}")
  set(FFTW_FOUND TRUE)
endif()

#mark_as_advanced(
#   FFTW_INCLUDE_DIR
#   FFTW_LIBRARIES
#   FFTW_FOUND
#)


