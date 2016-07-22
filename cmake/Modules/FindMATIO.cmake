# Find the MATIO headers and library.
#
# MATIO_INCLUDE_DIRS - where to find matio.h, etc.
# MATIO_LIBRARIES - List of libraries.
# MATIO_FOUND - True if matio found.
# Look for the header file.
find_path(MATIO_INCLUDE_DIR NAMES matio.h
		PATHS /usr/include )
mark_as_advanced(MATIO_INCLUDE_DIR)
# Look for the library.
find_library(MATIO_LIBRARY NAMES matio
			   PATH_SUFFIXES lib64 libs lib
		           PATHS  /usr/local /usr /usr/lib /usr/lib/x86_64-linux-gnu/ )
mark_as_advanced(MATIO_LIBRARY)
# handle the QUIETLY and REQUIRED arguments and set MATIO_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MATIO DEFAULT_MSG MATIO_LIBRARY MATIO_INCLUDE_DIR)
if(MATIO_FOUND)
set(MATIO_LIBRARIES ${MATIO_LIBRARY} ${HDF5_LIBRARIES})
set(MATIO_INCLUDE_DIRS ${MATIO_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})
else()
set(MATIO_LIBRARIES)
set(MATIO_INCLUDE_DIRS)
endif()
