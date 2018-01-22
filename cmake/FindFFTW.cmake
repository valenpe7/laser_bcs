# Find the native FFTW includes and library
#
# FFTW_INCLUDES    - where to find fftw3.h
# FFTW_LIBRARIES   - List of libraries when using FFTW.
# FFTW_FOUND       - True if FFTW found.

if(FFTW_INCLUDES)
  # Already in cache, be silent
  set(FFTW_FIND_QUIETLY TRUE)
endif(FFTW_INCLUDES)

if(UNIX)
	find_path(FFTW_INCLUDES NAMES fftw3.h HINTS $ENV{FFTW_INC})
	find_library(FFTW_LIBRARIES NAMES libfftw3.a HINTS $ENV{FFTW_LIB})
endif()

if(WIN32)
	find_path(FFTW_INCLUDES NAMES fftw3.h HINTS $ENV{FFTW_INC})
	find_library(FFTW_LIBRARIES NAMES libfftw3-3.lib HINTS $ENV{FFTW_LIB})
endif()

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced(FFTW_LIBRARIES FFTW_INCLUDES)
