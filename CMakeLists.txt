#CMake project file
#author: P. Valenta (valenpe7@email.cz), 2017

cmake_minimum_required(VERSION 3.3)
project(lbcs)
enable_language(CXX)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(BOOST_ROOT $ENV{BOOST_HOME})

find_package(MPI REQUIRED)
find_package(FFTW REQUIRED)
find_package(Boost REQUIRED)

include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${FFTW_INCLUDES})
include_directories(${Boost_INCLUDE_DIRS})

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release', Debug mode was not specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -no-prec-div -ansi-alias -qopt-prefetch=4 -unroll-aggressive -m64")
  set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -g -traceback -mp1 -fp-trap=common -fp-model strict")
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -msse4 -mtune=native -march=native -funroll-loops -fno-math-errno -ffast-math")
  set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -g -pedantic -Wall -Wextra -Wno-unused")
else()
  message(STATUS "No optimized compiler flags are known")
  message(STATUS "Compiler full path: " ${CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O2")
  set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0")
endif()

option(OPENMP "Enable multithreading using OpenMP directives." OFF)

if(OPENMP)
  message(STATUS "Option 'OPENMP' enabled")
  if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -openmp")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -openmp")
  elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp")
  else()
    message(STATUS "C++ OpenMP compiler flag not known.")
  endif()
endif()

add_subdirectory(2d)
add_subdirectory(3d)
