#CMake project file
#author: P. Valenta (valenpe7@email.cz), 2017

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/2d/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/2d/lib)

set(SOURCES
  ${CMAKE_SOURCE_DIR}/2d/src/main.cpp
  ${CMAKE_SOURCE_DIR}/2d/src/domain_param.cpp
  ${CMAKE_SOURCE_DIR}/2d/src/laser_param.cpp
  ${CMAKE_SOURCE_DIR}/2d/src/global.cpp
  ${CMAKE_SOURCE_DIR}/2d/src/laser_bcs.cpp
)

add_library(lbcs2d ${SOURCES})
target_link_libraries(lbcs2d LINK_PRIVATE ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARIES})
set_target_properties(lbcs2d PROPERTIES LINKER_LANGUAGE CXX)
