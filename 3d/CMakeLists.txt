#CMake project file
#author: P. Valenta (valenpe7@email.cz), 2017

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/3d/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/3d/lib)

set(SOURCES
  ${CMAKE_SOURCE_DIR}/3d/src/main.cpp
  ${CMAKE_SOURCE_DIR}/3d/src/domain_param.cpp
  ${CMAKE_SOURCE_DIR}/3d/src/laser_param.cpp
  ${CMAKE_SOURCE_DIR}/3d/src/global.cpp
  ${CMAKE_SOURCE_DIR}/3d/src/laser_bcs.cpp
)

add_library(lbcs3d ${SOURCES})
target_link_libraries(lbcs3d LINK_PRIVATE ${MPI_CXX_LIBRARIES} ${FFTW_LIBRARIES})
set_target_properties(lbcs3d PROPERTIES LINKER_LANGUAGE CXX)
