cmake_minimum_required(VERSION 2.8)
add_definitions(-std=c++11)
project( FiberMpp )
find_package( OpenCV REQUIRED )
include_directories(${OpenCV_INCLUDE_DIRS})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(FiberMppLibs STATIC
								 CylinderModel3D.cpp
								 model3D_types.cpp
								 randlib.cpp
								 )

add_executable( FiberMpp main.cpp)
target_link_libraries(FiberMpp FiberMppLibs)
target_link_libraries( FiberMpp ${OpenCV_LIBS} )
