cmake_minimum_required(VERSION 3.0.2)
project(lkh_tsp_solver)

set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_FLAGS "-O3 -Wall -std=c++17") # For C++ compilation, skip C files

# add_compile_options(-std=c++17) # For C++ Compiler Standards
find_package(Eigen3 REQUIRED)
find_package(catkin REQUIRED COMPONENTS roscpp random_numbers nav_msgs geometry_msgs)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS roscpp
  LIBRARIES lkh_tsp_solver
)

include_directories(SYSTEM include ${catkin_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

# Build library
file(GLOB SRCS src/*.c) # message("${SRCS}")
add_library(lkh_tsp_solver ${SRCS} src/lkh_interface.cpp)
target_link_libraries(lkh_tsp_solver -lm)

# Build main program executable
add_executable(lkh_tsp ${SRCS})
target_link_libraries(lkh_tsp ${catkin_LIBRARIES} -lm)

# Build library use example
add_executable(lkh_tsp_node src/lkh_tsp_node.cpp)
target_link_libraries(lkh_tsp_node ${catkin_LIBRARIES} lkh_tsp_solver -lm)
