# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 2.6)

add_definitions(-std=c++17)

project(FOAM)

include_directories(tetgen)

# Add an executable to the project using the specified source files.
add_executable(tetgen ./tetgen/tetgen.cxx ./tetgen/predicates.cxx)

#Add a library to the project using the specified source files. 
# In Linux/Unix, it will creates the libtet.a
add_library(tet STATIC ./tetgen/tetgen.cxx ./tetgen/predicates.cxx)

#Set properties on a target. 
#We use this here to set -DTETLIBRARY for when compiling the
#library
set_target_properties(tet PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)

add_executable(foam main.cpp CreateDxf.cpp)
target_link_libraries(foam tet)
