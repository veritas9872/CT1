cmake_minimum_required(VERSION 3.14)
project(CT1)

set(CMAKE_CXX_STANDARD 11)

add_executable(CT1 forward_projection.cpp)

find_package(Matlab REQUIRED)

target_include_directories(CT1 PUBLIC ${Matlab_INCLUDE_DIRS})

# This must come after add_executable.
matlab_add_mex(NAME forward_projection SRC forward_projection.cpp R2018a)  # Assuming MATLAB version is 2018a or higher.