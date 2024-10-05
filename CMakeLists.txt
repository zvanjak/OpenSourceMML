cmake_minimum_required(VERSION 3.22)

set (CMAKE_CXX_STANDARD 20)

project(MinimalMathLibrary) 

add_subdirectory(src)

include(CTest)
enable_testing()
add_subdirectory(tests)