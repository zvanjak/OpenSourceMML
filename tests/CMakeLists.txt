include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.0)
FetchContent_MakeAvailable(Catch2)

file(GLOB MML_TESTS 
                    base/vector_tests.cpp
                )

add_executable("MML_Tests" ${MML_TESTS} )
target_link_libraries("MML_Tests" PRIVATE Catch2::Catch2WithMain)
target_include_directories(MML_Tests PRIVATE "../include")

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

include(CTest)
include(Catch)
catch_discover_tests("MML_Tests")
