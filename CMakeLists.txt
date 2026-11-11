cmake_minimum_required(VERSION 3.20)
project(esm_deinterleaver CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # clang-tidy reads this later

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# I add each module's .cpp here as I write it
add_library(esm
  src/pdw_io.cpp
  src/clusterer.cpp
)
target_include_directories(esm PUBLIC include)
target_compile_options(esm PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Werror)

add_executable(esm_cli apps/esm_cli.cpp)
target_link_libraries(esm_cli PRIVATE esm)
target_compile_options(esm_cli PRIVATE -Wall -Wextra -Wpedantic -Werror)

enable_testing()
find_package(GTest REQUIRED)
add_executable(esm_tests
  tests/test_pdw_io.cpp
  tests/test_clusterer.cpp
)
target_link_libraries(esm_tests PRIVATE esm GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(esm_tests)
