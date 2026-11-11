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
  src/pri_analyzer.cpp
  src/pipeline.cpp
  src/emitter_library.cpp
  src/pulse_detector.cpp
)
target_include_directories(esm PUBLIC include)
target_compile_options(esm PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Werror)

add_executable(esm_cli apps/esm_cli.cpp)
target_link_libraries(esm_cli PRIVATE esm)
target_compile_options(esm_cli PRIVATE -Wall -Wextra -Wpedantic -Werror)

enable_testing()
# Use the system GoogleTest if it's installed (libgtest-dev). If not, download it.
# CI takes the download path so it doesn't depend on Ubuntu's apt mirrors.
find_package(GTest QUIET)
if(NOT GTest_FOUND)
  if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)  # stops a warning about extracted file timestamps
  endif()
  include(FetchContent)
  FetchContent_Declare(googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)  # I only use gtest
  FetchContent_MakeAvailable(googletest)
endif()
add_executable(esm_tests
  tests/test_pdw_io.cpp
  tests/test_clusterer.cpp
  tests/test_pri_analyzer.cpp
  tests/test_emitter_library.cpp
  tests/test_pulse_detector.cpp
)
target_link_libraries(esm_tests PRIVATE esm GTest::gtest_main)
# lets a test find config/emitter_library.csv no matter where ctest runs from
target_compile_definitions(esm_tests PRIVATE ESM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
include(GoogleTest)
gtest_discover_tests(esm_tests)
