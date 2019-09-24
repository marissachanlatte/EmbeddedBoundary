# Embedded Boundary
#
# CMakeLists derived from BART example file by J.S. Rehak
#

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
PROJECT(EmbeddedBoundary)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Use C++11
set(CMAKE_CXX_STANDARD 11)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

### DEPENDENCIES #####################################################
FIND_PACKAGE(Eigen3 3.3)

### GTEST/GMOCK  ######################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(GTEST_ROOT "/usr")
find_package(GTest REQUIRED 1.8.0)

### EXECUTABLES AND SOURCES ##########################################

# Recursively find source files

file(GLOB_RECURSE sources "src/[a-zA-Z]*.cpp" "src/[a-zA-Z]*.cc")
set(testing_sources ${sources})
list(FILTER sources EXCLUDE REGEX ".*/tests/.*")
list(FILTER sources EXCLUDE REGEX ".*/test_helpers/.*")
list(FILTER sources EXCLUDE REGEX ".*/test_main.cc$")
list(FILTER testing_sources EXCLUDE REGEX ".*/main.cc$")

# Include directories
include_directories(${GTEST_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/src)

# Add executables
ADD_EXECUTABLE(makeGeometry ${sources})
ADD_EXECUTABLE(geometry_test ${testing_sources})

set_target_properties(geometry_test
    PROPERTIES
    COMPILE_FLAGS "-ftest-coverage -fprofile-arcs"
    LINK_FLAGS    "--coverage"
)

# Add testing definition and library to geometry_test
target_compile_definitions(geometry_test PUBLIC -DTEST)

target_link_libraries(geometry_test ${GTEST_BOTH_LIBRARIES})

# Link Eigen to both executables
target_link_libraries(makeGeometry Eigen3::Eigen)
target_link_libraries(geometry_test Eigen3::Eigen)

### TEST FILES ##################################################
# Create copies of .gold files in the src directory for gtest

add_custom_target(copy_gtest_gold_files ALL)

add_custom_command(TARGET copy_gtest_gold_files PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory test_data)

file(GLOB_RECURSE gtest_files "src/*.gold")
foreach(test_file ${gtest_files})
  get_filename_component(file_name ${test_file} NAME)
  add_custom_command(
    TARGET copy_gtest_gold_files PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${test_file}
            ${CMAKE_CURRENT_BINARY_DIR}/test_data/${file_name})
endforeach()
add_dependencies(geometry_test copy_gtest_gold_files)

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE)
