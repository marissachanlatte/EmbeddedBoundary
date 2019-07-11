# Embedded Boundary
#
# CMakeLists derived from BART example file by J.S. Rehak
#

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(EmbeddedBoundary)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

### DEPENDENCIES #####################################################

### GTEST/GMOCK  ######################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(GTEST_ROOT "/usr")
find_package(GTest REQUIRED 1.8.0)

### EXECUTABLES AND SOURCES ##########################################

# Recursively find source files

file(GLOB_RECURSE sources "src/[a-zA-Z]*.cpp" "src/[a-zA-Z]*.cc")
list(APPEND sources ${PROTO_SRCS} ${PROTO_HDRS})
set(testing_sources ${sources})
list(FILTER sources EXCLUDE REGEX ".*/tests/.*")
list(FILTER sources EXCLUDE REGEX ".*/test_helpers/.*")
list(FILTER sources EXCLUDE REGEX ".*/test_main.cc$")
list(FILTER testing_sources EXCLUDE REGEX ".*/main.cc$")

# Include directories
include_directories(${GTEST_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/inc
  ${CMAKE_SOURCE_DIR}/src
  ${PROTO_INC_DIR})

# Add executables
ADD_EXECUTABLE(makeGeometry ${sources})
ADD_EXECUTABLE(geometry_test ${testing_sources})

# Add testing definition and library to geometry_test
target_compile_definitions(geometry_test PUBLIC -DTEST)

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

# Creates copies of .output files in the BART/test directory for CTEST

add_custom_target(copy_ctest_output_files ALL)

file(GLOB_RECURSE ctest_files "tests/*.output")
foreach(test_file ${ctest_files})
  string(REGEX REPLACE ".*BART/tests" "${CMAKE_CURRENT_BINARY_DIR}/tests"
    file_name ${test_file})
  add_custom_command(
    TARGET copy_ctest_output_files PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${test_file} ${file_name})
endforeach()
add_dependencies(geometry_test copy_ctest_output_files)
