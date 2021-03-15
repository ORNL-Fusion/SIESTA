cmake_minimum_required (VERSION 3.14)

project (siesta C CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

add_library (siesta)
add_dependencies (siesta bmw)

target_link_libraries (siesta bmw)
target_include_directories (siesta PUBLIC $<TARGET_PROPERTY:bmw,BINARY_DIR>)

add_executable (xsiesta)
add_dependencies (xsiesta siesta)
target_link_libraries (xsiesta siesta)

add_subdirectory (Sources)

################################################################################
#  Testing                                                                     #
################################################################################

#  Build test utilities.
add_executable (xsiesta_test)
add_dependencies (xsiesta_test stell)

target_link_libraries (xsiesta_test stell)

add_subdirectory (Testing)
