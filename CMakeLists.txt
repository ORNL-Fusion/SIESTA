cmake_minimum_required (VERSION 3.14)

project (siesta C CXX Fortran)

add_library (siesta)

target_link_libraries (siesta PUBLIC bmw)
target_include_directories (siesta

                            PUBLIC

                            $<TARGET_PROPERTY:siesta,BINARY_DIR>
)

add_executable (xsiesta)

target_link_libraries (xsiesta PUBLIC siesta)

add_subdirectory (Sources)

################################################################################
#  Testing                                                                     #
################################################################################

#  Build test utilities.
add_executable (xsiesta_test)

target_link_libraries (xsiesta_test PUBLIC stell)
target_compile_features (xsiesta_test PRIVATE cxx_std_11)

add_subdirectory (Testing)
