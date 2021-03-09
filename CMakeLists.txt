cmake_minimum_required (VERSION 3.14)

project (siesta C CXX Fortran)

add_library (siesta)
add_dependencies (siesta stell)

target_link_libraries (siesta stell)
target_include_directories (siesta PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_executable (xsiesta)
add_dependencies (xsiesta siesta)
target_link_libraries (xsiesta siesta)

add_subdirectory (Sources)

