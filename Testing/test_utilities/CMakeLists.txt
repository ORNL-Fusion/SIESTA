target_sources(xsiesta_test
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vmec_quantities.hpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/siesta_test.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/siesta_quantities.hpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/grid_quantity.hpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/parity.hpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/grid_quantity.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/commandline_parser.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/commandline_parser.hpp>
)

target_sources(xsiesta_unit_tests
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/siesta_unit_tests.f>
)
