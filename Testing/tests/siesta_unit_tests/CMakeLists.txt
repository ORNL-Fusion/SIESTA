#  Run SIESTA unit tests.
add_test (NAME    siesta_unit_tests
          COMMAND $<TARGET_PROPERTY:xsiesta_unit_tests,BINARY_DIR>/xsiesta_unit_tests)
