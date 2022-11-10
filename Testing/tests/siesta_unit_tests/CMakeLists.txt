#  Run SIESTA unit tests.
add_test (NAME    siesta_unit_tests
          COMMAND $<TARGET_PROPERTY:xiesta_unit_tests,BINARY_DIR>/xiesta_unit_tests)
