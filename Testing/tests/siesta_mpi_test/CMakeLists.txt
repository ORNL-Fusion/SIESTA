#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/siesta.jcf
                ${CMAKE_CURRENT_BINARY_DIR}/siesta.jcf
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test.vmec
                COPYONLY)

#  Run VMEC and SIESTA the compare the results.
add_test (NAME    siesta_mpi_vmec_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test.vmec)
set_tests_properties (siesta_mpi_vmec_test
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

macro (siesta_mpi_tests nproc)
    math (EXPR last_nproc "${nproc} - 1")

    add_test (NAME    siesta_mpi_siesta_${nproc}_test
              COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> ${nproc} $<TARGET_PROPERTY:xsiesta,BINARY_DIR>/xsiesta)
    set_tests_properties (siesta_mpi_siesta_${nproc}_test
                          PROPERTIES
                          DEPENDS $<IF:$<BOOL:${last_nproc}>,"siesta_mpi_vmec_test;siesta_mpi_${last_nproc}_test;siesta_mpi_${last_nproc}_curtor_test;siesta_mpi_${last_nproc}_chipf_test;siesta_mpi_${last_nproc}_phipf_test;siesta_mpi_${last_nproc}_r_1_test;siesta_mpi_${last_nproc}_r_2_test;siesta_mpi_${last_nproc}_r_3_test;siesta_mpi_${last_nproc}_r_4_test;siesta_mpi_${last_nproc}_drdu_1_test;siesta_mpi_${last_nproc}_drdu_2_test;siesta_mpi_${last_nproc}_drdu_3_test;siesta_mpi_${last_nproc}_drdu_4_test;siesta_mpi_${last_nproc}_drdv_test;siesta_mpi_${last_nproc}_z_1_test;siesta_mpi_${last_nproc}_z_2_test;siesta_mpi_${last_nproc}_z_3_test;siesta_mpi_${last_nproc}_z_4_test;siesta_mpi_${last_nproc}_dzdu_1_test;siesta_mpi_${last_nproc}_dzdu_2_test;siesta_mpi_${last_nproc}_dzdu_3_test;siesta_mpi_${last_nproc}_dzdu_4_test;siesta_mpi_${last_nproc}_dzdv_test;siesta_mpi_${last_nproc}_jacobian_test;siesta_mpi_${last_nproc}_jbsups_test;siesta_mpi_${last_nproc}_jbsupu_test;siesta_mpi_${last_nproc}_jbsupv_test;siesta_mpi_${last_nproc}_siesta_pressure_test;siesta_mpi_${last_nproc}_jksups_test;siesta_mpi_${last_nproc}_jksupu_test;siesta_mpi_${last_nproc}_jksupv_test","siesta_mpi_vmec_test">
                          PROCESSORS ${nproc})

    add_test (NAME    siesta_mpi_${nproc}_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=test -min=0.0 -max=1.0 -tol=5.0E-5)
    set_tests_properties (siesta_mpi_${nproc}_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_curtor_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=curtor -tol=0.69)
    set_tests_properties (siesta_mpi_${nproc}_curtor_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_chipf_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=chipf -min=0.0 -max=1.0 -tol=1.2E-6)
    set_tests_properties (siesta_mpi_${nproc}_chipf_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_phipf_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=phipf -min=0.0 -max=1.0 -tol=1.1E-17)
    set_tests_properties (siesta_mpi_${nproc}_phipf_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_r_1_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.0 -max=0.01 -tol=5.1E-3)
    set_tests_properties (siesta_mpi_${nproc}_r_1_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_r_2_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.01 -max=0.04 -tol=3.7E-4)
    set_tests_properties (siesta_mpi_${nproc}_r_2_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_r_3_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.04 -max=0.2 -tol=6.7E-5)
    set_tests_properties (siesta_mpi_${nproc}_r_3_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_r_4_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.2 -max=1.0 -tol=6.9E-6)
    set_tests_properties (siesta_mpi_${nproc}_r_4_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_drdu_1_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.0 -max=0.01 -tol=5.1E-3)
    set_tests_properties (siesta_mpi_${nproc}_drdu_1_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_drdu_2_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.01 -max=0.04 -tol=3.6E-4)
    set_tests_properties (siesta_mpi_${nproc}_drdu_2_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_drdu_3_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.04 -max=0.2 -tol=6.7E-5)
    set_tests_properties (siesta_mpi_${nproc}_drdu_3_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_drdu_4_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.2 -max=1.0 -tol=6.9E-6)
    set_tests_properties (siesta_mpi_${nproc}_drdu_4_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_drdv_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdv -min=0.0 -max=1.0 -tol=1.0E-20)
    set_tests_properties (siesta_mpi_${nproc}_drdv_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_z_1_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.0 -max=0.01 -tol=5.2E-3)
    set_tests_properties (siesta_mpi_${nproc}_z_1_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_z_2_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.01 -max=0.04 -tol=3.7E-4)
    set_tests_properties (siesta_mpi_${nproc}_z_2_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_z_3_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.04 -max=0.2 -tol=6.8E-5)
    set_tests_properties (siesta_mpi_${nproc}_z_3_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_z_4_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.2 -max=1.0 -tol=7.0E-6)
    set_tests_properties (siesta_mpi_${nproc}_z_4_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_dzdu_1_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.0 -max=0.01 -tol=5.2E-3)
    set_tests_properties (siesta_mpi_${nproc}_dzdu_1_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_dzdu_2_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.01 -max=0.04 -tol=3.7E-4)
    set_tests_properties (siesta_mpi_${nproc}_dzdu_2_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_dzdu_3_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.04 -max=0.2 -tol=6.8E-5)
    set_tests_properties (siesta_mpi_${nproc}_dzdu_3_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_dzdu_4_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.2 -max=1.0 -tol=7.0E-6)
    set_tests_properties (siesta_mpi_${nproc}_dzdu_4_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_dzdv_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdv -min=0.0 -max=1.0 -tol=1.0E-20)
    set_tests_properties (siesta_mpi_${nproc}_dzdv_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_jacobian_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jacobian -min=0.0 -max=1.0 -tol=6.0E-5)
    set_tests_properties (siesta_mpi_${nproc}_jacobian_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_jbsups_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=1.2E-6)
    set_tests_properties (siesta_mpi_${nproc}_jbsups_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_jbsupu_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=5.2E-5)
    set_tests_properties (siesta_mpi_${nproc}_jbsupu_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_jbsupv_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupv -min=0.0 -max=1.0 -tol=1.8E-4)
    set_tests_properties (siesta_mpi_${nproc}_jbsupv_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_siesta_pressure_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=pressure -min=0.0 -max=1.0 -tol=0.57)
    set_tests_properties (siesta_mpi_${nproc}_siesta_pressure_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")

    add_test (NAME    siesta_mpi_${nproc}_jksups_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=1.2E-6)
    set_tests_properties (siesta_mpi_${nproc}_jksups_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_jksupu_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=0.99 -relative)
    set_tests_properties (siesta_mpi_${nproc}_jksupu_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
    add_test (NAME    siesta_mpi_${nproc}_jksupv_test
              COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupv -min=0.0 -max=1.0 -tol=0.99 -relative)
    set_tests_properties (siesta_mpi_${nproc}_jksupv_test
                          PROPERTIES
                          DEPENDS "siesta_mpi_vmec_test;siesta_mpi_siesta_${nproc}_test")
endmacro ()

foreach (i RANGE 1 ${MPIEXEC_MAX_NUMPROCS})
    siesta_mpi_tests(${i})
endforeach ()
