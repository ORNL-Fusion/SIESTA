#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/siesta.jcf
                ${CMAKE_CURRENT_BINARY_DIR}/siesta.jcf
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test.vmec
                COPYONLY)

#  Run VMEC and SIESTA the compare the results.
add_test (NAME    siesta_sparse_vmec_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test.vmec)
set_tests_properties (siesta_sparse_vmec_test
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)
add_test (NAME    siesta_sparse_siesta_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xsiesta,BINARY_DIR>/xsiesta)
set_tests_properties (siesta_sparse_siesta_test
                      PROPERTIES
                      DEPENDS siesta_sparse_vmec_test
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

add_test (NAME    siesta_sparse_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=test -min=0.0 -max=1.0 -tol=5.0E-5)
set_tests_properties (siesta_sparse_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_curtor_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=curtor -tol=13.0)
set_tests_properties (siesta_sparse_curtor_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_chipf_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=chipf -min=0.0 -max=1.0 -tol=1.2E-6)
set_tests_properties (siesta_sparse_chipf_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_phipf_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=phipf -min=0.0 -max=1.0 -tol=1.1E-17)
set_tests_properties (siesta_sparse_phipf_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_r_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.0 -max=0.01 -tol=5.4E-3)
set_tests_properties (siesta_sparse_r_1_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_r_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.01 -max=0.04 -tol=3.8E-4)
set_tests_properties (siesta_sparse_r_2_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_r_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.04 -max=0.2 -tol=7.1E-5)
set_tests_properties (siesta_sparse_r_3_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_r_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.2 -max=1.0 -tol=7.2E-6)
set_tests_properties (siesta_sparse_r_4_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_drdu_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.0 -max=0.01 -tol=5.4E-3)
set_tests_properties (siesta_sparse_drdu_1_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_drdu_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.01 -max=0.04 -tol=3.8E-4)
set_tests_properties (siesta_sparse_drdu_2_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_drdu_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.04 -max=0.2 -tol=7.1E-5)
set_tests_properties (siesta_sparse_drdu_3_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_drdu_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.2 -max=1.0 -tol=7.2E-6)
set_tests_properties (siesta_sparse_drdu_4_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_drdv_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdv -min=0.0 -max=0.1 -tol=5.6E-4)
set_tests_properties (siesta_sparse_drdv_1_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_drdv_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdv -min=0.1 -max=0.4 -tol=2.1E-6)
set_tests_properties (siesta_sparse_drdv_2_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_drdv_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdv -min=0.4 -max=1.0 -tol=7.3E-7)
set_tests_properties (siesta_sparse_drdv_3_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_z_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.0 -max=0.01 -tol=5.2E-3)
set_tests_properties (siesta_sparse_z_1_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_z_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.01 -max=0.04 -tol=3.7E-4)
set_tests_properties (siesta_sparse_z_2_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_z_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.04 -max=0.2 -tol=6.8E-5)
set_tests_properties (siesta_sparse_z_3_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_z_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.2 -max=1.0 -tol=7.6E-6)
set_tests_properties (siesta_sparse_z_4_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_dzdu_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.0 -max=0.01 -tol=5.2E-3)
set_tests_properties (siesta_sparse_dzdu_1_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_dzdu_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.01 -max=0.04 -tol=3.7E-4)
set_tests_properties (siesta_sparse_dzdu_2_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_dzdu_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.04 -max=0.2 -tol=6.8E-5)
set_tests_properties (siesta_sparse_dzdu_3_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_dzdu_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.2 -max=1.0 -tol=7.1E-6)
set_tests_properties (siesta_sparse_dzdu_4_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_dzdv_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdv -min=0.0 -max=0.1 -tol=2.6E-5)
set_tests_properties (siesta_sparse_dzdv_1_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_dzdv_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdv -min=0.1 -max=0.4 -tol=1.5E-7)
set_tests_properties (siesta_sparse_dzdv_2_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_dzdv_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdv -min=0.4 -max=1.0 -tol=3.6E-7)
set_tests_properties (siesta_sparse_dzdv_3_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_jacobian_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jacobian -min=0.0 -max=1.0 -tol=6.0E-5)
set_tests_properties (siesta_sparse_jacobian_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_jbsups_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=7.8E-6)
set_tests_properties (siesta_sparse_jbsups_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_jbsupu_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=5.3E-5)
set_tests_properties (siesta_sparse_jbsupu_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_jbsupv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupv -min=0.0 -max=1.0 -tol=1.7E-4)
set_tests_properties (siesta_sparse_jbsupv_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_siesta_pressure_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=pressure -min=0.0 -max=1.0 -tol=0.77)
set_tests_properties (siesta_sparse_siesta_pressure_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")

add_test (NAME    siesta_sparse_jksups_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=7.8E-6)
set_tests_properties (siesta_sparse_jksups_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_jksupu_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=0.99 -relative)
set_tests_properties (siesta_sparse_jksupu_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
add_test (NAME    siesta_sparse_jksupv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupv -min=0.0 -max=1.0 -tol=0.97 -relative)
set_tests_properties (siesta_sparse_jksupv_test
                      PROPERTIES
                      DEPENDS "siesta_sparse_vmec_test;siesta_sparse_siesta_test")
