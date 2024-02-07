#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/siesta.jcf
                ${CMAKE_CURRENT_BINARY_DIR}/siesta.jcf
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test.vmec
                COPYONLY)

#  Run VMEC and SIESTA the compare the results.
add_test (NAME    siesta_init_vmec_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test.vmec)
set_tests_properties (siesta_init_vmec_test
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)
add_test (NAME    siesta_init_siesta_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xsiesta,BINARY_DIR>/xsiesta)
set_tests_properties (siesta_init_siesta_test
                      PROPERTIES
                      DEPENDS siesta_init_vmec_test
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

add_test (NAME    siesta_init_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=test -min=0.0 -max=1.0 -tol=5.0E-5)
set_tests_properties (siesta_init_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_curtor_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=curtor -tol=0.18)
set_tests_properties (siesta_init_curtor_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_chipf_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=chipf -min=0.0 -max=1.0 -tol=1.2E-6)
set_tests_properties (siesta_init_chipf_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_phipf_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=phipf -min=0.0 -max=1.0 -tol=8.7E-18)
set_tests_properties (siesta_init_phipf_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_r_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.0 -max=0.01 -tol=5.1E-3)
set_tests_properties (siesta_init_r_1_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_r_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.01 -max=0.04 -tol=3.7E-4)
set_tests_properties (siesta_init_r_2_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_r_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.04 -max=0.2 -tol=6.7E-5)
set_tests_properties (siesta_init_r_3_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_r_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.2 -max=1.0 -tol=6.9E-6)
set_tests_properties (siesta_init_r_4_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_drdu_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.0 -max=0.01 -tol=5.1E-3)
set_tests_properties (siesta_init_drdu_1_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_drdu_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.01 -max=0.04 -tol=3.6E-4)
set_tests_properties (siesta_init_drdu_2_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_drdu_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.04 -max=0.2 -tol=6.7E-5)
set_tests_properties (siesta_init_drdu_3_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_drdu_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.2 -max=1.0 -tol=6.9E-6)
set_tests_properties (siesta_init_drdu_4_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_drdv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdv -min=0.0 -max=1.0 -tol=1.0E-20)
set_tests_properties (siesta_init_drdv_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_z_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.0 -max=0.01 -tol=5.2E-3)
set_tests_properties (siesta_init_z_1_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_z_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.01 -max=0.04 -tol=3.7E-4)
set_tests_properties (siesta_init_z_2_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_z_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.04 -max=0.2 -tol=6.8E-5)
set_tests_properties (siesta_init_z_3_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_z_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.2 -max=1.0 -tol=7.0E-6)
set_tests_properties (siesta_init_z_4_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_dzdu_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.0 -max=0.01 -tol=5.2E-3)
set_tests_properties (siesta_init_dzdu_1_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_dzdu_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.01 -max=0.04 -tol=3.7E-4)
set_tests_properties (siesta_init_dzdu_2_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_dzdu_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.04 -max=0.2 -tol=6.8E-5)
set_tests_properties (siesta_init_dzdu_3_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_dzdu_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.2 -max=1.0 -tol=7.0E-6)
set_tests_properties (siesta_init_dzdu_4_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_dzdv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdv -min=0.0 -max=1.0 -tol=1.0E-20)
set_tests_properties (siesta_init_dzdv_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_jacobian_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jacobian -min=0.0 -max=1.0 -tol=6.0E-5)
set_tests_properties (siesta_init_jacobian_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_jbsups_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=1.0E-20)
set_tests_properties (siesta_init_jbsups_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_jbsupu_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=4.9E-5)
set_tests_properties (siesta_init_jbsupu_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_jbsupv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupv -min=0.0 -max=1.0 -tol=1.8E-4)
set_tests_properties (siesta_init_jbsupv_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_vmec_divb_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=vmec_divb -min=0.0 -max=1.0 -tol=3.4E-7)
set_tests_properties (siesta_init_vmec_divb_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_siesta_divb_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=siesta_divb -min=0.0 -max=1.0 -tol=1.0E-20)
set_tests_properties (siesta_init_siesta_divb_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_siesta_pressure_test
    COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=pressure -min=0.0 -max=1.0 -tol=5.3E-1)
set_tests_properties (siesta_init_siesta_pressure_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_bsups_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsups -min=0.0 -max=1.0 -tol=1.0E-20)
set_tests_properties (siesta_init_bsups_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_bsupu_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsupu -min=0.005 -max=0.99 -tol=3.5E-3 -relative)
set_tests_properties (siesta_init_bsupu_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_bsupv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsupv -min=0.015 -max=0.99 -tol=5.5E-3 -relative)
set_tests_properties (siesta_init_bsupv_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_bsubu_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsubu -min=0.005 -max=1.0 -tol=5.5E-3 -relative)
set_tests_properties (siesta_init_bsubu_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_bsubv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsubv -min=0.005 -max=1.0 -tol=5.2E-6 -relative)
set_tests_properties (siesta_init_bsubv_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")

add_test (NAME    siesta_init_jksups_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=1.0E-20)
set_tests_properties (siesta_init_jksups_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_jksupu_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=0.99 -relative)
set_tests_properties (siesta_init_jksupu_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
add_test (NAME    siesta_init_jksupv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupv -min=0.0 -max=1.0 -tol=0.95 -relative)
set_tests_properties (siesta_init_jksupv_test
                      PROPERTIES
                      DEPENDS "siesta_init_vmec_test;siesta_init_siesta_test")
