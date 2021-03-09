#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/siesta.jcf
                ${CMAKE_CURRENT_BINARY_DIR}/siesta.jcf
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test.vmec
                COPYONLY)

#  Run VMEC and SIESTA the compare the results.
add_test (NAME    siesta_final_vmec_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test.vmec)
add_test (NAME    siesta_final_siesta_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xsiesta,BINARY_DIR>/xsiesta)

add_test (NAME    siesta_final_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=test -min=0.0 -max=1.0 -tol=1.0E-3)

add_test (NAME    siesta_final_curtor_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=curtor -tol=10.0)

add_test (NAME    siesta_final_chipf_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=chipf -min=0.0 -max=1.0 -tol=2.0E-6)
add_test (NAME    siesta_final_phipf_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=phipf -min=0.0 -max=1.0 -tol=1.0E-17)

add_test (NAME    siesta_final_r_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.0 -max=0.01 -tol=1.0E-2)
add_test (NAME    siesta_final_r_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.01 -max=0.04 -tol=1.0E-3)
add_test (NAME    siesta_final_r_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.04 -max=0.2 -tol=1.0E-4)
add_test (NAME    siesta_final_r_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.2 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_final_drdu_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.0 -max=0.01 -tol=1.0E-2)
add_test (NAME    siesta_final_drdu_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.01 -max=0.04 -tol=1.0E-3)
add_test (NAME    siesta_final_drdu_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.04 -max=0.2 -tol=1.0E-4)
add_test (NAME    siesta_final_drdu_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.2 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_final_drdv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdv -min=0.0 -max=1.0 -tol=1.0E-20)

add_test (NAME    siesta_final_z_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.0 -max=0.01 -tol=1.0E-2)
add_test (NAME    siesta_final_z_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.01 -max=0.04 -tol=1.0E-3)
add_test (NAME    siesta_final_z_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.04 -max=0.2 -tol=1.0E-4)
add_test (NAME    siesta_final_z_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.2 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_final_dzdu_1_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.0 -max=0.01 -tol=1.0E-2)
add_test (NAME    siesta_final_dzdu_2_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.01 -max=0.04 -tol=1.0E-3)
add_test (NAME    siesta_final_dzdu_3_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.04 -max=0.2 -tol=1.0E-4)
add_test (NAME    siesta_final_dzdu_4_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.2 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_final_dzdv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdv -min=0.0 -max=1.0 -tol=1.0E-20)

add_test (NAME    siesta_final_jacobian_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jacobian -min=0.0 -max=1.0 -tol=1.0E-4)

add_test (NAME    siesta_final_jbsups_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=2.1E-4)
add_test (NAME    siesta_final_jbsupu_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=6.0E-4)
add_test (NAME    siesta_final_jbsupv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupv -min=0.0 -max=1.0 -tol=2.0E-4)

#add_test (NAME    siesta_final_vmec_divb_test
#          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=vmec_divb -min=0.0 -max=1.0 -tol=2.0E-5)
#add_test (NAME    siesta_final_siesta_divb_test
#          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=siesta_divb -min=0.0 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_final_siesta_pressure_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=pressure -min=0.0 -max=1.0 -tol=26.0)

#add_test (NAME    siesta_final_bsups_test
#          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsups -min=0.0 -max=1.0 -tol=0.004)
#add_test (NAME    siesta_final_bsupu_test
#          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsupu -min=0.005 -max=0.99 -tol=1.0E-3 -relative)
#add_test (NAME    siesta_final_bsupv_test
#          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsupv -min=0.015 -max=0.99 -tol=5.5E-4)

#add_test (NAME    siesta_final_bsubu_test
#          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsubu -min=0.005 -max=1.0 -tol=2.0E-1 -relative)
#add_test (NAME    siesta_final_bsubv_test
#          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsubv -min=0.005 -max=1.0 -tol=1.0E-5 -relative)

add_test (NAME    siesta_final_jksups_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=0.0002)
add_test (NAME    siesta_final_jksupu_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=2.0 -relative)
add_test (NAME    siesta_final_jksupv_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=2.0 -relative)
