if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.16)
  project(amg_c_prec_extensions_test C CXX Fortran)
  find_package(amg4psblas CONFIG REQUIRED)
  find_package(MPI QUIET COMPONENTS C Fortran)
  enable_testing()
  set(amg_prec_test_library amg4psblas::amgcbind)
else()
  set(amg_prec_test_library amgcbind)
endif()

add_executable(amg_c_prec_extensions prec_extensions.c)
target_link_libraries(amg_c_prec_extensions PRIVATE ${amg_prec_test_library} m)
set_target_properties(amg_c_prec_extensions PROPERTIES LINKER_LANGUAGE Fortran)
target_compile_features(amg_c_prec_extensions PRIVATE c_std_99)

if(MPIEXEC_EXECUTABLE AND NOT PSB_SERIAL_MPI)
  foreach(nprocs IN ITEMS 1 2 4)
    add_test(NAME amg_c_prec_extensions_${nprocs}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nprocs}
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:amg_c_prec_extensions>
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties(amg_c_prec_extensions_${nprocs}
      PROPERTIES PROCESSORS ${nprocs} TIMEOUT 120)
  endforeach()
else()
  add_test(NAME amg_c_prec_extensions_1 COMMAND amg_c_prec_extensions)
endif()
