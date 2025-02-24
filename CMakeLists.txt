cmake_minimum_required(VERSION 3.10)
project(amg4psblas VERSION 1.0 LANGUAGES C Fortran)


set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Check for the installation path for psblas
if(NOT DEFINED PSBLAS_INSTALL_DIR)
    message(FATAL_ERROR "Please specify the path to the psblas installation directory using -DPSBLAS_INSTALL_DIR=<path>")
endif()

# Set the include and library directories based on the provided path
set(INSTALLDIR "${PSBLAS_INSTALL_DIR}")
set(INCDIR "${INSTALLDIR}/include")
set(MODDIR "${INSTALLDIR}/modules")
set(LIBDIR "${INSTALLDIR}/lib")


message(STATUS "PSBLAS DIRECTORY INC ${INCDIR}; MOD ${MODDIR}; LIB ${LIBDIR};")

# Find the psblas package
find_package(psblas REQUIRED PATHS ${INSTALLDIR})

if(NOT psblas_FOUND)
    message(FATAL_ERROR "PSBLAS not found!")
else()
    message(STATUS "Found PSBLAS: ${psblas_LIBRARIES}")
endif()

# Include directories for the project
include_directories(${PSBLAS_INSTALL_DIR})


# Include directories for the Fortran compiler
include_directories(${INCDIR} ${MODDIR} ${LIBDIR})












#TODO: what follow should be taken from psblas. We need to set in PSBLAS CMakeLists.txt the right variable to be used here
#------------------------------------
# Fortran name mangling introspection
#------------------------------------
include("${CMAKE_CURRENT_LIST_DIR}/cmake/CapitalizeString.cmake")
include(FortranCInterface)
CapitalizeString(${FortranCInterface_GLOBAL__CASE} fc_case)
message(STATUS "Name mangling capitalization: ${fc_case}")
message(STATUS "Name mangling fortran global suffix underscore: ${FortranCInterface_GLOBAL__SUFFIX}")
if(FortranCInterface_GLOBAL__SUFFIX STREQUAL "")
  add_compile_options("-D${fc_case}Case")
elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "_")
  add_compile_options("-D${fc_case}Underscore")
elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "__")
  add_compile_options("-D${fc_case}DoubleUnderscore")
else()
  message( FATAL_ERROR "Fortran name mangling suffix, \'${FortranCInterface_GLOBAL__SUFFIX}\', unknown to PSBLAS")
endif()


  message(STATUS "win? ${WIN32};")
if(TRUE )#NOT ${WIN32})
#previous check did not work if WIN32 is empty string
  #----------------------------------------------
  # Determine system endian-ness and pointer size
  #----------------------------------------------
  include(TestBigEndian)
  TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
  if(IS_BIG_ENDIAN)
    message( STATUS "System appears to be big endian.")
  else()
    message( STATUS "System appears to be little endian.")
    add_compile_options(-DLittleEndian)
  endif()
  include(CheckTypeSize)
  CHECK_TYPE_SIZE("void *" VOID_P_SIZE LANGUAGE C)
  if(${VOID_P_SIZE} EQUAL 8)
    add_compile_options(-DPtr64Bits)
  endif()
  message(STATUS "Have 64bit pointers")
  
  #add define values for integer size (IPKx) and long size (LPKx)
  CHECK_TYPE_SIZE("int" INT_SIZE LANGUAGE C)
  CHECK_TYPE_SIZE("long" LONG_SIZE LANGUAGE C)
  message(STATUS "INT SIZE ${INT_SIZE}")

  add_compile_options(-DIPK${INT_SIZE})
  add_compile_options(-DLPK${LONG_SIZE})

  
endif()



  











include(${CMAKE_CURRENT_LIST_DIR}/amgprec/CMakeLists.txt)  # include amgprec_source_files and amgprec_source_C_files source files

add_library(amgprec_C OBJECT ${amgprec_source_C_files})

target_link_libraries(amgprec_C
  #PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  #PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES}
  psblas::util psblas::linsolve psblas::prec psblas::ext psblas::cbind psblas::base)  #TODO  check actual libraries needed

add_library(amgprec ${amgprec_source_files} $<TARGET_OBJECTS:amgprec_C>)

set_target_properties(amgprec
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"  
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_amgprec
  LINKER_LANGUAGE Fortran
  )

target_include_directories(amgprec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules> 
  $<INSTALL_INTERFACE:modules>)

message(STATUS "include dir := ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")

#target_include_directories(base PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

target_include_directories(amgprec PUBLIC ${INCDIR} ${MODDIR})
target_link_libraries(amgprec
  #PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  #PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES}
  psblas::util psblas::linsolve psblas::prec psblas::ext psblas::cbind psblas::base)  #TODO  check actual libraries needed







include(${CMAKE_CURRENT_LIST_DIR}/cbind/CMakeLists.txt)  # include amgprec_source_files and amgprec_source_C_files source files


foreach(path IN LISTS amgcbind_header_C_files)
  # Copy the header file to the include directory
  file(COPY "${path}" DESTINATION "${CMAKE_BINARY_DIR}/include")
endforeach()


add_library(amgcbind_C OBJECT ${amgprec_source_C_files})

target_link_libraries(amgcbind_C
  #PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  #PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES}
  psblas::util psblas::linsolve psblas::prec psblas::ext psblas::cbind psblas::base)  #TODO  check actual libraries needed

add_library(amgcbind ${amgprec_source_files} $<TARGET_OBJECTS:amgcbind_C>)

set_target_properties(amgcbind
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"  
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_amgcbind
  LINKER_LANGUAGE Fortran
  )

target_include_directories(amgcbind PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules> 
  $<INSTALL_INTERFACE:modules>)

message(STATUS "include dir := ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")

#target_include_directories(base PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

target_include_directories(amgcbind PUBLIC ${INCDIR} ${MODDIR})
target_link_libraries(amgcbind
  #PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  #PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES}
  PUBLIC amgprec psblas::util psblas::linsolve psblas::prec psblas::ext psblas::cbind psblas::base)  #TODO  check actual libraries needed








# Install the library
install(TARGETS amgprec amgcbind
    DESTINATION lib
)

# Optionally, you can install the headers
install(DIRECTORY include/
    DESTINATION include
)
