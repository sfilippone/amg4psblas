cmake_minimum_required(VERSION 3.10)
project(amg4psblas VERSION 1.0 LANGUAGES C Fortran)


set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Check for the installation path for psblas
if(NOT DEFINED PSBLAS_INSTALL_DIR)
    message(FATAL_ERROR "Please specify the path to the psblas installation directory using -DPSBLAS_INSTALL_DIR=<path>")
endif()
message(STATUS "psblas directory is ${PSBLAS_INSTALL_DIR};;")
# Set the include and library directories based on the provided path
#set(TEST_INSTALLDIR "${PSBLAS_INSTALL_DIR}")
set(INCDIR "${PSBLAS_INSTALL_DIR}/include")
set(MODDIR "${PSBLAS_INSTALL_DIR}/modules")
set(LIBDIR "${PSBLAS_INSTALL_DIR}/lib")


message(STATUS "PSBLAS DIRECTORY INC ${INCDIR}; MOD ${MODDIR}; LIB ${LIBDIR};")

# Find the psblas package
find_package(psblas REQUIRED PATHS ${PSBLAS_INSTALL_DIR})

if(NOT psblas_FOUND)
    message(FATAL_ERROR "PSBLAS not found!")
else()
    message(STATUS "Found PSBLAS: ${psblas_LIBRARIES}")
endif()

# Include directories for the project
include_directories(${PSBLAS_INSTALL_DIR})


# Include directories for the Fortran compiler
include_directories(${INCDIR} ${MODDIR} ${LIBDIR})



# Specify the installation directory

#set(${CMAKE_INSTALL_LIBDIR} "lib")
#message(STATUS "\t\t install libdir ${CMAKE_INSTALL_LIBDIR};")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${${CMAKE_PROJECT_NAME}_dist_string}-tests")
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${${CMAKE_PROJECT_NAME}_dist_string}-tests")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
#set(CMAKE_INSTALL_LIBDIR "lib" CACHE STRING "Library install directory")
#set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE STRING "Include directory")
#set(CMAKE_INSTALL_MODULDIR "modules" CACHE STRING "Module directory")

message(STATUS "Initial CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
set(AMG_CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

if(NOT AMG_CMAKE_INSTALL_LIBDIR)
    message(STATUS "CMAKE_INSTALL_LIBDIR is set to default value lib")
    #set(CMAKE_INSTALL_LIBDIR "lib" CACHE STRING "Library install directory" FORCE)
    set(CMAKE_INSTALL_LIBDIR ${PSB_CMAKE_INSTALL_LIBDIR})
    set(AMG_CMAKE_INSTALL_LIBDIR ${PSB_CMAKE_INSTALL_LIBDIR})
else()
    set(CMAKE_INSTALL_LIBDIR ${AMG_CMAKE_INSTALL_LIBDIR})
    message(STATUS "CMAKE_INSTALL_LIBDIR is set to: ${CMAKE_INSTALL_LIBDIR}")
endif()

if(NOT AMG_CMAKE_INSTALL_INCLUDEDIR)
    message(STATUS "CMAKE_INSTALL_INCLUDEDIR is set to default value lib")
    #set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE STRING "Include directory" FORCE)
    set(CMAKE_INSTALL_INCLUDEDIR ${PSB_CMAKE_INSTALL_INCLUDEDIR})
    set(AMG_CMAKE_INSTALL_INCLUDEDIR ${PSB_CMAKE_INSTALL_INCLUDEDIR})
else()
    set(CMAKE_INSTALL_INCLUDEDIR ${AMG_CMAKE_INSTALL_INCLUDEDIR})
    message(STATUS "CMAKE_INSTALL_INCLUDEDIR is set to: ${CMAKE_INSTALL_INCLUDEDIR}")
endif()

if(NOT AMG_CMAKE_INSTALL_MODULDIR)
    message(STATUS "CMAKE_INSTALL_MODULDIR is set to default value lib")
    #set(CMAKE_INSTALL_MODULDIR "modules" CACHE STRING "Modules directory" FORCE)
    set(CMAKE_INSTALL_MODULDIR ${PSB_CMAKE_INSTALL_MODULDIR})
    set(AMG_CMAKE_INSTALL_LIBDIR ${PSB_CMAKE_INSTALL_MODULDIR})
else()
    set(CMAKE_INSTALL_MODULDIR ${AMG_CMAKE_INSTALL_MODULDIR})
    message(STATUS "CMAKE_INSTALL_MODULDIR is set to: ${CMAKE_INSTALL_MODULDIR}")
endif()







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



  


#-----------------------------------------------------
# Publicize installed location to other CMake projects
#-----------------------------------------------------
#install(EXPORT ${CMAKE_PROJECT_NAME}-targets
#  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
#)


message(STATUS "NAME project ${CMAKE_PROJECT_NAME};")

install(EXPORT ${CMAKE_PROJECT_NAME}-targets
  FILE ${CMAKE_PROJECT_NAME}Config.cmake
  NAMESPACE ${CMAKE_PROJECT_NAME}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
)


include(CMakePackageConfigHelpers) # standard CMake module
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${amg4psblas_VERSION}"
  COMPATIBILITY SameMajorVersion
  )

configure_file("${CMAKE_SOURCE_DIR}/cmake/pkg/${CMAKE_PROJECT_NAME}Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${CMAKE_PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Targets.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}"
)







#------------------------------------------
# Add portable unistall command to makefile
#------------------------------------------
# Adapted from the CMake Wiki FAQ
configure_file ( "${CMAKE_SOURCE_DIR}/cmake/uninstall.cmake.in" "${CMAKE_BINARY_DIR}/uninstall.cmake"
  @ONLY)

add_custom_target ( uninstall
  COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/uninstall.cmake" )

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)
# See JSON-Fortran's CMakeLists.txt file to find out how to get the check target to depend
# on the test executables

#----------------------------------
# Determine if we're using Open MPI
#---------------------------------
if(MPI_FOUND)
  execute_process(COMMAND ${MPIEXEC} --version
    OUTPUT_VARIABLE mpi_version_out)
  if (mpi_version_out MATCHES "[Oo]pen[ -][Mm][Pp][Ii]")
    message( STATUS "OpenMPI detected")
    set ( openmpi true )
  endif()
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






install(DIRECTORY "${CMAKE_BINARY_DIR}/include" DESTINATION "include"
  FILES_MATCHING PATTERN "*.h")

install(DIRECTORY "${CMAKE_BINARY_DIR}/modules" DESTINATION "modules"
  FILES_MATCHING PATTERN "*.mod")






# Install the library
install(TARGETS amgprec amgcbind
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)



if(WIN32) #TODO
#  install(TARGETS psb_base_C
#    EXPORT ${CMAKE_PROJECT_NAME}-targets
#    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
#    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
#    )
#  if(METIS_FOUND)
#    install(TARGETS psb_util_C
#      EXPORT ${CMAKE_PROJECT_NAME}-targets
#      DESTINATION "${CMAKE_INSTALL_LIBDIR}"
#      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
#      )
#  endif()
endif()


message(STATUS "install directory is ${CMAKE_INSTALL_LIBDIR};;;")

# Step 2: Create the configuration file from the template
#configure_package_config_file(
#  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/amg4psblasConfig.cmake.in"
#  "${CMAKE_CURRENT_BINARY_DIR}/amg4psblasConfig.cmake"
#  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/amg4psblas"
#)

# Step 3: Install the generated config files
#install(FILES 
#  "${CMAKE_CURRENT_BINARY_DIR}/amg4psblasConfig.cmake"
#  "${CMAKE_CURRENT_BINARY_DIR}/amg4psblasConfigVersion.cmake"
#  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/amg4psblas"
#)

# Step 4: Export targets so that the build directory can be used directly
#export(
#  EXPORT ${CMAKE_PROJECT_NAME}-targets
#  FILE "${CMAKE_CURRENT_BINARY_DIR}/amg4psblasTargets.cmake"
#    NAMESPACE psblas::
#)
export(
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/psblasTargets.cmake"
    NAMESPACE ${CMAKE_PROJECT_NAME}::
)

export(
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Targets.cmake"
    NAMESPACE ${CMAKE_PROJECT_NAME}::
)

		
 

# Optionally, you can install the headers
#install(DIRECTORY include/
#    DESTINATION include
#)



message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} - ${PSB_CMAKE_INSTALL_PREFIX};")
message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR} - ${PSB_CMAKE_INSTALL_LIBDIR};")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR} - ${PSB_CMAKE_INSTALL_INCLUDEDIR};")
message(STATUS "CMAKE_INSTALL_MODULDIR: ${CMAKE_INSTALL_MODULDIR} - ${PSB_CMAKE_INSTALL_MODULDIR};")

