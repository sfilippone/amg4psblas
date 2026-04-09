cmake_minimum_required(VERSION 3.10)
project(amg4psblas VERSION 1.0 LANGUAGES C CXX Fortran)


set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


set(PSBLAS_INSTALL_DIR "" CACHE PATH "Path to the PSBLAS installation
directory")
if(PSBLAS_INSTALL_DIR STREQUAL "")
     message(FATAL_ERROR "Please specify the path to the PSBLAS installation directory using -DPSBLAS_INSTALL_DIR=<path> or set it in ccmake.")
endif()


# Check for the installation path for psblas
message(STATUS "psblas directory is ${PSBLAS_INSTALL_DIR};")

message(STATUS "PSBLAS DIRECTORY INC ${INCDIR}; MOD ${MODDIR}; LIB ${LIBDIR};")

#set(CMAKE_CXX_STANDARD 17) # Set cxx standard for the c++ part of the library


# Find the psblas package
find_package(psblas REQUIRED PATHS ${PSBLAS_INSTALL_DIR})

if(NOT psblas_FOUND)
    message(FATAL_ERROR "PSBLAS not found!")
else()
    message(STATUS "Found PSBLAS: ${PSBLAS_LIBRARIES}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Add -g to the Fortran compiler flags.
    # We use STRING(APPEND) to ensure we don't overwrite other important flags.
    string(APPEND CMAKE_Fortran_FLAGS " -g")
    string(APPEND CMAKE_CXX_FLAGS " -g")
    message(STATUS "Fortran and CXX debug flags added: -g")
endif()

string(APPEND CMAKE_Fortran_FLAGS " -O2")
string(APPEND CMAKE_CXX_FLAGS " -O2")
message(STATUS "Fortran and CXX optimization flags added: -O2")



# Set the include and library directories based on the provided path
#set(TEST_INSTALLDIR "${PSBLAS_INSTALL_DIR}")
set(INCDIR "${PSBLAS_INSTALL_DIR}/${PSB_CMAKE_INSTALL_INCLUDEDIR}")
set(MODDIR "${PSBLAS_INSTALL_DIR}/${PSB_CMAKE_INSTALL_MODULDIR}")
set(LIBDIR "${PSBLAS_INSTALL_DIR}/${PSB_CMAKE_INSTALL_LIBDIR}")


# Include directories for the project
include_directories(${PSBLAS_INSTALL_DIR} ${MPI_INCLUDE_PATH} )


# Include directories for the Fortran compiler
include_directories(${INCDIR} ${MODDIR} ${LIBDIR})



message(STATUS "Using IPK size: ${PSB_IPK_SIZE}")
message(STATUS "Using LPK size: ${PSB_LPK_SIZE}")

# Add PSB_IPK/LPK flag only for fortran files.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_IPK${PSB_IPK_SIZE}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_LPK${PSB_LPK_SIZE}")




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

configure_file("${CMAKE_SOURCE_DIR}/cmake/${CMAKE_PROJECT_NAME}Config.cmake.in"
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



find_package( MPI REQUIRED Fortran C CXX )

if(MPI_FOUND)
  #-----------------------------------------------
  # Work around an issue present on fedora systems
  #-----------------------------------------------
  if( (MPI_CXX_LINK_FLAGS MATCHES "noexecstack") OR (MPI_Fortran_LINK_FLAGS MATCHES "noexecstack") )
    message ( WARNING
      "The `noexecstack` linker flag was found in the MPI_<lang>_LINK_FLAGS variable. This is
known to cause segmentation faults for some Fortran codes. See, e.g.,
https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71729 or
https://github.com/sourceryinstitute/OpenCoarrays/issues/317.

`noexecstack` is being replaced with `execstack`"
      )
    string(REPLACE "noexecstack"
      "execstack" MPI_CXX_LINK_FLAGS_FIXED ${MPI_CXX_LINK_FLAGS})
    string(REPLACE "noexecstack"
      "execstack" MPI_C_LINK_FLAGS_FIXED ${MPI_C_LINK_FLAGS})
    string(REPLACE "noexecstack"
      "execstack" MPI_Fortran_LINK_FLAGS_FIXED ${MPI_Fortran_LINK_FLAGS})
    set(MPI_CXX_LINK_FLAGS "${MPI_CXX_LINK_FLAGS_FIXED}" CACHE STRING
      "MPI CXX linking flags" FORCE)
    set(MPI_C_LINK_FLAGS "${MPI_C_LINK_FLAGS_FIXED}" CACHE STRING
      "MPI C linking flags" FORCE)
    set(MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS_FIXED}" CACHE STRING
      "MPI Fortran linking flags" FORCE)
      
  endif()

    message(STATUS "Found MPI: ${MPI_C_LIBRARIES} - ${MPI_CXX_LIBRARIES} - ${MPI_Fortran_LIBRARIES}")
    
    
  #----------------
  # Setup MPI compilers
  #----------------  
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "C compiler" FORCE)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "C++ compiler" FORCE)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} CACHE FILEPATH "Fortran compiler" FORCE)
    
  #----------------
  # Setup MPI flags
  #----------------
  list(REMOVE_DUPLICATES MPI_Fortran_INCLUDE_PATH)
  set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
  set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
  set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
  set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
  include_directories(BEFORE ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
message(STATUS "${MPI_C_INCLUDE_PATH}; ${MPI_Fortran_INCLUDE_PATH};; ${CMAKE_Fortran_LINK_FLAGS} ;")
  if(MPI_Fortran_HAVE_F90_MODULE OR MPI_Fortran_HAVE_F08_MODULE)
    add_compile_options(-DPSB_MPI_MOD)
    message(STATUS "-DPSB_MPI_MOD")
    #add_compile_options(-DSERIAL_MPI) # Is it right??
    #message(STATUS "-DSERIAL_MPI")
  endif()
  set(PSB_SERIAL_MPI OFF)
  
else()
  message(STATUS "MPI  not found, serial ahead")
  add_compile_options(-DPSB_SERIAL_MPI)
  add_compile_options(-DPSB_MPI_MOD)
  set(PSB_SERIAL_MPI ON)
  set(CSERIALMPI "#define PSB_SERIAL_MPI")
endif()

add_compile_options(-O3)
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-frecursive>)

if(MPI_FOUND)
  execute_process(COMMAND ${MPIEXEC} --version
    OUTPUT_VARIABLE mpi_version_out)
  if (mpi_version_out MATCHES "[Oo]pen[ -][Mm][Pp][Ii]")
    message( STATUS "OpenMPI detected")
    set ( openmpi true )
  endif()
  
  
  set(MPI_H_COPIED FALSE)
  set(MPI_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include") # Define the include directory

  # Create the include directory if it doesn't exist
  file(MAKE_DIRECTORY "${MPI_INCLUDE_DIR}")

  foreach(path IN LISTS MPI_INCLUDE_PATH)
    # Construct the full path to the mpi.h file
    set(mpi_h_path "${path}/mpi.h")
    
    # Check if the mpi.h file exists
    if(EXISTS "${mpi_h_path}")
        # Copy the mpi.h file to the include directory
        file(COPY "${mpi_h_path}" DESTINATION "${MPI_INCLUDE_DIR}")
        message(STATUS "Copied mpi.h from ${mpi_h_path} to ${MPI_INCLUDE_DIR}")
        set(MPI_H_COPIED TRUE)
        break()  # Exit the loop once we've copied the file
    endif()
  endforeach()

  if(NOT MPI_H_COPIED)
    message(WARNING "mpi.h not found in any of the specified paths: ${MPI_INCLUDE_PATH}")
  endif()

  # Add the created include directory to the project's include directories
  #include_directories("${MPI_INCLUDE_DIR}")
endif()

# Find AMG constants
include(${CMAKE_CURRENT_LIST_DIR}/cmake/readAMGConst.cmake)
_amg_read_const()

if ("${PSB_LPK_SIZE}" EQUAL 8) 
   set(CXXMATCHBOXBIT "#define AMG_MATCHBOXP_BIT64")
endif()


#------------------------------------------
# Configure the amg_config.h file
#------------------------------------------

message(STATUS "bin dir ${CMAKE_CURRENT_BINARY_DIR}; source dir ${CMAKE_CURRENT_SOURCE_DIR};;")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/amgprec/amg_config.h.in 
    ${CMAKE_CURRENT_BINARY_DIR}/include/amg_config.h 	
    @ONLY  # Replace variables only
)



#---------------------------------------
# Add the AMG libraries
#---------------------------------------

# In your CMakeLists.txt
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-256")

message(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")



include(${CMAKE_CURRENT_LIST_DIR}/amgprec/CMakeLists.txt)  # include amgprec_source_files and amgprec_source_C_files source files

include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")

add_library(amgprec_C OBJECT ${amgprec_source_C_files})
add_library(amgprec_CPP OBJECT ${amgprec_source_CPP_files})

target_link_libraries(amgprec_C
  #PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  #PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES}
  psblas::util psblas::linsolve psblas::prec psblas::ext psblas::cbind psblas::base
  #${MPI_C_LIBRARIES}
  )  #TODO  check actual libraries needed

target_link_libraries(amgprec_CPP
  #PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  #PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES}
  psblas::util psblas::linsolve psblas::prec psblas::ext psblas::cbind psblas::base
  stdc++
  ${MPI_CXX_LIBRARIES})  #TODO  check actual libraries needed

add_library(amgprec ${amgprec_source_files} $<TARGET_OBJECTS:amgprec_CPP> $<TARGET_OBJECTS:amgprec_C> )

set_target_properties(amgprec
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"  
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME amg_prec
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
  psblas::util psblas::linsolve psblas::prec psblas::ext psblas::cbind psblas::base
  #${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES}
  )  #TODO  check actual libraries needed







include(${CMAKE_CURRENT_LIST_DIR}/cbind/CMakeLists.txt)  # include amgprec_source_files and amgprec_source_C_files source files


foreach(path IN LISTS amgcbind_header_C_files)
  # Copy the header file to the include directory
  file(COPY "${path}" DESTINATION "${CMAKE_BINARY_DIR}/include")
endforeach()


add_library(amgcbind_C OBJECT ${amgcbind_source_C_files})

target_link_libraries(amgcbind_C
  #PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  #PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES}
  psblas::util psblas::linsolve psblas::prec psblas::ext psblas::cbind psblas::base)  #TODO  check actual libraries needed

add_library(amgcbind ${amgcbind_source_files} $<TARGET_OBJECTS:amgcbind_C>)

set_target_properties(amgcbind
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"  
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME amg_cbind
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
  PUBLIC amgprec psblas::util psblas::linsolve psblas::prec
  psblas::ext psblas::cbind psblas::base)
  #TODO  check actual libraries needed






install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.h")

install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/ DESTINATION "${CMAKE_INSTALL_MODULDIR}"
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
  FILE "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Targets.cmake"
    NAMESPACE ${CMAKE_PROJECT_NAME}::
)

#export(
#  EXPORT ${CMAKE_PROJECT_NAME}-targets
#  FILE "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Targets.cmake"
#    NAMESPACE ${CMAKE_PROJECT_NAME}::
#)

		
 

# Set the installation directory for the test files
set(INSTALL_TEST_DIR "${CMAKE_INSTALL_PREFIX}/samples" CACHE PATH "Installation directory for sample files")

function(install_directory_recursive source_dir install_base_dir) # Function to install a directory and its subdirectories recursively
    file(GLOB_RECURSE ALL_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/${source_dir}" "${source_dir}/*")

    foreach(FILE_PATH IN LISTS ALL_FILES)
        # Construct the full source and destination paths
        set(FULL_SOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${source_dir}/${FILE_PATH}")
        set(FULL_INSTALL_PATH "${install_base_dir}/${FILE_PATH}")

        # Check if it's a directory
        if(IS_DIRECTORY "${FULL_SOURCE_PATH}")
            # Create the directory in the install destination
            file(MAKE_DIRECTORY "${FULL_INSTALL_PATH}")
        else()
            # Install the file
            install(FILES "${FULL_SOURCE_PATH}" DESTINATION "${install_base_dir}" RENAME "${FILE_PATH}")
        endif()
    endforeach()
endfunction()



# Install test/fileread directory
install_directory_recursive(samples/simple "${INSTALL_TEST_DIR}/simple")

# Install test/pdegen directory
install_directory_recursive(samples/advanced "${INSTALL_TEST_DIR}/advanced")



message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} - ${PSB_CMAKE_INSTALL_PREFIX};")
message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR} - ${PSB_CMAKE_INSTALL_LIBDIR};")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR} - ${PSB_CMAKE_INSTALL_INCLUDEDIR};")
message(STATUS "CMAKE_INSTALL_MODULDIR: ${CMAKE_INSTALL_MODULDIR} - ${PSB_CMAKE_INSTALL_MODULDIR};")

