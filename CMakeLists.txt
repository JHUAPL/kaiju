cmake_minimum_required(VERSION 3.20.2)
project(Kaiju Fortran)

# use the new version of CMP0074, this tells cmake to use <PACKAGE>_ROOT environment
#   variables when looking for packages
cmake_policy(SET CMP0074 NEW)

# add and search for pfunit (fingers crossed)
list(APPEND CMAKE_PREFIX_PATH "./external")
find_package(PFUNIT QUIET)
if(PFUNIT_FOUND)
	enable_testing()
endif()

#Set hostname & configuration defaults before going to compilers.cmake
cmake_host_system_information(RESULT HOST QUERY HOSTNAME)

#Get the current working branch
execute_process(COMMAND git rev-parse --abbrev-ref HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE GIT_BRANCH
	OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND git log -1 --format=%h WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE GIT_COMMIT_HASH
	OUTPUT_STRIP_TRAILING_WHITESPACE)

# Generate git_info.F90 with GITBRANCH value
set(GIT_INFO_FILE "${CMAKE_SOURCE_DIR}/src/base/git_info.F90")
file(WRITE ${GIT_INFO_FILE} "module git_info\n")
file(APPEND ${GIT_INFO_FILE} "    implicit none\n")
file(APPEND ${GIT_INFO_FILE} "    character(len=*), parameter :: gitB = '${GIT_BRANCH}'\n")
file(APPEND ${GIT_INFO_FILE} "    character(len=*), parameter :: gitH = '${GIT_COMMIT_HASH}'\n")
file(APPEND ${GIT_INFO_FILE} "end module git_info\n")

#Set some configuration defaults
option(ENABLE_OMP "Enable OMP parallelization" ON )
option(ENABLE_MPI "Enable MPI parallelization" OFF)
option(ENABLE_MKL "Enable Intel MKL usage" OFF)
option(ENABLE_COLORTEXT "Enable ascii color in text output" ON)
option(ENABLE_CODECOV "Enable Code Coverage tracking" OFF)
option(ALLOW_INVALID_COMPILERS "Allow Invalid Compiler Versions" OFF)

#Vector blocking
set(bricksize 16 CACHE STRING "Vector brick size")
set_property(CACHE bricksize PROPERTY STRINGS 8 16 32)
add_compile_definitions(BRICKSIZE=${bricksize})
mark_as_advanced(bricksize)

#Hard-coded buffer for holding flux tubes
set(MAXTUBESIZE 5000 CACHE STRING "Size of array buffer for holding flux tubes")
add_compile_definitions(MAXTUBESIZE=${MAXTUBESIZE})

#RCM grid size (required to be defined at compile time)
set(RCMSIZEI 180 CACHE STRING "RCM Grid: Number of Colatitude Cells")
set(RCMSIZEJ 361 CACHE STRING "RCM Grid: Number of Longitude Cells")
set(RCMSIZEK 160 CACHE STRING "RCM Grid: Number of Lambda Cells")
set(RCMWRAPJ 3   CACHE STRING "RCM Grid: Number of periodic wrap cells")
mark_as_advanced(RCMWRAPJ)
add_compile_definitions(RCMSIZEI=${RCMSIZEI})
add_compile_definitions(RCMSIZEJ=${RCMSIZEJ})
add_compile_definitions(RCMSIZEK=${RCMSIZEK})
add_compile_definitions(RCMWRAPJ=${RCMWRAPJ})

#-------------
#Set default IC files
file(GLOB GAMICDEF src/gamera/ICs/null.F90)
set(GAMIC ${GAMICDEF} CACHE FILEPATH "Gamera IC file")

file(GLOB VOLTICDEF src/voltron/ICs/earthcmi.F90)
set(VOLTIC ${VOLTICDEF} CACHE FILEPATH "Voltron IC file")

file(GLOB TPICDEF src/chimp/tpICs/tpICstd.F90)
set(TPIC ${TPICDEF} CACHE FILEPATH "CHIMP TP IC file")

file(GLOB EBICDEF src/chimp/ebICs/ebICstd.F90)
set(EBIC ${EBICDEF} CACHE FILEPATH "CHIMP EB IC file")

file(GLOB GAMHELICDEF src/gamera/ICs/helio/wsa.F90)
set(GAMHELIC ${GAMHELICDEF} CACHE FILEPATH "Gamera Helio IC file")

#-------------
#Include user cmake file if present
include(cmake/user.cmake OPTIONAL)

if(ENABLE_MPI)
    if(MPI_BASE_ADDR_SIZE)
        add_compile_definitions(MPI_BASE_ADDR_SIZE=${MPI_BASE_ADDR_SIZE})
        message("Over-riding default MPI Address Size with new value of ${MPI_BASE_ADDR_SIZE} for basic MPI functions")
    else()
        message("Using default MPI Address Size of MPI_ADDRESS_KIND for basic MPI functions")
    endif()
    if(MPI_AN_ADDR_SIZE)
        add_compile_definitions(MPI_AN_ADDR_SIZE=${MPI_AN_ADDR_SIZE})
        message("Over-riding default MPI Address Size with new value of ${MPI_AN_ADDR_SIZE} for neighborhood MPI functions")
    else()
        message("Using default MPI Address Size of MPI_ADDRESS_KIND for neighborhood MPI functions")
    endif()
endif()

if(ENABLE_MKL)
    add_compile_definitions(USEMKL)
endif()
if(ENABLE_COLORTEXT)
    add_compile_definitions(USECOLORTEXT)
endif()

#Include compile flag generator
include(cmake/compilers.cmake)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#-------------
#Print info
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)

message("-------------------------")
message("Configuration summary ...")
message("System: ${HOST}")
message("OS: ${CMAKE_SYSTEM_NAME}")
message("Processor: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message("Compiler: ${CMAKE_Fortran_COMPILER_ID} / ${CMAKE_Fortran_COMPILER_VERSION}")
message("          ${CMAKE_Fortran_COMPILER}")
message("HDF5 Wrapper: ${HDF5_Fortran_COMPILER_EXECUTABLE}")
message("Version: ${GIT_COMMIT_HASH} / ${GIT_BRANCH}")
message("Build Type: ${CMAKE_BUILD_TYPE}")
message("Base Flags: ${CMAKE_Fortran_FLAGS}")
message("Build Flags: ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
message("-------------------------\n")
#message("\tBinary directory: ${PROJECT_BINARY_DIR}")

#-------------
#Start doing real work

#Start w/ base
add_subdirectory(src/base)
add_custom_target(base ALL)
add_executable(kaitoy.x src/drivers/kaitoyx.F90)
target_link_libraries(kaitoy.x baselib)
add_dependencies(base kaitoy.x)

#-------------
#Kaiju: CHIMP
message("Adding CHIMP module ...")
add_subdirectory(src/chimp)
#Chimp executables, target "chimp"
add_custom_target(chimp ALL)
message("\tEB IC file is ${EBIC}")
message("\tTP IC file is ${TPIC}")
set(chimpX "project;psd;push;slice;chop;trace;sctrack;calcdb;wpicheck")
foreach (EXE ${chimpX})
	message("\tAdding executable ${EXE}.x")
	add_executable(${EXE}.x src/drivers/${EXE}x.F90)
	target_link_libraries(${EXE}.x chimplib baselib)
	add_dependencies(chimp ${EXE}.x)
endforeach()

#-------------
#Kaiju: Gamera
message("Adding Gamera module ...")

#Add source
add_subdirectory(src/gamera)

#Print gamera info
message("\tBricksize is ${bricksize}")
message("\tIC file is ${GAMIC}")
add_custom_target(gamera ALL)
message("\tAdding executable gamera.x")
add_executable(gamera.x src/drivers/gamerax.F90 ${GAMIC})
target_link_libraries(gamera.x gamlib baselib)
add_dependencies(gamera gamera.x)

#-------------
#Kaiju: Gamera helio
message("Adding Gamera Helio module ...")

#Add source
#add_subdirectory(src/gamera)

#Print gamera helio info
#message("\tBricksize is ${bricksize}")
message("\tIC file is ${GAMHELIC}")
add_custom_target(gamhelio ALL)
message("\tAdding executable gamhelio.x")
add_executable(gamhelio.x src/drivers/gamerax.F90 ${GAMHELIC})
target_link_libraries(gamhelio.x gamlib baselib)
add_dependencies(gamera gamhelio.x)

#-------------
#Kaiju: ReMIX
message("Adding ReMIX module ...")
#Add source
add_subdirectory(src/remix)
add_custom_target(remix ALL)
message("\tAdding executable remix.x")
add_executable(remix.x src/drivers/remix.F90)
add_executable(remix2remix.x src/drivers/remix2remix.F90)
add_executable(remix2rcm.x src/drivers/remix2rcm.F90)
target_link_libraries(remix.x       remixlib baselib)
target_link_libraries(remix2remix.x remixlib baselib)
target_link_libraries(remix2rcm.x   remixlib baselib)
add_dependencies(remix remix.x remix2remix.x remix2rcm.x baselib)

#-------------
# ShellGrid Tester
message("\tAdding executable testNewRM.x")
add_executable(testNewRMS.x src/drivers/testNewRMS.F90)
target_link_libraries(testNewRMS.x baselib)
add_dependencies(testNewRMS.x baselib)

#-------------
#Kaiju: RCM
message("Adding RCM module ...")
message("\tRCM Grid is of size ${RCMSIZEI} x ${RCMSIZEJ} x ${RCMSIZEK}")
#Add source
add_subdirectory(src/rcm)
add_custom_target(rcm ALL)
message("\tAdding executable rcm.x")
add_executable(rcm.x src/drivers/rcmx.F90)
target_link_libraries(rcm.x rcmlib baselib)
add_dependencies(rcm rcm.x baselib)

#-------------
#Kaiju: Voltron
message("Adding Voltron module ...")
message("\tIC file is ${VOLTIC}")
add_subdirectory(src/voltron)
add_custom_target(voltron ALL)
message("\tAdding executable voltron.x")
add_executable(voltron.x src/drivers/voltronx.F90)
target_link_libraries(voltron.x baselib voltlib gamlib remixlib chimplib rcmlib)
add_dependencies(voltron voltron.x)

if(ENABLE_MPI)
    #-------------
    #Kaiju: Base MPI
    message("Adding Base MPI module ...")
    add_subdirectory(src/base/mpi)
    add_custom_target(base_mpi ALL)
    message("\tAdding executable kaitoy_mpi.x")
    add_executable(kaitoy_mpi.x src/drivers/kaitoy_mpix.F90)
    target_link_libraries(kaitoy_mpi.x baselib basempilib)
    add_dependencies(base_mpi kaitoy_mpi.x)

    #-------------
    #Kaiju: Gamera MPI
    message("Adding Gamera MPI module ...")
    add_subdirectory(src/gamera/mpi)
    add_custom_target(gamera_mpi ALL)
    message("\tAdding executable gamera_mpi.x")
    add_executable(gamera_mpi.x src/drivers/gamera_mpix.F90 ${GAMIC})
    target_link_libraries(gamera_mpi.x baselib gamlib basempilib gammpilib chimplib)
    add_dependencies(gamera_mpi gamera_mpi.x)

    #-------------
    #Kaiju: Voltron MPI
    message("Adding Voltron MPI module ...")
    add_subdirectory(src/voltron/mpi)
    add_custom_target(voltron_mpi ALL)
    message("\tAdding executable voltron_mpi.x")
    add_executable(voltron_mpi.x src/drivers/voltron_mpix.F90)
    target_link_libraries(voltron_mpi.x voltmpilib voltlib gammpilib gamlib remixlib chimplib rcmlib basempilib baselib)
    add_dependencies(voltron_mpi voltron_mpi.x)

    #-------------
    #Kaiju: Gamera Helio MPI
    message("Adding Gamera Helio MPI module ...")
    add_custom_target(gamhelio_mpi ALL)
    message("\tAdding executable gamhelio_mpi.x")
    add_executable(gamhelio_mpi.x src/drivers/gamera_mpix.F90 ${GAMHELIC})
    target_link_libraries(gamhelio_mpi.x baselib gamlib basempilib gammpilib)
    add_dependencies(gamhelio_mpi gamhelio_mpi.x)
    
endif()

if(PFUNIT_FOUND)
	#-------------
	#Kaiju: tests
	message("Adding tests module ...")
	#Add source
	add_subdirectory(tests)
	add_dependencies(allTests baselib gamlib remixlib chimplib voltlib)
endif()

# copy git hook scripts into .git/hooks
file(GLOB HOOK_SCRIPTS gitHookScripts/*)
if (HOOK_SCRIPTS)
    file(COPY ${HOOK_SCRIPTS} DESTINATION ${CMAKE_SOURCE_DIR}/.git/hooks) 
endif (HOOK_SCRIPTS)

