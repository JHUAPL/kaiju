#Handle compiler/defaults

#-------------
#Handle libraries (HDF5/OMP/MPI)
find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)
set(CMAKE_REQUIRED_INCLUDES ${HDF5_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
# h5 compiler helper can have issues with some preprocessor commands, just link libraries
link_libraries(${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})

find_package(OpenMP COMPONENTS Fortran)
set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_FLAGS})
set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_Fortran_LIBRARIES})

if (ENABLE_MPI)
    #mpi is a nightmare
    #try to explicitly find intel mpi first
    set(MPI_Fortran_COMPILER mpiifort)
    find_package(MPI COMPONENTS Fortran QUIET)
    find_program(MPIEXE mpiifort)
    if(NOT MPI_FOUND OR NOT MPI_Fortran_HAVE_F08_MODULE OR NOT MPIEXE)
        #try to find a fortran 2008 specific wrapper
        set(MPI_Fortran_COMPILER mpif08)
    	find_package(MPI COMPONENTS Fortran QUIET)
        unset(MPIEXE CACHE)
        unset(MPIEXE-NOTFOUND CACHE)
        find_program(MPIEXE mpif08)
        if(NOT MPI_FOUND OR NOT MPI_Fortran_HAVE_F08_MODULE OR NOT MPIEXE)
            #just look for whatever
            unset(MPI_Fortran_COMPILER)
            find_package(MPI REQUIRED COMPONENTS Fortran)
        endif()
    endif()
    if(MPI_FOUND AND MPI_Fortran_HAVE_F08_MODULE)
        message("-- Found MPI")
    else()
        message(FATAL_ERROR "Could not find an MPI Library that supports the F08 interface")
    endif()
endif()

#-------------
#Set default build to release
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

#-------------
#Set base release options
set(CMAKE_DEFOPT "-O3") #Default optimization

#-------------
#Set minimum compiler versions
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 17.0)
		message("Fortran compiler too old!  What, were you gonna use punch cards?")
		message(FATAL_ERROR "ifort > 17.0 required")
	elseif( (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 17.0) AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 18.0) )
		message(WARNING "Compiler has incomplete F2008 features, Git hash/compiler information won't be saved to H5 files")
		add_compile_definitions(__INTEL_COMPILER_OLD)
	elseif( (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 20.0) )
		message(WARNING "Setting default optimization to O2 to avoid certain Intel compiler bugs")
    	set(CMAKE_DEFOPT "-O2")
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 8.0)
		message("Fortran compiler too old!  What, were you gonna use punch cards?")
		message(FATAL_ERROR "gfortran > 8.0 required")
	endif()
endif()

if(FAST_YOLO)
    message("*** RUNNING WITH IMPRECISE HIGH SPEED OPTIONS ***")
    set(CMAKE_DEFOPT "-O3") #over-ride default optimization to maximum
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_DEFOPT}")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_DEFOPT} -g")

#Do compiler specific options
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	set(dialect "-free -implicitnone")
	#Base
	string(APPEND CMAKE_Fortran_FLAGS " -fPIC -fpconstant")
	#Production
    set(PROD "-align array64byte -align rec32byte -no-prec-div -fast-transcendentals")
    #Production with Debug Info
    set(PRODWITHDEBUGINFO "-traceback -debug all -align array64byte -align rec32byte -no-prec-div -fast-transcendentals")
	#Debug
        set(DEBUG "-traceback -debug all -gen-interfaces")
        if(NOT DISABLE_DEBUG_BOUNDS_CHECKS)
                string(APPEND DEBUG " -check bounds -check uninit -fp-stack-check -fstack-security-check")
        endif()
        if(NOT PFUNIT_FOUND)
                # Known bug with warning interfaces and PFUNIT
                string(APPEND DEBUG " -warn interfaces")
        endif()

	#Now do OS-dep options
	if (CMAKE_SYSTEM_NAME MATCHES Darwin)
		string(APPEND CMAKE_Fortran_FLAGS " -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 -xHost")
	else()
		#If we're not doing Mac, then add IPO
		string(APPEND PROD " -ipo")
		string(APPEND PRODWITHDEBUGINFO " -ipo")
	endif()

	#Handle individual hosts
	if (HOST MATCHES cheyenne)
		string(APPEND PROD " -march=corei7 -axCORE-AVX2")
		string(APPEND PRODWITHDEBUGINFO " -march=corei7 -axCORE-AVX2")
    elseif(HOST MATCHES pfe)
		string(APPEND PROD " -march=corei7 -axCORE-AVX2")
		string(APPEND PRODWITHDEBUGINFO " -march=corei7 -axCORE-AVX2")
	elseif (HOST MATCHES derecho)
		message("You're on Derecho, good for you!")
		if (ENABLE_MKL)
			string(APPEND CMAKE_Fortran_FLAGS " -qmkl")
		endif()
		string(APPEND PROD " -march=core-avx2")
		string(APPEND PRODWITHDEBUGINFO " -march=core-avx2")
	endif()

	#Check Intel Fortran version
	if(NOT ALLOW_INVALID_COMPILERS AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "2021.9")
		message(FATAL_ERROR "Intel Fortran compilers newer than 2023 (version 2021.8) are not supported. Set the ALLOW_INVALID_COMPILERS variable to ON to force compilation at your own risk.")
	endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	set(dialect "-ffree-form -ffree-line-length-none -fimplicit-none")
	#Base
	string(APPEND CMAKE_Fortran_FLAGS " -fPIC")
	#Production
	set(PROD "-ffast-math")
	#Debug
    if(NOT DISABLE_DEBUG_BOUNDS_CHECKS)
    	set(DEBUG "-fbacktrace -g -Warray-temporaries -Wall -Wfatal-errors -finit-local-zero")
    else()
        set(DEBUG "-fbacktrace -g -Warray-temporaries -Wall -Wfatal-errors")
    endif()
	#Now do machine-dep options
	if (CMAKE_SYSTEM_NAME MATCHES Darwin)
		string(APPEND CMAKE_Fortran_FLAGS " -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000")
	endif()
endif()

string(APPEND CMAKE_Fortran_FLAGS " ${dialect}")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " ${DEBUG}")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " ${PROD}")
string(APPEND CMAKE_Fortran_FLAGS_RELWITHDEBINFO " ${PRODWITHDEBUGINFO}")

if(ENABLE_OMP)
	string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
endif()
if(ENABLE_MPI)
    add_compile_options(${MPI_Fortran_COMPILE_OPTIONS})
	add_definitions(${MPI_Fortran_COMPILE_DEFINITIONS})
	include_directories(${MPI_Fortran_INCLUDE_DIRS})
	link_directories(${MPI_Fortran_LIBRARIES})
    string(APPEND CMAKE_Fortran_FLAGS " ${MPI_Fortran_LINK_FLAGS}")

    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel AND MPI_Fortran_COMPILER MATCHES mpiifort)
        #Using Intel Compiler and Intel MPI, use thread safe mpi compiler flag
        string(APPEND CMAKE_Fortran_FLAGS " -mt_mpi")
    else()
        #use different MPI link command
        string(APPEND CMAKE_Fortran_FLAGS " -lmpi")
    endif()

	set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
	# we changed compiler, link HDF5 libraries
	link_libraries(${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
endif()
if(ENABLE_CODECOV)
    # track code coverage
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/codecov_prof)
    string(APPEND CMAKE_Fortran_FLAGS " -prof-gen=srcpos -prof-dir=${CMAKE_BINARY_DIR}/bin/codecov_prof")
endif()

