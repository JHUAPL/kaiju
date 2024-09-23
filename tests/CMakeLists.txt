# All Tests use the lfmbw IC for Gamera
file(GLOB TGIC "${PROJECT_SOURCE_DIR}/src/gamera/ICs/lfmbw.F90")
# compile this Test Gamera IC into a library to be linked into the tests
add_library(tgiclib ${TGIC})
add_dependencies(tgiclib baselib gamlib)

# helper files for all testing code
file(GLOB helperFiles "${PROJECT_SOURCE_DIR}/tests/helperCode/*.F90")

# NON-MPI TESTS GO UP HERE
file(GLOB caseTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/cases/*.pf")
set(caseTestLibs baselib gamlib remixlib voltlib tgiclib)
add_pfunit_ctest (caseTests TEST_SOURCES ${caseTestFiles} OTHER_SOURCES ${helperFiles} LINK_LIBRARIES ${caseTestLibs})
add_dependencies(caseTests ${caseTestLibs})

file(GLOB gamTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/gamera/*.pf")
set(gamTestLibs baselib gamlib tgiclib)
add_pfunit_ctest (gamTests TEST_SOURCES ${gamTestFiles} OTHER_SOURCES ${helperFiles} LINK_LIBRARIES ${gamTestLibs})
add_dependencies(gamTests ${gamTestLibs})

file(GLOB mixTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/remix/*.pf")
set(mixTestLibs baselib gamlib remixlib voltlib tgiclib)
add_pfunit_ctest (mixTests TEST_SOURCES ${mixTestFiles} OTHER_SOURCES ${helperFiles} LINK_LIBRARIES ${mixTestLibs})
add_dependencies(mixTests ${mixTestLibs})

file(GLOB voltTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/voltron/*.pf")
set(voltTestLibs baselib gamlib remixlib voltlib tgiclib)
add_pfunit_ctest (voltTests TEST_SOURCES ${voltTestFiles} OTHER_SOURCES ${helperFiles} LINK_LIBRARIES ${voltTestLibs})
add_dependencies(voltTests ${voltTestLibs})

file(GLOB rcmTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/rcm/*.pf")
set(rcmTestLibs baselib gamlib voltlib remixlib rcmlib tgiclib)
add_pfunit_ctest(rcmTests TEST_SOURCES ${rcmTestFiles} OTHER_SOURCES ${helperFiles} LINK_LIBRARIES ${rcmTestLibs})
add_dependencies(rcmTests ${rcmTestLibs})

file(GLOB shgrTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/shellgrid/*.pf")
set(shgrTestLibs baselib)
add_pfunit_ctest (shgrTests TEST_SOURCES ${shgrTestFiles} OTHER_SOURCES ${helperFiles} LINK_LIBRARIES ${shgrTestLibs})
add_dependencies(shgrTests ${shgrTestLibs})

# copy input xml files and grid h5 files to binary folder for running
file(GLOB_RECURSE inputFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.xml")
foreach(inputFile ${inputFiles})
  get_filename_component(inputFileJustName ${inputFile} NAME)
  add_custom_command(
    OUTPUT "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${inputFileJustName}"
    COMMAND cmake -E copy
      "${inputFile}"
      "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${inputFileJustName}"
    DEPENDS "${inputFile}"
  )
  list(APPEND inputFilesDest "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${inputFileJustName}")
endforeach(inputFile)
add_custom_target(CopyInputs DEPENDS ${inputFilesDest})

file(GLOB_RECURSE gridFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h5")
foreach(gridFile ${gridFiles})
  get_filename_component(gridFileJustName ${gridFile} NAME)
  add_custom_command(
    OUTPUT "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${gridFileJustName}"
    COMMAND cmake -E copy
      "${gridFile}"
      "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${gridFileJustName}"
    DEPENDS "${gridFile}"
  )
  list(APPEND gridFilesDest "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${gridFileJustName}")
endforeach(gridFile)
add_custom_target(CopyGrids DEPENDS ${gridFilesDest})

if(ENABLE_MPI)
    #MPI TESTS GO IN HERE

    #helper files for mpi testing code
    file(GLOB helperFilesMpi "${PROJECT_SOURCE_DIR}/tests/helperCode_mpi/*.F90")

    file(GLOB caseMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/cases_mpi/*.pf")
    set(caseMpiTestLibs baselib gamlib remixlib voltlib basempilib gammpilib tgiclib)
    add_pfunit_ctest (caseMpiTests TEST_SOURCES ${caseMpiTestFiles} OTHER_SOURCES ${helperFiles} ${helperFilesMpi} LINK_LIBRARIES ${caseMpiTestLibs} MAX_PES 64)
    add_dependencies(caseMpiTests ${caseMpiTestLibs})

    file(GLOB baseMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/base_mpi/*.pf")
    set(baseMpiTestLibs baselib basempilib)
    add_pfunit_ctest (baseMpiTests TEST_SOURCES ${baseMpiTestFiles} OTHER_SOURCES ${helperFiles} ${helperFilesMpi} LINK_LIBRARIES ${baseMpiTestLibs} MAX_PES 64)
    add_dependencies(baseMpiTests ${baseMpiTestLibs})

    file(GLOB gamMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/gamera_mpi/*.pf")
    set(gamMpiTestLibs baselib gamlib basempilib gammpilib tgiclib)
    add_pfunit_ctest (gamMpiTests TEST_SOURCES ${gamMpiTestFiles} OTHER_SOURCES ${helperFiles} ${helperFilesMpi} LINK_LIBRARIES ${gamMpiTestLibs} MAX_PES 64)
    add_dependencies(gamMpiTests ${gamMpiTestLibs})

    file(GLOB voltMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/voltron_mpi/*.pf")
    set(voltMpiTestLibs baselib gamlib basempilib gammpilib voltlib voltmpilib remixlib tgiclib)
    add_pfunit_ctest (voltMpiTests TEST_SOURCES ${voltMpiTestFiles} OTHER_SOURCES ${helperFiles} ${helperFilesMpi} LINK_LIBRARIES ${voltMpiTestLibs} MAX_PES 17)
    add_dependencies(voltMpiTests ${voltMpiTestLibs})

endif()

# all tests in one
if(ENABLE_MPI)
    # include MPI tests
    set(allTestFiles ${caseTestFiles} ${gamTestFiles} ${mixTestFiles} ${voltTestFiles} ${rcmTestFiles} ${shgrTestFiles} ${caseMpiTestFiles} ${baseMpiTestFiles} ${gamMpiTestFiles} ${voltMpiTestFiles})
    set(allTestLibs baselib gamlib remixlib voltlib rcmlib basempilib gammpilib voltmpilib tgiclib)
    set(allTestBins caseTests gamTests mixTests voltTests rcmTests shgrTests gamMpiTests baseMpiTests caseMpiTests voltMpiTests)
    add_pfunit_ctest (allTests TEST_SOURCES ${allTestFiles} OTHER_SOURCES ${helperFiles} ${helperFilesMpi} LINK_LIBRARIES ${allTestLibs} MAX_PES 64)
else()
    # exclude MPI tests
    set(allTestFiles ${caseTestFiles} ${gamTestFiles} ${mixTestFiles} ${voltTestFiles} ${rcmTestFiles} ${shgrTestFiles})
    set(allTestLibs baselib gamlib remixlib voltlib rcmlib tgiclib)
    set(allTestBins caseTests gamTests mixTests voltTests rcmTests shgrTests)
    add_pfunit_ctest (allTests TEST_SOURCES ${allTestFiles} OTHER_SOURCES ${helperFiles} LINK_LIBRARIES ${allTestLibs})
endif()
add_dependencies(allTests CopyInputs CopyGrids ${allTestBins} ${allTestLibs})

