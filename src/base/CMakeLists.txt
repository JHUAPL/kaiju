file(GLOB base_srcs *.F90 defs/*.F90 types/*.F90 shellutils/*.F90 3rd-party/*.F90)
set(CMAKE_SUPPRESS_REGENERATION ON)
# ---- Generate git_info.F90 at build time in this directory's build tree ----
set(GIT_INFO_FILE "${CMAKE_CURRENT_BINARY_DIR}/git_info.F90")
message(STATUS "Git info will be generated at: ${GIT_INFO_FILE}")
message(STATUS "Git info will be generated from: ${CMAKE_CURRENT_BINARY_DIR}")

add_custom_command(
  OUTPUT "${GIT_INFO_FILE}"
  COMMAND "${CMAKE_COMMAND}"
          -DSRC_DIR=${CMAKE_SOURCE_DIR}
          -DOUT_FILE=${GIT_INFO_FILE}.tmp
          -P ${CMAKE_SOURCE_DIR}/cmake/git_info.cmake
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          ${GIT_INFO_FILE}.tmp ${GIT_INFO_FILE}
  COMMAND "${CMAKE_COMMAND}" -E rm -f ${GIT_INFO_FILE}.tmp
  # Re-run on commits/checkouts (covers branch switches + new commits)
  DEPENDS
    "${CMAKE_SOURCE_DIR}/.git/HEAD"
    "${CMAKE_SOURCE_DIR}/.git/logs/HEAD"
    "${CMAKE_SOURCE_DIR}/.git/index"
  VERBATIM
)

add_library(baselib ${base_srcs} "${GIT_INFO_FILE}")
set_source_files_properties("${GIT_INFO_FILE}" PROPERTIES GENERATED TRUE)