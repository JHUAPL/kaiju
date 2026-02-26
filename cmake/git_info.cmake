# cmake/git_info.cmake

# Defaults
set(GIT_BRANCH "nogit")
set(GIT_COMMIT_HASH "nogit")

# Try git only if it looks like a git checkout
if(EXISTS "${SRC_DIR}/.git")
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY "${SRC_DIR}"
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _rc1
    ERROR_QUIET
  )
  execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY "${SRC_DIR}"
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _rc2
    ERROR_QUIET
  )
  if(NOT _rc1 EQUAL 0)
    set(GIT_BRANCH "nogit")
  endif()
  if(NOT _rc2 EQUAL 0)
    set(GIT_COMMIT_HASH "nogit")
  endif()
endif()

# Write out the git info as a Fortran module
file(WRITE  "${OUT_FILE}" "module git_info\n")
file(APPEND "${OUT_FILE}" "    implicit none\n")
file(APPEND "${OUT_FILE}" "    character(len=*), parameter :: gitB = '${GIT_BRANCH}'\n")
file(APPEND "${OUT_FILE}" "    character(len=*), parameter :: gitH = '${GIT_COMMIT_HASH}'\n")
file(APPEND "${OUT_FILE}" "end module git_info\n")