file(GLOB base_srcs *.F90 defs/*.F90 transforms/*.F90 types/*.F90 shellutils/*.F90)
add_library(baselib ${base_srcs})
