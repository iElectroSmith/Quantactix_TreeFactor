
cmake_minimum_required(VERSION 3.10)
project(PTreeCpp)

set(CMAKE_CXX_STANDARD 11)



aux_source_directory(. DIR_SRCS)
add_executable( PTreeCpp  ${DIR_SRCS})

target_link_libraries(PTreeCpp -larmadillo)


