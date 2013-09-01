cmake_minimum_required (VERSION 2.8)
project (genops)

add_definitions( -std=c++11 )


set(GENOPS_SRC differ.cpp
               generator.cpp
               population.cpp)


add_library(genops STATIC ${GENOPS_SRC})

