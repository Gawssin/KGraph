cmake_minimum_required(VERSION 3.12)
project(ForkGraph)

set(CMAKE_CXX_STANDARD 14)

include_directories(ligra)
include_directories(utils)
add_definitions("-DWEIGHTED")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -I ./ligra -march=native -O3 -DOPENMP -DBYTERLE")


               
add_executable(kgraph-sssp-sparse
                apps/kgraph-sssp-sparse.C) 

add_executable(kgraph-sssp-dense
                apps/kgraph-sssp-dense.C) 

add_executable(kgraph-sswp-sparse
                apps/kgraph-sswp-sparse.C) 

add_executable(kgraph-sswp-dense
                apps/kgraph-sswp-dense.C) 

add_executable(kgraph-ssr-sparse
                apps/kgraph-ssr-sparse.C) 

add_executable(kgraph-ssr-dense
                apps/kgraph-ssr-dense.C) 

add_executable(kgraph-deepwalk
                apps/kgraph-deepwalk.C) 

add_executable(kgraph-ppr
                apps/kgraph-ppr.C) 
