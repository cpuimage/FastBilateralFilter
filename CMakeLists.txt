cmake_minimum_required(VERSION 3.9)
project(fbf)

add_executable(fbf fastbf_main.c mt19937ar.c affinity.c arrayalloc.c noisycomputations.c imageio.c maxfilter.c young.c deriche_o3opt.c gnuplot_i.c fastbf.c)
target_link_libraries(fbf -lm)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()