cmake_minimum_required(VERSION 3.20)
project(LICO)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )

# add_compile_definitions(HUGEPAGE) # Notably, to enable hugepage, you need check if your system supports hugepage and set up hugepage properly. You can check the hugepage usage with "cat /proc/meminfo | grep Huge" command
# add_compile_definitions(PERF_OUTPUT) # Notably, to enable perf test, you need run this script with sudo, otherwise the perf test will fail due to insufficient permissions
add_compile_definitions(SIMD_512_2D1S)
# add_compile_definitions(SIMD_256_2D1S)

add_compile_definitions(OptimalPLA)
#add_compile_definitions(SwingFilterPLA)
#add_compile_definitions(LinearSpline)


set(CMAKE_CXX_FLAGS_NDEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")

# simd intersection
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/external/Vp2SetOperation/include)

# sdsl-lite
include_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/include)
link_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/lib)

# FastPFor
include_directories(${PROJECT_SOURCE_DIR}/external/FastPFor/headers)
link_directories(${PROJECT_SOURCE_DIR}/external/FastPFor/build)

# libdivide
include_directories(${PROJECT_SOURCE_DIR}/external/libdivide)

# mm_file
include_directories(${PROJECT_SOURCE_DIR}/external/mm_file/include/mm_file)


function(add_variants target_base source_file)
  set(target_z "${target_base}")
  add_executable(${target_z} ${source_file})
  target_compile_definitions(${target_z} PRIVATE RESIDUAL_COMPRESS=1)

  target_compile_options(${target_z} PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:
      -mlzcnt
      -mbmi2
      -mavx512f
      -mavx512vl
      -mavx512cd
      -mavx512dq
    >
  )

  target_link_libraries(${target_z} sdsl FastPFOR)

  set(target_p "${target_base}_none")
  add_executable(${target_p} ${source_file})
  target_compile_definitions(${target_p} PRIVATE RESIDUAL_COMPRESS=0)

  target_compile_options(${target_p} PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:
      -mlzcnt
      -mbmi2
      -mavx512f
      -mavx512vl
      -mavx512cd
      -mavx512dq
    >
  )

  target_link_libraries(${target_p} sdsl)
endfunction()


function(add_uncompress target_base source_file)
  set(target_p "${target_base}")
  add_executable(${target_p} ${source_file})
  set_target_properties(${target_p} PROPERTIES
          COMPILE_DEFINITIONS "RESIDUAL_COMPRESS=0" )
  target_link_libraries(${target_p} sdsl FastPFOR)
endfunction()


add_variants(lico_build lico_build.cpp)

add_variants(lico_decode lico_decode.cpp)

add_variants(lico_query lico_query.cpp)

#add_uncompress(other other_index.cpp)