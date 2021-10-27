include(FetchContent)

## Add malloc_count
FetchContent_Declare(
    malloc_count
    GIT_REPOSITORY https://github.com/bingmann/malloc_count
)
  
FetchContent_GetProperties(malloc_count)
if(NOT malloc_count_POPULATED)
    FetchContent_Populate(malloc_count)

    add_library(malloc_count OBJECT ${malloc_count_SOURCE_DIR}/malloc_count.c ${malloc_count_SOURCE_DIR}/malloc_count.h)
    target_link_libraries(malloc_count dl)
    target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")

    add_library(memprofile OBJECT ${malloc_count_SOURCE_DIR}/memprofile.h)
     target_include_directories(memprofile PUBLIC "${malloc_count_SOURCE_DIR}")
endif()

## Add Big-BWT
FetchContent_Declare(
    bigbwt
    GIT_REPOSITORY https://github.com/alshai/Big-BWT.git
)
  
FetchContent_GetProperties(bigbwt)
if(NOT bigbwt_POPULATED)
    FetchContent_Populate(bigbwt)
    add_subdirectory(${bigbwt_SOURCE_DIR} ${bigbwt_BINARY_DIR})
endif()


  
## Add gsacak
FetchContent_Declare(
    gsacak
    GIT_REPOSITORY https://github.com/felipelouza/gsa-is.git
)

FetchContent_GetProperties(gsacak)
if(NOT gsacak_POPULATED)
    FetchContent_Populate(gsacak)
    add_library(gsacak OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak PUBLIC "${gsacak_SOURCE_DIR}")
    
    add_library(gsacak64 OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak64 PUBLIC "${gsacak_SOURCE_DIR}")
    target_compile_options(gsacak64 PUBLIC -DM64)
endif()


## Google benchmark
FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        master
)
  
 FetchContent_GetProperties(benchmark)
 if(NOT benchmark_POPULATED)
    FetchContent_Populate(benchmark)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable testing of the benchmark library." FORCE)
    add_subdirectory(${benchmark_SOURCE_DIR} ${benchmark_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()


## Add sdsl
FetchContent_Declare(
    sdsl
    # GIT_REPOSITORY https://github.com/elarielcl/sdsl-lite.git
    GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
    FetchContent_Populate(sdsl)

    set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")
  
    add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()