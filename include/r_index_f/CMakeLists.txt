set(RIF  r-index-f.hpp)

add_library(r-index-f OBJECT ${MS_SOURCES})
target_link_libraries(r-index-f common sdsl)
target_include_directories(r-index-f PUBLIC "${r-index_SOURCE_DIR}/internal")