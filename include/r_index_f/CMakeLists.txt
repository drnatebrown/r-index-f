set(RIF_SOURCES  r_index_f.hpp)

add_library(r_index_f OBJECT ${RIF_SOURCES})
target_link_libraries(r_index_f common sdsl)