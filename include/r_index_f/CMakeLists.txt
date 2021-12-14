set(RIF_SOURCES r_index_f.hpp LF_table.hpp) 

add_library(r_index_f OBJECT ${RIF_SOURCES})
target_link_libraries(r_index_f block_compression common sdsl)