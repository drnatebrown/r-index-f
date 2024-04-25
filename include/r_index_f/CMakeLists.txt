add_subdirectory(ds)

set(RIF_SOURCES r_index_f.hpp columnar_table.hpp block_table.hpp LF_table.hpp alt_table.hpp hybrid_table.hpp)

add_library(r_index_f OBJECT ${RIF_SOURCES})
target_link_libraries(r_index_f ds common sdsl)