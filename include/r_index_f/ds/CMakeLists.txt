set(DS_SOURCES interval_block.hpp interval_pos.hpp idx_bit_vector.hpp idx_list.hpp symbol_map.hpp intervals_rank_w.hpp base_bv.hpp heads_wt_w.hpp)
add_library(ds OBJECT ${DS_SOURCES})
target_link_libraries(ds common sdsl) 