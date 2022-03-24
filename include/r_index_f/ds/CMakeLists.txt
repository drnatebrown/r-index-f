set(DS_SOURCES interval_block.hpp interval_pos.hpp idx_bit_vector.hpp idx_list.hpp symbol_map.hpp intervals_rank_w.hpp base_bv.hpp base_interpolate.hpp base_sample.hpp heads_wt_w.hpp ACGT_map.hpp, heads_bv_w.hpp scan_list.hpp)
add_library(ds OBJECT ${DS_SOURCES})
target_link_libraries(ds common sdsl) 