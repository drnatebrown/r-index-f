set(BLOCK_SOURCES  block_table.hpp inteval_block.hpp interval_pos.hpp)

add_library(block_compression OBJECT ${BLOCK_SOURCES})
target_link_libraries(block_compression common sdsl)