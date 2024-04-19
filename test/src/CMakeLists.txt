add_executable(build_rif build_rif.cpp)
target_link_libraries(build_rif sdsl malloc_count divsufsort divsufsort64)
target_include_directories(build_rif PUBLIC    "../../include/r_index_f" 
                                        "../../include/common" 
                                        "../../include/block_compression" 
                                        )
target_compile_options(build_rif PUBLIC "-std=c++17")

add_executable(rif_tests rif_tests.cpp)
target_link_libraries(rif_tests sdsl divsufsort divsufsort64 malloc_count)
target_include_directories(rif_tests PUBLIC    "../../include/r_index_f" 
                                        "../../include/common"
                                        "../../include/block_compression" 
                                        )
target_compile_options(rif_tests PUBLIC "-std=c++17")

add_executable(crazy_train crazy_train.cpp)
target_link_libraries(crazy_train sdsl divsufsort divsufsort64 malloc_count)
target_include_directories(crazy_train PUBLIC    "../../include/r_index_f" 
                                        "../../include/common"
                                        "../../include/block_compression" 
                                        )
target_compile_options(crazy_train PUBLIC "-std=c++17")