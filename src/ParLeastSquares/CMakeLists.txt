add_library(ParLeastSquares SHARED par_ls.cpp)
target_include_directories(ParLeastSquares PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ParLeastSquares PUBLIC Eigen3::Eigen)
