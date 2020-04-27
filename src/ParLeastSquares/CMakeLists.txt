configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/par_ls_defs.hpp.in 
    ${CMAKE_CURRENT_BINARY_DIR}/par_ls_defs.hpp)

add_library(ParLeastSquares SHARED par_ls.cpp)

target_include_directories(ParLeastSquares PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(ParLeastSquares PUBLIC Eigen3::Eigen)
