configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/defs.hpp.in 
  ${CMAKE_CURRENT_BINARY_DIR}/defs.hpp)

add_library(ParLeastSquares SHARED ParLeastSquares.cpp)

target_include_directories(ParLeastSquares PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(ParLeastSquares PUBLIC Eigen3::Eigen)
