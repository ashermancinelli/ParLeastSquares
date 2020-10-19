configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/defs.hpp.in 
  ${CMAKE_CURRENT_BINARY_DIR}/defs.hpp)

add_library(ParLeastSquares STATIC
  src/interface.cpp
  src/LMFunctor.cpp
  src/LMFunctorAnalytical.cpp
  src/utilities.cpp)

target_include_directories(ParLeastSquares PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(ParLeastSquares
  PRIVATE warnings options
  PUBLIC eigen_dep)

install(TARGETS ParLeastSquares DESTINATION ${Python3_SITELIB})
