set(PYBIND11_CPP_STANDARD -std=c++14)
include(FindPythonInterp)

add_library(__levmar_eigen_C MODULE module.cpp)

message(STATUS "Linking against ${Python3_LIBRARIES}")
target_link_libraries(__levmar_eigen_C PRIVATE 
  pybind11_dep 
  eigen_dep 
  ParLeastSquares
  ${Python3_LIBRARIES})

set_target_properties(
  __levmar_eigen_C
  PROPERTIES
  PREFIX "${PYTHON_MODULE_PREFIX}"
  SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS __levmar_eigen_C DESTINATION ${Python3_SITELIB})
