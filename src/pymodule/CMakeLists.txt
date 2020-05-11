set(PYBIND11_CPP_STANDARD -std=c++14)
include(FindPythonInterp)

add_library(pstep MODULE module.cpp)

target_link_libraries(pstep PRIVATE
  pybind11::module
  Eigen3::Eigen
  ParLeastSquares)

set_target_properties(
    pstep
    PROPERTIES
        PREFIX "${PYTHON_MODULE_PREFIX}"
        SUFFIX "${PYTHON_MODULE_EXTENSION}")
