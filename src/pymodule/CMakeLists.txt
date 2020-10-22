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
  SUFFIX "${PYTHON_MODULE_EXTENSION}"
  INSTALL_RPATH ${Python3_SITELIB})
add_custom_command(
  TARGET __levmar_eigen_C
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    $<TARGET_FILE:__levmar_eigen_C>
    ${PROJECT_SOURCE_DIR}/levmar_eigen/$<TARGET_FILE_NAME:__levmar_eigen_C>)

install(TARGETS __levmar_eigen_C DESTINATION ${Python3_SITELIB})
