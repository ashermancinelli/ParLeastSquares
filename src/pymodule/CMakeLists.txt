set(PYBIND11_CPP_STANDARD -std=c++14)
include(FindPythonInterp)

add_library(__levmar_eigen_C MODULE module.cpp)
target_include_directories(__levmar_eigen_C
  PRIVATE
  ${Python3_INCLUDE_DIRS}
  )

target_link_libraries(__levmar_eigen_C PRIVATE 
  pybind11_dep 
  eigen_dep 
  ParLeastSquares
  # ${Python3_LIBRARIES}
  )

set_target_properties(
  __levmar_eigen_C
  PROPERTIES
  PREFIX "${PYTHON_MODULE_PREFIX}"
  SUFFIX "${PYTHON_MODULE_EXTENSION}"
  )

if(APPLE)
  set_target_properties(
    __levmar_eigen_C
    PROPERTIES
    COMPILE_FLAGS "-undefined dynamic_lookup"
    LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

add_custom_command(
  TARGET __levmar_eigen_C
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    $<TARGET_FILE:__levmar_eigen_C>
    ${PROJECT_SOURCE_DIR}/levmar_eigen/$<TARGET_FILE_NAME:__levmar_eigen_C>)
