
add_library(pybind11_dep INTERFACE)

find_package(pybind11)
if(NOT TARGET pybind11::module)
  target_include_directories(pybind11_dep INTERFACE
    ${EXTERNALS_PATH}/pybind11/include)
  include(${EXTERNALS_PATH}/pybind11/tools/FindPythonLibsNew.cmake)
else()
  target_link_libraries(pybind11_dep INTERFACE pybind11::module)
endif()

target_include_directories(pybind11_dep INTERFACE ${Python_INCLUDE_DIRS})
