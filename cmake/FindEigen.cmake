
add_library(eigen_dep INTERFACE)

find_package(Eigen3 3.3 NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  target_include_directories(eigen_dep INTERFACE ${EXTERNALS_PATH}/eigen)
else()
  target_link_libraries(eigen_dep INTERFACE Eigen3::Eigen)
endif()
