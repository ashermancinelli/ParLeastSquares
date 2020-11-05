add_library(TestClass OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestDriver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/test_utilities.cpp)
target_include_directories(TestClass
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  )

target_link_libraries(TestClass PRIVATE warnings options)
target_link_libraries(TestClass PUBLIC ParLeastSquares)

add_executable(LeastSquares least_squares.cpp)
target_link_libraries(LeastSquares PRIVATE TestClass warnings options)

foreach(maxfevs RANGE 1000 5000 1000)
  foreach(xtol_power RANGE 5 10 20)
    foreach(method numerical analytical)
      add_test(
        NAME "Test_Maxfevs_${maxfevs}_XTol_1e-${xtol_power}_${method}"
        COMMAND $<TARGET_FILE:LeastSquares> ${maxfevs} 1e-${xtol_power} ${method} ${DATA_PATH})
    endforeach()
  endforeach()
endforeach()
