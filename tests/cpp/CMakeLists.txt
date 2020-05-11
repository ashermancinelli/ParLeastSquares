add_library(TestClass OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestDriver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/test_utilities.cpp)
target_include_directories(TestClass PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(TestClass PRIVATE warnings options)
target_link_libraries(TestClass PUBLIC ParLeastSquares)

add_executable(TestDriver driver.cpp)
target_link_libraries(TestDriver PRIVATE TestClass warnings options)

add_test(
  NAME TestAll
  COMMAND $<TARGET_FILE:TestDriver> ${DATA_PATH})
