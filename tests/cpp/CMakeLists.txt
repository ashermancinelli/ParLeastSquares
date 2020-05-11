add_library(TestClass OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/tests/TestDriver.cpp)
target_include_directories(TestClass PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(TestClass PUBLIC
  ParLeastSquares
  warnings
  options)

add_executable(TestDriver driver.cpp)
target_link_libraries(TestDriver PRIVATE TestClass)

add_test(
  NAME TestAll
  COMMAND $<TARGET_FILE:TestDriver>)
