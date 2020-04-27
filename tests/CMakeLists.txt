add_executable(test_driver driver.cpp)

target_link_libraries(test_driver 
    PRIVATE options warnings ParLeastSquares)
