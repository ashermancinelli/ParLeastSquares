# add_executable(test_driver driver.cpp)

# target_link_libraries(test_driver 
# PRIVATE options warnings ParLeastSquares)

# Add built python module to ld path
set(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:${CMAKE_SOURCE_DIR}/build/src/pymodule")

add_test(
    NAME TestImport
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_import.py)

add_test(
    NAME TestCall_Pickle
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_call.py ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cases/pickle)

add_test(
    NAME TestCall_MatrixMarket
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_call.py ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cases/matrix-market)
