add_executable(TestInterface test_call_interface.cpp)
add_executable(TestOperator test_call_op.cpp)
add_executable(TestMinimize test_minimize.cpp)

target_compile_options(TestInterface PRIVATE -g)
target_compile_options(TestOperator PRIVATE -g)
target_compile_options(TestMinimize PRIVATE -g)

target_link_libraries(TestOperator 
    PRIVATE options warnings ParLeastSquares)

target_link_libraries(TestMinimize 
    PRIVATE options warnings ParLeastSquares)

target_link_libraries(TestInterface 
    PRIVATE options warnings ParLeastSquares)

add_test(
    NAME TestInterface
    COMMAND $<TARGET_FILE:TestInterface> ${CMAKE_SOURCE_DIR}/tests/test_cases/matrix_market)

add_test(
    NAME TestOperator
    COMMAND $<TARGET_FILE:TestOperator> ${CMAKE_SOURCE_DIR}/tests/test_cases/matrix_market)

add_test(
    NAME TestMinimize
    COMMAND $<TARGET_FILE:TestMinimize> ${CMAKE_SOURCE_DIR}/tests/test_cases/matrix_market)

# 
# add_test(
#     NAME TestImport
#     COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_import.py)
# 
# add_test(
#     NAME TestCall_Pickle
#     COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_call.py ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cases/pickle)
# 
# add_test(
#     NAME TestCall_MatrixMarket
#     COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_call.py ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cases/matrix-market)
