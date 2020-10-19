
set(CMAKE_FIND_FRAMEWORK NEVER)

find_package(Python3 3.4 COMPONENTS Interpreter Development REQUIRED)

message(STATUS "Python include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python3_LIBRARIES}")
message(STATUS "Python site-packages directory: ${Python3_SITELIB}")

include_directories(${Python3_INCLUDE_DIRS})
