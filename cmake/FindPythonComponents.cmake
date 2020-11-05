
set(CMAKE_FIND_FRAMEWORK NEVER)

if(NOT "$ENV{Python3_INCLUDE_DIRS}" STREQUAL "")
  set(Python3_INCLUDE_DIRS $ENV{Python3_INCLUDE_DIRS})
else()
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
endif()

message(STATUS "Python include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python3_LIBRARIES}")
message(STATUS "Python site-packages directory: ${Python3_SITELIB}")

include_directories(${Python3_INCLUDE_DIRS})
