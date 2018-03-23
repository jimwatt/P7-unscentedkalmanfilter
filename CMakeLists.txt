project(ExtendedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp src/tools.cpp src/FusionEKF.cpp src/FusionUKF.cpp src/extended_kalman_filter.cpp src/unscented_kalman_filter.cpp src/utility.cpp src/tools.h src/FusionEKF.h src/FusionUKF.h src/FusionKF.h src/extended_kalman_filter.h src/unscented_kalman_filter.h src/utility.h)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1*/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 


add_executable(UnscentedKF ${sources})

target_link_libraries(UnscentedKF z ssl uv uWS)
