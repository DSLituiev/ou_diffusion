cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=c++0x")

# include_directories (sambamio)

add_executable(simulator main.cpp)
 					
TARGET_LINK_LIBRARIES(simulator sqlite3pp)

