cmake_minimum_required(VERSION 3.16)
project(wiener_filter_reference CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wiener_ref reference/wiener_filter.cpp)
target_compile_options(wiener_ref PRIVATE -Wall -Wextra)
