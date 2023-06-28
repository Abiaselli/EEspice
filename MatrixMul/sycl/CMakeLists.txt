cmake_minimum_required(VERSION 3.20.5)


project(timebenchmark LANGUAGES CXX)

find_package(IntelSYCL REQUIRED)
find_package(benchmark REQUIRED)


add_executable(matrix_mul_sycl matrix_mul_sycl.cpp)

target_link_libraries(matrix_mul_sycl benchmark::benchmark)
add_sycl_to_target(TARGET matrix_mul_sycl SOURCES matrix_mul_sycl.cpp)
