cmake_minimum_required(VERSION 3.20.5)


project(timebenchmark LANGUAGES CXX)

find_package(IntelSYCL REQUIRED)
find_package(benchmark REQUIRED)
find_package(MKL CONFIG REQUIRED)


add_executable(matrix_mul_onemkl matrix_mul_onemkl.cpp)

target_link_libraries(matrix_mul_onemkl benchmark::benchmark)
add_sycl_to_target(TARGET matrix_mul_onemkl SOURCES matrix_mul_onemkl.cpp)

target_compile_options(matrix_mul_onemkl PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(matrix_mul_onemkl PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(matrix_mul_onemkl $<LINK_ONLY:MKL::MKL>)
