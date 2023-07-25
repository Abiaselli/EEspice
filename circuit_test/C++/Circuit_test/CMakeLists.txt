cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

project (Simulation)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_compile_options(-O2 -std=c++17 -lopenblas)

# add_executable(Transient_code Transient_code.cpp Transient_code.h)
add_executable(RC_Test RC_test.cpp Transient_code.h)
# add_executable(MOSFET MOSFET_IV.cpp Transient_code.h)

# target_link_libraries(Transient_code ${ARMADILLO_LIBRARIES})
target_link_libraries(RC_Test ${ARMADILLO_LIBRARIES})
# target_link_libraries(MOSFET ${ARMADILLO_LIBRARIES})

