cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

project (Simulation)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_compile_options(-O2 -std=c++17 -lopenblas)

add_executable(Oscillator Transient_code.cpp Transient_code.h)
add_executable(RC_Test RC_test.cpp Transient_code.h)
add_executable(MOSFET MOSFET_IV.cpp Transient_code.h)
add_executable(Diodo_Trans Diodo_circuit.cpp Transient_code.h)
add_executable(Diodo_IV Diodo_IV.cpp Transient_code.h)
add_executable(Voltage_doubler Voltage_doubler.cpp Transient_code.h)



target_link_libraries(Oscillator ${ARMADILLO_LIBRARIES})
target_link_libraries(RC_Test ${ARMADILLO_LIBRARIES})
target_link_libraries(MOSFET ${ARMADILLO_LIBRARIES})
target_link_libraries(Diodo_Trans ${ARMADILLO_LIBRARIES})
target_link_libraries(Diodo_IV ${ARMADILLO_LIBRARIES})
target_link_libraries(Voltage_doubler ${ARMADILLO_LIBRARIES})


