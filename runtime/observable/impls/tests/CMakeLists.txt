link_directories(${XACC_ROOT}/lib)

add_executable(HubbardTester HubbardTester.cpp)
add_test(NAME qcor_HubbardTester COMMAND HubbardTester)
target_include_directories(HubbardTester PRIVATE ../ ../../ ../../../ ${XACC_ROOT}/include/gtest)
target_link_libraries(HubbardTester ${XACC_TEST_LIBRARIES} xacc::xacc xacc::quantum_gate qcor)
