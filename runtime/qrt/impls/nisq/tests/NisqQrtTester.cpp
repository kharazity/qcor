#include "qcor.hpp"
#include <gtest/gtest.h>
#include "xacc_service.hpp"

TEST(NisqQrtTester, checkExpInst) {
  ::quantum::initialize("qpp", "empty");
  const std::string obs_str =
      "5.907 - 2.1433 X0X1 - 2.1433 Y0Y1 + .21829 Z0 - 6.125 Z1";
  auto observable = qcor::createObservable(obs_str);
  auto qreg = qalloc(2);
  qreg.setName("q");
  ::quantum::exp(qreg, 1.0, observable);
  std::cout << "HOWDY\n"
            << ::quantum::qrt_impl->get_current_program()->toString() << "\n";
  // Get the XACC reference implementation
  auto exp = std::dynamic_pointer_cast<xacc::quantum::Circuit>(
      xacc::getService<xacc::Instruction>("exp_i_theta"));
  EXPECT_TRUE(exp->expand({{ "pauli", obs_str}}));
  auto evaled = exp->operator()({0.5});
  std::cout << "HELLO\n" << evaled->toString() << "\n";
  EXPECT_EQ(evaled->toString(), ::quantum::qrt_impl->get_current_program()->toString());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  return ret;
}
