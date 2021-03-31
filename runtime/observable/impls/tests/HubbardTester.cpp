#include "qcor.hpp"
#include "observable_provider.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include <gtest/gtest.h>

TEST(HubbardTester, checkSimple){
    using qcor::ObservableProvider;
    using qcor::FermionOperator;
    double coulomb = 7.5;
    double chemical_ptnl = coulomb/2;
    auto op = qcor::createOperator("hubbard", {{"x-dimension", 2}, 
                                                {"y-dimension", 1}, 
                                                {"coulomb", coulomb},
                                                {"chemical-potential", chemical_ptnl},
                                                {"spinless", false}});
    std::cout<<op->toString()<<"\n";
    EXPECT_TRUE(FermionOperator(op->toString()) == FermionOperator("(-3.75, 0) 0^ 0 + \
                                                                    (7.5, 0) 0^ 0 1^ 1 + \
                                                                    (-1.0, 0) 0^ 2 + \
                                                                    (-3.75, 0) 1^ 1 + \
                                                                    (-1.0, 0) 1^ 3 + \
                                                                    (-1.0, 0) 2^ 0 + \
                                                                    (-3.75, 0) 2^ 2 + \
                                                                    (7.5, 0) 2^ 2 3^ 3 + \
                                                                    (-1.0, 0) 3^ 1 + \
                                                                    (-3.75, 0) 3^ 3"));   

}

TEST(HubbardTester, check1x3Spinless){
    int x_dim = 1;
    int y_dim = 3;
    double coulomb = 4.;
    double chem_ptnl = 0.5;
    bool spinless = true;
    auto op = qcor::createOperator("hubbard", {{"x-dimension", x_dim},
                                                {"y-dimension", y_dim},
                                                {"coulomb", coulomb},
                                                {"chemical-potential", chem_ptnl},
                                                {"spinless", spinless}});
    std::cout<<op->toString()<<"\n";
    EXPECT_TRUE(qcor::FermionOperator(op->toString()) == qcor::FermionOperator("(-0.5, 0) 0^ 0 + \
                                                                          (4.0, 0) 0^ 0 1^ 1 + \
                                                                          (-1.0, 0) 0^ 1 + \
                                                                          (-1.0, 0) 0^ 2 + \
                                                                          (-1.0, 0) 1^ 0 + \
                                                                          (-0.5, 0) 1^ 1 + \
                                                                          (4.0, 0) 1^ 1 2^ 2 + \
                                                                          (-1.0, 0) 1^ 2 + \
                                                                          (-1.0, 0) 2^ 0 + \
                                                                          (-1.0, 0) 2^ 1 + \
                                                                          (-0.5, 0) 2^ 2 + \
                                                                          (4.0, 0) 2^ 2 0^ 0"));

}

TEST(HubbardTester, check3x1Spinless){
    int x_dim = 3;
    int y_dim = 1;
    double coulomb = 4.;
    double chem_ptnl = 0.5;
    bool spinless = true;
    auto op = qcor::createOperator("hubbard", {{"x-dimension", x_dim},
                                                {"y-dimension", y_dim},
                                                {"coulomb", coulomb},
                                                {"chemical-potential", chem_ptnl},
                                                {"spinless", spinless}});
    std::cout<<op->toString()<<"\n";
    EXPECT_TRUE(qcor::FermionOperator(op->toString()) == qcor::FermionOperator("(-0.5, 0) 0^ 0 + \
                                                                               (4.0, 0) 0^ 0 1^ 1 + \
                                                                               (-1.0, 0) 0^ 1 + \
                                                                               (-1.0, 0) 0^ 2 + \
                                                                               (-1.0, 0) 1^ 0 + \
                                                                               (-0.5, 0) 1^ 1 + \
                                                                               (4.0, 0) 1^ 1 2^ 2 + \
                                                                               (-1.0, 0) 1^ 2 + \
                                                                               (-1.0, 0) 2^ 0 + \
                                                                               (-1.0, 0) 2^ 1 + \
                                                                               (-0.5, 0) 2^ 2 + \
                                                                               (4.0, 0) 2^ 2 0^ 0"));
}



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  return ret;
}
