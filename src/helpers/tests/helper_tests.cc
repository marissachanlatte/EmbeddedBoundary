#include "helpers/math_helpers.cc"

#include "gtest/gtest.h"

class HelperTest : public ::testing::Test{
  protected:
    HelperTest();
};

TEST(HelperTest, Factorial){
  EXPECT_EQ(1, boundary::helpers::Factorial(0));
  EXPECT_EQ(6, boundary::helpers::Factorial(3));
  EXPECT_EQ(362880, boundary::helpers::Factorial(9));
  EXPECT_EQ(39916800, boundary::helpers::Factorial(11));
}

TEST(HelperTest, MultiIndexFactorial){
  std::vector<int> alpha{1, 2};
  EXPECT_EQ(2, boundary::helpers::MultiIndexFactorial(alpha));

  std::vector<int> beta{4, 5, 7};
  EXPECT_EQ(14515200, boundary::helpers::MultiIndexFactorial(beta));
}

TEST(HelperTest, MultiIndexBinomial){
  std::vector<int> alpha{4, 5};
  std::vector<int> beta{1, 2};
  EXPECT_EQ(40, boundary::helpers::MultiIndexBinomial(alpha, beta));

  std::vector<int> gamma{2, 5, 6};
  std::vector<int> delta{1, 3, 4};
  EXPECT_EQ(300, boundary::helpers::MultiIndexBinomial(gamma, delta));
}
