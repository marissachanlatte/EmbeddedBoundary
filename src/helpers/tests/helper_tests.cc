#include "helpers/math_helpers.h"

#include "gtest/gtest.h"

class HelperTest : public ::testing::Test{
  protected:
    HelperTest();
};

TEST(HelperTest, Factorial){
  boundary::helpers::MathHelper helper;
  EXPECT_EQ(1, helper.Factorial(0));
  EXPECT_EQ(6, helper.Factorial(3));
  EXPECT_EQ(362880, helper.Factorial(9));
  EXPECT_EQ(39916800, helper.Factorial(11));
}

TEST(HelperTest, MultiIndexFactorial){
  boundary::helpers::MathHelper helper;
  std::vector<int> alpha{1, 2};
  EXPECT_EQ(2, helper.MultiIndexFactorial(alpha));

  std::vector<int> beta{4, 5, 7};
  EXPECT_EQ(14515200, helper.MultiIndexFactorial(beta));
}

TEST(HelperTest, MultiIndexBinomial){
  boundary::helpers::MathHelper helper;
  std::vector<int> alpha{4, 5};
  std::vector<int> beta{1, 2};
  EXPECT_EQ(40, helper.MultiIndexBinomial(alpha, beta));

  std::vector<int> gamma{2, 5, 6};
  std::vector<int> delta{1, 3, 4};
  EXPECT_EQ(300, helper.MultiIndexBinomial(gamma, delta));
}
