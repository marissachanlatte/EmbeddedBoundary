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
