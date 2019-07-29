#include "helpers/math_helpers.cc"
#include "helpers/geometry_objects.h"

#include "gtest/gtest.h"

class HelperTest : public ::testing::Test{
  protected:
    HelperTest();
};

namespace boundary {

namespace helpers {

TEST(HelperTest, Factorial){
  EXPECT_EQ(1, Factorial(0));
  EXPECT_EQ(6, Factorial(3));
  EXPECT_EQ(362880, Factorial(9));
  EXPECT_EQ(39916800, Factorial(11));
}

TEST(HelperTest, MultiIndexFactorial){
  std::vector<int> alpha{1, 2};
  EXPECT_EQ(2, MultiIndexFactorial(alpha));

  std::vector<int> beta{4, 5, 7};
  EXPECT_EQ(14515200, MultiIndexFactorial(beta));
}

TEST(HelperTest, MultiIndexBinomial){
  std::vector<int> alpha{4, 5};
  std::vector<int> beta{1, 2};
  EXPECT_EQ(40, MultiIndexBinomial(alpha, beta));

  std::vector<int> gamma{2, 5, 6};
  std::vector<int> delta{1, 3, 4};
  EXPECT_EQ(300, MultiIndexBinomial(gamma, delta));
}

TEST(HelperTest, PointAddition){
  Point alpha = Point(7, 3);
  Point beta = Point(4, 8);
  Point sum = alpha + beta;

  EXPECT_EQ(11, sum.x_val);
  EXPECT_EQ(11, sum.y_val);
}

}

}
