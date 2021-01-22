#include "helpers/math_helpers.cc"
#include "helpers/geometry_objects.h"

#include "gtest/gtest.h"
#include <cmath>

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
  std::vector<double> alpha_coords{7, 3};
  Point alpha = Point(alpha_coords);
  std::vector<double> beta_coords{4, 8};
  Point beta = Point(beta_coords);
  Point sum = alpha + beta;

  EXPECT_EQ(11, sum.value(0));
  EXPECT_EQ(11, sum.value(1));
}


TEST(HelperTest, InitializePoint){
  Point test = Point();
  EXPECT_TRUE(std::isnan(test.value(0)));
}


TEST(HelperTest, NormalizeVector){
  std::vector<double> test_vector{4.5, 1, .4};
  std::vector<double> maxes{10, 10, 10};
  std::vector<double> mins{-10, -10, -10};
  std::vector<double> correct{0.725, 0.55, 0.52};
  std::vector<double> normalized = NormalizeVector(test_vector, maxes, mins);
  EXPECT_EQ(normalized[0], correct[0]);
  EXPECT_EQ(normalized[1], correct[1]);
  EXPECT_EQ(normalized[2], correct[2]);
}


TEST(HelperTest, IntegerMap){
  std::vector<double> test_point{.3, .6};
  std::vector<int> zero{0, 0};
  EXPECT_EQ(IntegerMap(test_point, 0), zero);
  std::vector<int> zero_one{0, 1};
  EXPECT_EQ(IntegerMap(test_point, 1), zero_one);
  std::vector<int> one_two{1, 2};
  EXPECT_EQ(IntegerMap(test_point, 2), one_two);
}


TEST(HelperTest, MortonKey){
  std::vector<double> zero_vector{0, 0, 0};
  std::vector<double> one_max{1, 1, 1};
  int key = MortonKey(zero_vector, 0, one_max, zero_vector);
  EXPECT_EQ(key, 1);
  std::vector<double> test_point{.6, .2};
  int key1 = MortonKey(test_point, 1, one_max, zero_vector);
  EXPECT_EQ(key1, 110);
  int key2 = MortonKey(test_point, 2, one_max, zero_vector);
  EXPECT_EQ(key2, 11000);
  int key3 = MortonKey(test_point, 3, one_max, zero_vector);
  EXPECT_EQ(key3, 1100001);
}

}

}
