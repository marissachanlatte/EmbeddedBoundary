#include "inputs/line/line.h"

#include "gtest/gtest.h"

class LineInputTest : public ::testing::Test{
  protected:
    LineInputTest();
};

TEST(LineInputTest, Boundary){
  boundary::inputs::LineGeometry line;
  EXPECT_EQ(1, line.BoundaryFunction(1));
  EXPECT_EQ(5, line.BoundaryFunction(5));
  EXPECT_EQ(-1, line.BoundaryFunction(-1));
  EXPECT_EQ(8.3, line.BoundaryFunction(8.3));
}

TEST(LineInputTest, Derivative){
  boundary::inputs::LineGeometry line;
  double point1[2] = {1, 1};
  double point2[2] = {3, 3};
  double point3[2] = {4, 4};
  double point4[2] = {6, 6};

  int order1[2] = {1, 0};
  int order2[2] = {0, 1};
  int order3[2] = {0, 0};
  int order4[2] = {4, 5};

  EXPECT_EQ(-1, line.BoundaryDerivatives(point1, order1));
  EXPECT_EQ(1, line.BoundaryDerivatives(point2, order2));
  EXPECT_EQ(0, line.BoundaryDerivatives(point3, order3));
  EXPECT_EQ(0, line.BoundaryDerivatives(point4, order4));
}

TEST(LineInputTest, Inverse){
  boundary::inputs::LineGeometry line;
  EXPECT_EQ(1, line.BoundaryInverse(1));
  EXPECT_EQ(5, line.BoundaryInverse(5));
  EXPECT_EQ(-1, line.BoundaryInverse(-1));
  EXPECT_EQ(8.3, line.BoundaryInverse(8.3));
}
