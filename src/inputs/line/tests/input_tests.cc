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
  EXPECT_EQ(1, line.BoundaryDerivatives(1, 1));
  EXPECT_EQ(1, line.BoundaryDerivatives(3, 1));
  EXPECT_EQ(4, line.BoundaryDerivatives(4, 0));
  EXPECT_EQ(0, line.BoundaryDerivatives(6, 4));
}

TEST(LineInputTest, Inverse){
  boundary::inputs::LineGeometry line;
  EXPECT_EQ(1, line.BoundaryInverse(1));
  EXPECT_EQ(5, line.BoundaryInverse(5));
  EXPECT_EQ(-1, line.BoundaryInverse(-1));
  EXPECT_EQ(8.3, line.BoundaryInverse(8.3));
}
