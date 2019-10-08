#include "inputs/line/line.h"
#include "helpers/geometry_objects.h"

#include "gtest/gtest.h"
#include <vector>

class LineInputTest : public ::testing::Test{
  protected:
    LineInputTest();
};

TEST(LineInputTest, Boundary){
  boundary::inputs::LineGeometry line;
  EXPECT_EQ(1, line.BoundaryFunction(1)[0]);
  EXPECT_EQ(5, line.BoundaryFunction(5)[0]);
  EXPECT_EQ(-1, line.BoundaryFunction(-1)[0]);
  EXPECT_EQ(8.3, line.BoundaryFunction(8.3)[0]);
}

TEST(LineInputTest, Derivative){
  boundary::inputs::LineGeometry line;

  boundary::helpers::Point point1 = boundary::helpers::Point(1, 1);
  boundary::helpers::Point point2 = boundary::helpers::Point(3, 3);
  boundary::helpers::Point point3 = boundary::helpers::Point(4, 4);
  boundary::helpers::Point point4 = boundary::helpers::Point(6, 6);

  std::vector<int> order1{1, 0};
  std::vector<int> order2{0, 1};
  std::vector<int> order3{0, 0};
  std::vector<int> order4{4, 5};

  EXPECT_EQ(-1, line.BoundaryDerivatives(point1, order1));
  EXPECT_EQ(1, line.BoundaryDerivatives(point2, order2));
  EXPECT_EQ(0, line.BoundaryDerivatives(point3, order3));
  EXPECT_EQ(0, line.BoundaryDerivatives(point4, order4));
}

TEST(LineInputTest, Inverse){
  boundary::inputs::LineGeometry line;
  EXPECT_EQ(1, line.BoundaryInverse(1)[0]);
  EXPECT_EQ(5, line.BoundaryInverse(5)[0]);
  EXPECT_EQ(-1, line.BoundaryInverse(-1)[0]);
  EXPECT_EQ(8.3, line.BoundaryInverse(8.3)[0]);
}
