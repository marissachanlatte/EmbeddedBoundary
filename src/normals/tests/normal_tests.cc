#include "normals/normals.h"
#include "inputs/geometries/line/line.h"
#include "inputs/geometries/circle/circle.h"
#include "helpers/geometry_objects.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <iostream>
#include <array>

class NormalTest : public ::testing::Test{
  protected:
    NormalTest();
};


TEST(NormalTest, ComputeNormal){
  boundary::inputs::LineGeometry line;
  std::vector<double> a_point_coords{1, 1};
  boundary::helpers::Point a_point = boundary::helpers::Point(a_point_coords);
  double test_array[2] = {-std::sqrt(2)/2, std::sqrt(2)/2};
  std::array<double, 2> normal = boundary::normals::Normal::ComputeNormal(a_point, &line);

  EXPECT_THAT(test_array[0], testing::DoubleNear(normal[0], 1e-10));
  EXPECT_THAT(test_array[1], testing::DoubleNear(normal[1], 1e-10));

  boundary::inputs::CircleGeometry circle;
  std::vector<double> second_point_coords{1, 0};
  boundary::helpers::Point second_point = boundary::helpers::Point(second_point_coords);
  std::array<double, 2> circle_normal = boundary::normals::Normal::ComputeNormal(second_point, &circle);
  EXPECT_FLOAT_EQ(1, circle_normal[0]);
  EXPECT_FLOAT_EQ(0, circle_normal[1]);
}

TEST(NormalTest, NormalDerivative){
  boundary::inputs::LineGeometry line;
  std::vector<double> a_point_coords{1, 1};
  boundary::helpers::Point a_point = boundary::helpers::Point(a_point_coords);
  std::vector<int> p_order1{0, 0};
  double derivative1 = boundary::normals::Normal::NormalDerivative(p_order1, 1, a_point, &line);
  EXPECT_THAT(-std::sqrt(2)/2, testing::DoubleNear(derivative1, 1e-10));

  std::vector<int> p_order2{1, 0};
  double derivative2 = boundary::normals::Normal::NormalDerivative(p_order2, 1, a_point, &line);
  EXPECT_EQ(0, derivative2);
}
