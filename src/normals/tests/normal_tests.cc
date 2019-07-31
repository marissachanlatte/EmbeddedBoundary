#include "normals/normals.h"
#include "inputs/line/line.h"
#include "helpers/geometry_objects.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <iostream>

class NormalTest : public ::testing::Test{
  protected:
    NormalTest();
};

TEST(NormalTest, ComputeNormal){
  boundary::normals::Normal test_normal;
  boundary::inputs::LineGeometry line;
  boundary::helpers::Point a_point = boundary::helpers::Point(1, 1);
  double test_array[2] = {-std::sqrt(2)/2, std::sqrt(2)/2};
  double* normal = test_normal.ComputeNormal(a_point, &line);

  EXPECT_THAT(test_array[0], testing::DoubleNear(normal[0], 1e-10));
  EXPECT_THAT(test_array[1], testing::DoubleNear(normal[1], 1e-10));
}

TEST(NormalTest, NormalDerivative){
  boundary::normals::Normal test_normal;
  boundary::inputs::LineGeometry line;
  boundary::helpers::Point a_point = boundary::helpers::Point(1, 1);
  std::vector<int> p_order1{0, 0};
  double derivative1 = test_normal.NormalDerivative(p_order1, 1, a_point, &line);
  EXPECT_THAT(-std::sqrt(2)/2, testing::DoubleNear(derivative1, 1e-10));

  std::vector<int> p_order2{1, 0};
  double derivative2 = test_normal.NormalDerivative(p_order2, 1, a_point, &line);
  EXPECT_EQ(0, derivative2);
}
