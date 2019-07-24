#include "normals/normals.h"
#include "inputs/line/line.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>


class NormalTest : public ::testing::Test{
  protected:
    NormalTest();
};

TEST(NormalTest, ComputeNormal){
  boundary::normals::Normal test_normal;
  boundary::inputs::LineGeometry line;
  double x_value[2] = {1, 1};
  double test_array[2] = {-std::sqrt(2)/2, std::sqrt(2)/2};
  double* normal = test_normal.ComputeNormal(x_value, &line);

  EXPECT_THAT(test_array[0], testing::DoubleNear(normal[0], 1e-10));
  EXPECT_THAT(test_array[1], testing::DoubleNear(normal[1], 1e-10));
}
