#include "geometry/boundary.h"
#include "inputs/line/line.h"

#include "gtest/gtest.h"
#include <array>
#include <iostream>

namespace boundary {

namespace geometry {

TEST(BoundaryTests, IsBoundaryCell){
  boundary::inputs::LineGeometry input;
  std::array<double, 2> lower_left = {-1, -1};
  std::array<double, 2> lower_right = {-.75, -1};
  std::array<double, 2> upper_right = {-.75, -.75};
  std::array<double, 2> upper_left = {-1, -.75};
  bool is_boundary = Boundary::IsBoundaryCell(lower_left, lower_right,
                                    upper_right, upper_left, &input);
  EXPECT_TRUE(is_boundary);

  std::array<double, 2> lower_left2 = {-1, .75};
  std::array<double, 2> lower_right2 = {-.75, .75};
  std::array<double, 2> upper_right2 = {-.75, 1};
  std::array<double, 2> upper_left2 = {-1, 1};
  bool is_boundary2 = Boundary::IsBoundaryCell(lower_left2, lower_right2,
                                    upper_right2, upper_left2, &input);
  EXPECT_FALSE(is_boundary2);

  std::array<double, 2> lower_left3 = {-.75, -1};
  std::array<double, 2> lower_right3 = {-.5, -1};
  std::array<double, 2> upper_right3 = {-.5, -.75};
  std::array<double, 2> upper_left3 = {-.75, -.75};
  bool is_boundary3 = Boundary::IsBoundaryCell(lower_left3, lower_right3,
                                    upper_right3, upper_left3, &input);
  EXPECT_FALSE(is_boundary3);
}

TEST(BoundaryTests, BoundaryCells){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  std::map<std::array<double, 2>, geo_info> boundary_cells = line_boundary.BoundaryCells();
  EXPECT_EQ(boundary_cells.size(), 8);
  std::array<double, 2> first_point = {-.875, -.875};
  EXPECT_TRUE(boundary_cells[first_point].irregular);
}

} // namespace geometry

} // namespace boundary
