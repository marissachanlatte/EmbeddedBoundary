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

TEST(BoundaryTests, Normals){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  std::map<std::array<double, 2>, geo_info> boundary_cells = line_boundary.BoundaryCells();
  std::array<double, 2> first_point = {-.875, -.875};
  std::vector<std::vector<std::vector<double>>> first_normals = boundary_cells[first_point].normal_derivatives;
  EXPECT_FLOAT_EQ(first_normals[0][0][0], -0.70710678);
  EXPECT_FLOAT_EQ(first_normals[0][0][1], 0.70710678);
  EXPECT_FLOAT_EQ(first_normals[0][1][1], 0);
  std::array<double, 2> second_point = {.125, .125};
  std::vector<std::vector<std::vector<double>>> second_normals = boundary_cells[second_point].normal_derivatives;
  EXPECT_FLOAT_EQ(second_normals[0][0][0], -0.70710678);
  EXPECT_FLOAT_EQ(second_normals[0][0][1], 0.70710678);
  EXPECT_FLOAT_EQ(second_normals[1][0][1], 0);
}


TEST(BoundaryTests, VolFrac1D){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  std::map<std::array<double, 2>, geo_info> boundary_cells = line_boundary.BoundaryCells();
  std::array<double, 2> first_point = {-.875, -.875};
  std::array<double, 4> first_vol_frac = boundary_cells[first_point].vol_frac_1d;
  EXPECT_EQ(first_vol_frac[0], 0);
  EXPECT_EQ(first_vol_frac[1], 0);
  EXPECT_EQ(first_vol_frac[2], .25);
  EXPECT_EQ(first_vol_frac[3], .25);
  std::array<double, 2> second_point = {.125, .125};
  std::array<double, 4> second_vol_frac = boundary_cells[second_point].vol_frac_1d;
  EXPECT_EQ(second_vol_frac[0], 0);
  EXPECT_EQ(second_vol_frac[1], 0);
  EXPECT_EQ(second_vol_frac[2], .25);
  EXPECT_EQ(second_vol_frac[3], .25);
}

TEST(BoundaryTests, VolumeMoments){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  std::map<std::array<double, 2>, geo_info> boundary_cells = line_boundary.BoundaryCells();
  double total_area = 0;
  for (std::map<std::array<double, 2>, geo_info>::iterator it=boundary_cells.begin();
       it != boundary_cells.end(); it++){
          total_area += it->second.volume_moments[0][0];
       }
  EXPECT_FLOAT_EQ(total_area, .25);
}

} // namespace geometry

} // namespace boundary
