#include "geometry/boundary.h"
#include "inputs/line/line.h"
#include "inputs/circle/circle.h"

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
  boundary::inputs::CircleGeometry circle;
  Boundary circle_boundary = Boundary(&circle);
  std::map<std::array<double, 2>, geo_info> circle_cells = circle_boundary.BoundaryCells();
  std::array<double, 2> first_circle_point = {-.625, -.875};
  std::array<double, 4> first_circle_frac = circle_cells[first_circle_point].vol_frac_1d;
  EXPECT_FLOAT_EQ(first_circle_frac[0], 0);
  EXPECT_FLOAT_EQ(first_circle_frac[1], 0.16143783);
  EXPECT_FLOAT_EQ(first_circle_frac[2], 0.1160254);
  EXPECT_FLOAT_EQ(first_circle_frac[3], 0);
  std::array<double, 2> second_circle_point = {.625, -.875};
  std::array<double, 4> second_circle_frac = circle_cells[second_circle_point].vol_frac_1d;
  EXPECT_FLOAT_EQ(second_circle_frac[0], 0.1160254);
  EXPECT_FLOAT_EQ(second_circle_frac[1], 0.16143783);
  EXPECT_FLOAT_EQ(second_circle_frac[2], 0);
  EXPECT_FLOAT_EQ(second_circle_frac[3], 0);

}

TEST(BoundaryTests, WhichValue){
  std::vector<double> values;
  values.push_back(-1);
  values.push_back(1);
  EXPECT_EQ(1, Boundary::WhichValue(values, .5, 1.75));
  EXPECT_EQ(-1, Boundary::WhichValue(values, -.5, -1.75));

}

TEST(BoundaryTests, NormalDerivatives){
  boundary::inputs::CircleGeometry circle;
  Boundary circle_boundary = Boundary(&circle);
  std::map<std::array<double, 2>, geo_info> boundary_cells = circle_boundary.BoundaryCells();
  std::array<double, 2> first_point = {.625, -.875};
  EXPECT_FLOAT_EQ(0.58123819, boundary_cells[first_point].normal_derivatives[0][0][0]);
  EXPECT_FLOAT_EQ(-0.81373347, boundary_cells[first_point].normal_derivatives[0][0][1]);
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

  boundary::inputs::CircleGeometry circle_input;
  Boundary circle_boundary = Boundary(&circle_input);
  std::map<std::array<double, 2>, geo_info> circle_boundary_cells = circle_boundary.BoundaryCells();
  double total_circle = 2;
  for (std::map<std::array<double, 2>, geo_info>::iterator it=circle_boundary_cells.begin();
       it != circle_boundary_cells.end(); it++){
          total_circle += it->second.volume_moments[0][0];
       }
  std::cout << total_circle << std::endl;
  EXPECT_NEAR(total_circle, 3.14, 5e-2);

}


TEST(BoundaryTests, DomainLimits){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  EXPECT_FLOAT_EQ(line_boundary.XMax(), 1.0);
  EXPECT_FLOAT_EQ(line_boundary.XMin(), -1.0);
  EXPECT_FLOAT_EQ(line_boundary.YMax(), 1.0);
  EXPECT_FLOAT_EQ(line_boundary.YMin(), -1.0);
}

} // namespace geometry

} // namespace boundary
