#include "geometry/boundary.h"
#include "inputs/geometries/line/line.h"
#include "inputs/geometries/circle/circle_test.h"
#include "inputs/geometries/ellipse/ellipse.h"

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
  std::map<double, geo_info> boundary_cells = line_boundary.BoundaryCells();
  EXPECT_EQ(boundary_cells.size(), 15);
  EXPECT_TRUE(boundary_cells[1].irregular);
  EXPECT_TRUE(boundary_cells[111].irregular);
  EXPECT_TRUE(boundary_cells[1000000].irregular);
  EXPECT_TRUE(boundary_cells[1000011].irregular);
  EXPECT_TRUE(boundary_cells[11100].irregular);
}

TEST(BoundaryTests, Normals){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  std::map<double, geo_info> boundary_cells = line_boundary.BoundaryCells();
  std::vector<std::vector<std::vector<double>>> first_normals = boundary_cells[1000000].normal_derivatives;
  EXPECT_FLOAT_EQ(first_normals[0][0][0], -0.70710678);
  EXPECT_FLOAT_EQ(first_normals[0][0][1], 0.70710678);
  EXPECT_FLOAT_EQ(first_normals[0][1][1], 0);
  std::vector<std::vector<std::vector<double>>> second_normals = boundary_cells[111].normal_derivatives;
  EXPECT_FLOAT_EQ(second_normals[0][0][0], -0.70710678);
  EXPECT_FLOAT_EQ(second_normals[0][0][1], 0.70710678);
  EXPECT_FLOAT_EQ(second_normals[1][0][1], 0);
}


TEST(BoundaryTests, VolFrac1D){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  std::map<double, geo_info> boundary_cells = line_boundary.BoundaryCells();
  std::array<double, 4> first_vol_frac = boundary_cells[1000000].vol_frac_1d;
  EXPECT_EQ(first_vol_frac[0], 0);
  EXPECT_EQ(first_vol_frac[1], 0);
  EXPECT_EQ(first_vol_frac[2], .25);
  EXPECT_EQ(first_vol_frac[3], .25);
  std::array<double, 4> parent_vol_frac = boundary_cells[10000].vol_frac_1d;
  EXPECT_EQ(parent_vol_frac[0], 0);
  EXPECT_EQ(parent_vol_frac[1], 0);
  EXPECT_EQ(parent_vol_frac[2], .5);
  EXPECT_EQ(parent_vol_frac[3], .5);
  std::array<double, 4> second_vol_frac = boundary_cells[1110000].vol_frac_1d;
  EXPECT_EQ(second_vol_frac[0], 0);
  EXPECT_EQ(second_vol_frac[1], 0);
  EXPECT_EQ(second_vol_frac[2], .25);
  EXPECT_EQ(second_vol_frac[3], .25);
}
TEST(BoundaryTests, VolFrac1dCircle){
  boundary::inputs::CircleTestGeometry circle;
  Boundary circle_boundary = Boundary(&circle);
  std::map<double, geo_info> circle_cells = circle_boundary.BoundaryCells();
  std::array<double, 4> first_circle_frac = circle_cells[1000010].vol_frac_1d;
  EXPECT_FLOAT_EQ(first_circle_frac[0], 0);
  EXPECT_FLOAT_EQ(first_circle_frac[1], 0.16143783);
  EXPECT_FLOAT_EQ(first_circle_frac[2], 0.1160254);
  EXPECT_FLOAT_EQ(first_circle_frac[3], 0);
  std::array<double, 4> second_circle_frac = circle_cells[1101000].vol_frac_1d;
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
  boundary::inputs::CircleTestGeometry circle;
  Boundary circle_boundary = Boundary(&circle);
  std::map<double, geo_info> boundary_cells = circle_boundary.BoundaryCells();
  EXPECT_FLOAT_EQ(0.58123819, boundary_cells[1101000].normal_derivatives[0][0][0]);
  EXPECT_FLOAT_EQ(-0.81373347, boundary_cells[1101000].normal_derivatives[0][0][1]);
}

TEST(BoundaryTests, VolumeMoments){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  std::map<double, geo_info> boundary_cells = line_boundary.BoundaryCells();
  EXPECT_FLOAT_EQ(boundary_cells[1].volume_moments[0][0], 2);

  boundary::inputs::CircleTestGeometry circle_input;
  Boundary circle_boundary = Boundary(&circle_input);
  std::map<double, geo_info> circle_boundary_cells = circle_boundary.BoundaryCells();
  EXPECT_NEAR(circle_boundary_cells[1].volume_moments[0][0], 3.14, 5e-2);

  boundary::inputs::EllipseGeometry ellipse_input;
  Boundary ellipse_boundary = Boundary(&ellipse_input);
  std::map<double, geo_info> ellipse_boundary_cells = ellipse_boundary.BoundaryCells();
  EXPECT_NEAR(ellipse_boundary_cells[1].volume_moments[0][0], std::sqrt(1.0/2)*3.14, 5e-2);
}

TEST(BoundaryTests, BoundaryMoments){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  std::map<double, geo_info> boundary_cells = line_boundary.BoundaryCells();
  EXPECT_FLOAT_EQ(boundary_cells[1].boundary_moments[0][0], std::sqrt(8));

  boundary::inputs::CircleTestGeometry circle_input;
  Boundary circle_boundary = Boundary(&circle_input);
  std::map<double, geo_info> circle_boundary_cells = circle_boundary.BoundaryCells();
  EXPECT_NEAR(circle_boundary_cells[1].boundary_moments[0][0], 2*3.14, 5e-2);
}

TEST(BoundaryTests, DomainLimits){
  boundary::inputs::LineGeometry input;
  Boundary line_boundary = Boundary(&input);
  EXPECT_FLOAT_EQ(line_boundary.XMax(), 1.0);
  EXPECT_FLOAT_EQ(line_boundary.XMin(), -1.0);
  EXPECT_FLOAT_EQ(line_boundary.YMax(), 1.0);
  EXPECT_FLOAT_EQ(line_boundary.YMin(), -1.0);
}


TEST(BoundaryTests, NeighborCell){
  boundary::inputs::CircleTestGeometry input;
  Boundary circle_boundary = Boundary(&input);
  // left
  EXPECT_EQ(circle_boundary.NeighborCell(2, 2, 0)[0], 2);
  EXPECT_EQ(circle_boundary.NeighborCell(2, 2, 0)[1], 1);
  // up
  EXPECT_EQ(circle_boundary.NeighborCell(2, 2, 1)[0], 3);
  EXPECT_EQ(circle_boundary.NeighborCell(2, 2, 1)[1], 2);
  // right
  EXPECT_EQ(circle_boundary.NeighborCell(2, 2, 2)[0], 2);
  EXPECT_EQ(circle_boundary.NeighborCell(2, 2, 2)[1], 3);
  // down
  EXPECT_EQ(circle_boundary.NeighborCell(2, 2, 3)[0], 1);
  EXPECT_EQ(circle_boundary.NeighborCell(2, 2, 3)[1], 2);
}


TEST(BoundaryTests, IJToGlobal){
  boundary::inputs::CircleTestGeometry input;
  Boundary circle_boundary = Boundary(&input);
  EXPECT_EQ(circle_boundary.IJToGlobal(2, 2, 3), 18);
  EXPECT_EQ(circle_boundary.IJToGlobal(3, 2, 3), 26);
  EXPECT_EQ(circle_boundary.IJToGlobal(2, 3, 3), 19);
}


TEST(BoundaryTests, InterpolationPair){
  boundary::inputs::CircleTestGeometry input;
  Boundary circle_boundary = Boundary(&input);
  std::array<std::array<int, 2>, 2> pair1 = circle_boundary.InterpolationPair(2, 2, -1, -1, 1);
  EXPECT_EQ(2, pair1[0][0]);
  EXPECT_EQ(3, pair1[0][1]);
  EXPECT_EQ(3, pair1[1][0]);
  EXPECT_EQ(3, pair1[1][1]);
  std::array<std::array<int, 2>, 2> pair2 = circle_boundary.InterpolationPair(2, 2, -1, -1, 2);
  EXPECT_EQ(3, pair2[0][0]);
  EXPECT_EQ(2, pair2[0][1]);
  EXPECT_EQ(3, pair2[1][0]);
  EXPECT_EQ(3, pair2[1][1]);
  std::array<std::array<int, 2>, 2> pair3 = circle_boundary.InterpolationPair(2, 5, 1, -1, 0);
  EXPECT_EQ(3, pair3[0][0]);
  EXPECT_EQ(5, pair3[0][1]);
  EXPECT_EQ(3, pair3[1][0]);
  EXPECT_EQ(4, pair3[1][1]);
}
} // namespace geometry

} // namespace boundary
