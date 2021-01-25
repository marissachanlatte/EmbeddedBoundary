#include "solvers/laplacian.h"
#include "inputs/circle/circle_test.h"
#include "geometry/boundary.h"

#include "gtest/gtest.h"
#include <iostream>

TEST(LaplaceTest, Symmetry){
  // Read in input
  boundary::inputs::CircleTestGeometry input;
  // Make geometry
  boundary::geometry::Boundary boundary = boundary::geometry::Boundary(&input);
  std::map<std::array<double, 2>, boundary::geometry::geo_info> boundary_cells = boundary.BoundaryCells();
  std::map<int, int> cell_map = boundary.CellMap();
  double cell_size = boundary.InitialCellSize();
  // Make laplacian
  boundary::solvers::Laplacian laplacian = boundary::solvers::Laplacian(boundary);
  Eigen::VectorXd solution = laplacian.solve();
  // Check symmetry
  int n = std::sqrt(solution.size());
  for (int i = 0; i < n; i ++){
    for (int j = 0; j < (i + 1); j++){
      EXPECT_NEAR(solution[boundary.IJToGlobal(i, j)], solution[boundary.IJToGlobal(j, i)], 1e-5);
      EXPECT_NEAR(solution[boundary.IJToGlobal(i, j)], solution[boundary.IJToGlobal(i, n - j - 1)], 1e-5);
      EXPECT_NEAR(solution[boundary.IJToGlobal(i, j)], solution[boundary.IJToGlobal(n - i - 1, j)], 1e-5);
    }
  }
}

