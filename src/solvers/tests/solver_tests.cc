#include "solvers/laplacian.h"
#include "inputs/geometries/circle/circle_test.h"
#include "geometry/boundary.h"

#include "gtest/gtest.h"
#include <iostream>

// TEST(LaplaceTest, Symmetry){
//   // Read in input
//   boundary::inputs::CircleTestGeometry input;
//   // Make geometry
//   boundary::geometry::Boundary boundary = boundary::geometry::Boundary(&input);
//   std::map<int, boundary::geometry::geo_info> boundary_cells = boundary.BoundaryCells();
//   std::map<int, int> cell_map = boundary.CellMap();
//   double cell_size = boundary.InitialCellSize();
//   // Make laplacian
//   boundary::solvers::Laplacian laplacian = boundary::solvers::Laplacian(&input, boundary);
//   Eigen::VectorXd solution = laplacian.solve();
//   // Check symmetry
//   int n = std::sqrt(solution.size());
//   int depth = input.MaxSolverDepth();
//   for (int i = 0; i < n; i ++){
//     for (int j = 0; j < (i + 1); j++){
//       EXPECT_NEAR(solution[boundary.IJToGlobal(i, j, depth)], solution[boundary.IJToGlobal(j, i, depth)], 1e-5);
//       EXPECT_NEAR(solution[boundary.IJToGlobal(i, j, depth)], solution[boundary.IJToGlobal(i, n - j - 1, depth)], 1e-5);
//       EXPECT_NEAR(solution[boundary.IJToGlobal(i, j, depth)], solution[boundary.IJToGlobal(n - i - 1, j, depth)], 1e-5);
//     }
//   }
// }

