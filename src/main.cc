#include "inputs/line/line.h"
#include "inputs/circle/circle.h"
#include "inputs/square/square.h"
#include "helpers/geometry_objects.h"
#include "normals/normals.h"
#include "geometry/boundary.h"
#include "solvers/laplacian.h"
#include "helpers/math_helpers.cc"

#include <stack>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

// TODO: Rewrite for adaptive mesh
void writeGeometry(std::map<int, int> cell_map, std::map<int, boundary::geometry::geo_info> boundary_cells,
                   boundary::geometry::Boundary boundary){
  std::ofstream output;
  output.open ("../outputs/output_boundary.txt");
  int depth = boundary.MaxSolverDepth();
  int num_x = std::pow(2, depth);
  double cell_size = boundary.InitialCellSize()/num_x;
  std::vector<double> maxes = boundary.Maxes();
  std::vector<double> mins = boundary.Mins();
  // Iterate through all cells
  for (int i = 0; i < num_x; i++){ // y-index
    for (int j = 0; j < num_x; j++){ // x-index
      int matrix_id = boundary.IJToGlobal(i, j, depth);
      // get cell center
      std::vector<double> center{boundary.XMin() + j*cell_size + cell_size/2, 
                                 boundary.YMin() + i*cell_size + cell_size/2};
      int key = boundary::helpers::MortonKey(center, depth, maxes, mins);
      if (cell_map[key] == 0){
        // exterior
        output << center[0] << " " << center[1] << " " << 0 << std::endl;
      }
      else if (cell_map[key] == 1){
        // interior
        output << center[0] << " " << center[1] << " " << std::pow(cell_size, 2) << std::endl;
      }
      else if (cell_map[key] == 2){
        // boundary
        output << center[0] << " " << center[1] << " " << boundary_cells[key].volume_moments[0][0] << std::endl;
      }   
    }
  }
}


// TODO: Rewrite for Adaptive Mesh
void writeSolution(Eigen::VectorXd solution, boundary::geometry::Boundary boundary){
  std::ofstream output;
  output.open ("../outputs/output_solution.txt");
  // for (int id = 0; id < solution.size(); id++){
  //   std::vector<double> point = boundary.IDtoCenter(id);
  //   output << point[0] << " " << point[1] << " " << solution[id] << std::endl;
  // }
  int depth = boundary.MaxSolverDepth();
  int num_x = std::pow(2, depth);
  double cell_size = boundary.InitialCellSize()/num_x;
  // Iterate through all cells
  for (int i = 0; i < num_x; i++){ // y-index
    for (int j = 0; j < num_x; j++){ // x-index
      int matrix_id = boundary.IJToGlobal(i, j, depth);
      // get cell center
      std::vector<double> cell_center{boundary.XMin() + j*cell_size + cell_size/2, 
                                      boundary.YMin() + i*cell_size + cell_size/2};
      output << cell_center[0] << " " << cell_center[1] << " " << solution[matrix_id] << std::endl;
    }
  }
  output.close();
}


Eigen::VectorXd makeLaplacian(boundary::geometry::Boundary boundary){
  boundary::solvers::Laplacian laplacian = boundary::solvers::Laplacian(boundary);
  Eigen::VectorXd solution = laplacian.solve();
  return solution;
}


int main(){
  // Read in input
  boundary::inputs::SquareGeometry input;
  // Make geometry
  boundary::geometry::Boundary boundary = boundary::geometry::Boundary(&input);
  std::map<int, boundary::geometry::geo_info> boundary_cells = boundary.BoundaryCells();
  std::map<int, int> cell_map = boundary.CellMap();
  // Make laplacian
  Eigen::VectorXd solution = makeLaplacian(boundary);
  // Write to file
  writeSolution(solution, boundary);
  writeGeometry(cell_map, boundary_cells, boundary);
  return 0;
}


