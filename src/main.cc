#include "inputs/geometries/line/line.h"
#include "inputs/geometries/circle/circle.h"
#include "inputs/geometries/ellipse/ellipse.h"
#include "inputs/geometries/ellipse/ellipse_flip.h"
#include "inputs/geometries/circle/circle_test.h"
#include "inputs/geometries/circle/circle_shift.h"
#include "inputs/geometries/square/square.h"
#include "inputs/equations/laplace_neumann.h"
#include "inputs/equations/cos_sin_neumann.h"
#include "helpers/geometry_objects.h"
#include "normals/normals.h"
#include "geometry/boundary.h"
#include "solvers/laplacian.h"
#include "solvers/laplace_operator.h"
#include "helpers/math_helpers.h"

#include <stack>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

void writeGeometry(std::map<int, int> cell_map, std::map<int, boundary::geometry::geo_info> boundary_cells,
                   boundary::geometry::Boundary boundary){
  // Points map
  std::map<std::vector<double>, int> points;
  // Cell vector
  std::vector<std::vector<int>> cells;
  // Iterate through all boundary cells
  std::map<int, boundary::geometry::geo_info>::iterator it;
  int point_count = 0;
  // Points file
  std::ofstream points_out;
  points_out.open ("../outputs/points.txt");
  for (it = boundary_cells.begin(); it != boundary_cells.end(); it++){
    std::vector<double> center = it->second.cell_center;
    double cell_size = it->second.cell_size;
    // Quad vector
    std::vector<int> quad;
    // Update points and cells
    int pos_neg[4][2] = {{1, 1}, {1, -1}, {-1, -1}, {-1, 1}};
    for (int i = 0; i < 4; i++){
      std::vector<double> corner{center[0] + pos_neg[i][0]*cell_size/2, 
                                 center[1] + pos_neg[i][1]*cell_size/2};
      if (!points.count(corner)){
        points.insert({corner, point_count});
        point_count += 1;
        // Write to file
        points_out << corner[0] << " " << corner[1] << std::endl;
      }
      quad.push_back(points[corner]);
    }
    // Add quad to cells
    cells.push_back(quad);
  }
  // Close points file
  points_out.close();

  // Write cell file
  std::ofstream cell_out;
  cell_out.open ("../outputs/cell.txt");
  for (int i = 0; i < cells.size(); i++){
    cell_out << cells[i][0] << " " << cells[i][1] << " " << cells[i][2] << " " << cells[i][3] << std::endl;
  }
  cell_out.close();
}


// TODO: Rewrite for Adaptive Mesh
void writeSolution(Eigen::VectorXd solution, boundary::geometry::Boundary boundary){
  std::ofstream output;
  output.open ("../outputs/output_solution.txt");
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

/// Checks input for suitability
void checkInput(boundary::inputs::GeometryInputBase* input){
  // Check domain is square
  if ((input->XMax() - input->XMin()) != (input->YMax() - input->YMin())){
    throw "Input Error: Domain must be square.";
  }
  if (input->MaxSolverDepth() > input->MaxDepth()){
    throw "Input Error: MaxSolverDepth can't be greater than MaxDepth.";
  }
}

int main(){
  // Read in input
  boundary::inputs::EllipseFlipGeometry geometry_input;

  try {
    checkInput(&geometry_input);
  }
  catch (const char* msg){
    std::cerr << msg << std::endl;
  }

  // Make geometry
  boundary::geometry::Boundary boundary = boundary::geometry::Boundary(&geometry_input);
  boundary::solvers::LaplaceOperator laplace_operator = boundary::solvers::LaplaceOperator(boundary);
  return 0;
}


