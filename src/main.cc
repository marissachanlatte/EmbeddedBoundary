#include "inputs/line/line.h"
#include "inputs/circle/circle.h"
#include "inputs/square/square.h"
#include "helpers/geometry_objects.h"
#include "normals/normals.h"
#include "geometry/boundary.h"
#include "solvers/laplacian.h"

#include <stack>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

void writeGeometry(std::map<int, int> cell_map, std::map<std::array<double, 2>, boundary::geometry::geo_info> boundary_cells,
                   boundary::geometry::Boundary boundary, double cell_size){
  std::ofstream output;
  output.open ("../outputs/output_boundary.txt");
  for (std::map<int, int>::iterator it=cell_map.begin(); it != cell_map.end(); it++){
    std::array<double, 2> center = boundary.IDtoCenter(it->first);
    if (it->second == 0){
      // exterior
      output << center[0] << " " << center[1] << " " << 0 << std::endl;
    }
    else if (it->second == 1){
      // interior
      output << center[0] << " " << center[1] << " " << std::pow(boundary.CellSize(), 2) << std::endl;
    }
    else if (it->second == 2){
      // boundary
      output << center[0] << " " << center[1] << " " << boundary_cells[center].volume_moments[0][0] << std::endl;
    }   
  }
  
}


void writeSolution(Eigen::VectorXd solution, boundary::geometry::Boundary boundary){
  std::ofstream output;
  output.open ("../outputs/output_solution.txt");
  for (int id = 0; id < solution.size(); id++){
    std::array<double, 2> point = boundary.IDtoCenter(id);
    output << point[0] << " " << point[1] << " " << solution[id] << std::endl;
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
  std::map<std::array<double, 2>, boundary::geometry::geo_info> boundary_cells = boundary.BoundaryCells();
  std::map<int, int> cell_map = boundary.CellMap();
  double cell_size = boundary.CellSize();
  // Make laplacian
  Eigen::VectorXd solution = makeLaplacian(boundary);
  // Write to file
  writeSolution(solution, boundary);
  // writeGeometry(cell_map, boundary_cells, boundary, cell_size);
  return 0;
}


