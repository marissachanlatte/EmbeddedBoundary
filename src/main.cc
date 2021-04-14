#include "inputs/geometries/line/line.h"
#include "inputs/geometries/circle/circle.h"
#include "inputs/geometries/ellipse/ellipse.h"
#include "inputs/geometries/ellipse/ellipse_flip.h"
#include "inputs/geometries/circle/circle_test.h"
#include "inputs/geometries/square/square.h"
#include "inputs/equations/laplace_neumann.h"
#include "inputs/equations/cos_sin_neumann.h"
#include "helpers/geometry_objects.h"
#include "normals/normals.h"
#include "geometry/boundary.h"
#include "solvers/laplacian.h"
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

/// Function Phi for Testing
double phi(std::vector<double> point){
  return 1.0/16*std::pow((std::pow(point[0], 2) + std::pow(point[1], 2)), 2);
}

/// BC for Phi for Testing
double neumannCondition(std::vector<double> point){
  double r = std::sqrt(std::pow(point[0], 2) + std::pow(point[1], 2));
  return 1.0/4*std::pow(r, 3);
  // return 1.0/4;
}


/// Outputs laplacian, volume fractions, boundary length, cell centers, boundary flag, and cell size as CSV file
void operatorTesting(boundary::geometry::Boundary geometry){
  int depth_ = geometry.MaxSolverDepth();
  int num_x_ = std::pow(2, depth_);
  // Loop through all cells at given depth
  double cell_size = geometry.InitialCellSize()/num_x_;
  // Get cell map
  std::map<double, int> cell_map = geometry.CellMap();
  // Get geometry information
  std::map<double, boundary::geometry::geo_info> geometry_info = geometry.BoundaryCells();
  // Output file
  std::ofstream laplace_out;
  laplace_out.open ("../outputs/laplace_out.txt");
  // Headers
  laplace_out << "Covered ID,Cell Size,CenterX,CenterY,Boundary Length,Boundary NormalX, Boundary NormalY,Volume Fraction,FirstX,FirstY,Laplacian" << std::endl;
  // Iterate through all cells
  for (int i = 0; i < num_x_; i++){ // y-index
    for (int j = 0; j < num_x_; j++){ // x-index
      // Cell Center
      std::vector<double> cell_center{geometry.XMin() + j*cell_size + cell_size/2, 
                                      geometry.YMin() + i*cell_size + cell_size/2};
      
      // Boundary Flag
      double key = boundary::helpers::MortonKey(cell_center, depth_, geometry.Maxes(), geometry.Mins());
      int covered_id = cell_map[key];

      // Boundary Length & Volume Fraction & Laplacian
      double boundary_length = 0;
      double volume_fraction = 0;
      double laplacian = 0;
      double first_x = 0;
      double first_y = 0;
      double boundary_normal_x = 0;
      double boundary_normal_y = 0;
      
      // Exterior stays 0
      // Interior
      if (covered_id == 1){
        volume_fraction = 1;
        std::vector<double> right_center = geometry.IJToCenter(i + 1, j, depth_);
        std::vector<double> left_center = geometry.IJToCenter(i - 1, j, depth_);
        std::vector<double> up_center = geometry.IJToCenter(i, j + 1, depth_);
        std::vector<double> down_center = geometry.IJToCenter(i, j - 1, depth_);
        laplacian += 1/std::pow(cell_size, 2)*(-4*phi(cell_center)
                                               + phi(right_center)
                                               + phi(left_center)
                                               + phi(up_center)
                                               + phi(down_center));
        }
      // Boundary
      else if (covered_id == 2){
        double volume_moment = geometry_info[key].volume_moments[0][0];
        // Set Volume Fraction
        volume_fraction = volume_moment/std::pow(cell_size, 2);
        // Set Boundary Length
        boundary_length = geometry_info[key].boundary_moments[0][0];
        // Calculate Laplacian
        // edge lengths
        std::array<double, 4> edge_lengths = geometry_info[key].vol_frac_1d;
        // normals
        std::vector<std::vector<std::vector<double>>> normal_derivatives = geometry_info[key].normal_derivatives;
        std::vector<double> normal = normal_derivatives[0][0];
        boundary_normal_x = normal[0];
        boundary_normal_y = normal[1];
        // first boundary moments
        first_x = geometry_info[key].boundary_moments[1][0];
        first_y = geometry_info[key].boundary_moments[0][1];
        std::vector<double> boundary_midpoint{first_x/boundary_length, first_y/boundary_length};
        // iterate through cell edges (left, up, right, down)
        for (int edge = 0; edge < 4; edge++){
          // find index of neighbor
          std::array<int, 2> neighbor_index = geometry.NeighborCell(i, j, edge);
          // find cell center of neighbor
          std::vector<double> neighbor_center = geometry.IJToCenter(neighbor_index[0], neighbor_index[1], depth_);
          // if full edge, then apply appropriate part of 5 pt stencil
          if (edge_lengths[edge] == cell_size){
            // add contributions to laplacian
            laplacian += phi(neighbor_center) - phi(cell_center);
          }
          // no edge, no flux 
          else if (edge_lengths[edge] == 0){}
          // partial edge, linearly interpolate 
          else {
            // get aperture
            double aperture = edge_lengths[edge]/cell_size;
            // interpolate with inside value
            std::array<std::array<int, 2>, 2> inter_pair = geometry.InterpolationPair(i, j, normal[0], normal[1], edge);
            // centers
            std::vector<double> first_pair_center = geometry.IJToCenter(inter_pair[0][0], inter_pair[0][1], depth_);
            std::vector<double> second_pair_center = geometry.IJToCenter(inter_pair[1][0], inter_pair[1][1], depth_);
            laplacian += aperture*((1 + aperture)/2*(phi(neighbor_center) - phi(cell_center))
                                  +(1 - aperture)/2*(phi(second_pair_center) - phi(first_pair_center)));
          }
        }
        // Boundary Flux
        // Neumann Condition should be applied at midpoint of front, not cell center
        laplacian += neumannCondition(boundary_midpoint)*boundary_length;
        
        // Scale by volume moment
        laplacian *= 1/volume_moment;
      }

      // Write to CSV
      laplace_out << covered_id << "," << cell_size << "," << cell_center[0] << "," << cell_center[1] << "," << boundary_length << "," << boundary_normal_x << "," << boundary_normal_y << "," << volume_fraction << "," << first_x << "," << first_y << "," << laplacian << std::endl;
    }
  }
  laplace_out.close();
}


Eigen::VectorXd makeLaplacian(boundary::inputs::SolverInputBase* input, 
                              boundary::geometry::Boundary boundary){
  boundary::solvers::Laplacian laplacian = boundary::solvers::Laplacian(input, boundary);
  Eigen::VectorXd solution = laplacian.solve();
  return solution;
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
  boundary::inputs::EllipseGeometry geometry_input;

  try {
    checkInput(&geometry_input);
  }
  catch (const char* msg){
    std::cerr << msg << std::endl;
  }

  // Make geometry
  boundary::geometry::Boundary boundary = boundary::geometry::Boundary(&geometry_input);
  operatorTesting(boundary);
  // std::map<int, boundary::geometry::geo_info> boundary_cells = boundary.BoundaryCells();
  // std::map<int, int> cell_map = boundary.CellMap();
  // // Make laplacian
  // boundary::inputs::LaplaceNeumann solver_input;
  // Eigen::VectorXd solution = makeLaplacian(&solver_input, boundary);
  // // Write to file
  // writeSolution(solution, boundary);
  // // writeGeometry(cell_map, boundary_cells, boundary);
  return 0;
}


