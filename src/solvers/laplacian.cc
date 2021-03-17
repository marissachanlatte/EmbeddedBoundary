#include "solvers/laplacian.h"
#include "helpers/math_helpers.h"

#define _USE_MATH_DEFINES
 
#include <iostream>
#include <cmath>

namespace boundary {

namespace solvers {

Laplacian::Laplacian(boundary::inputs::SolverInputBase* input, boundary::geometry::Boundary geometry){
  depth_ = geometry.MaxSolverDepth();
  num_x_ = std::pow(2, depth_);
  matrix_.resize(std::pow(num_x_, 2), std::pow(num_x_, 2));
  rhs_ = Eigen::VectorXd::Zero(std::pow(num_x_, 2));
  BuildMatrix(input, geometry);
};


void Laplacian::BuildMatrix(boundary::inputs::SolverInputBase* input, boundary::geometry::Boundary geometry){
  double cell_size = geometry.InitialCellSize()/std::pow(2, depth_);
  // Get cell map
  std::map<int, int> cell_map = geometry.CellMap();
  // Get geometry information
  std::map<int, geometry::geo_info> geometry_info = geometry.BoundaryCells();
  // Iterate through all cells
  for (int i = 0; i < num_x_; i++){ // y-index
    for (int j = 0; j < num_x_; j++){ // x-index
      int matrix_id = geometry.IJToGlobal(i, j, depth_);
      // get cell center
      std::vector<double> cell_center{geometry.XMin() + j*cell_size + cell_size/2, 
                                      geometry.YMin() + i*cell_size + cell_size/2};
      int key = helpers::MortonKey(cell_center, depth_, geometry.Maxes(), geometry.Mins());
      int covered_id = cell_map[key];
      // If boundary cell
      if (covered_id == 2){
        // get volume moment (note, actual volume, not fraction)
        double volume_moment = geometry_info[key].volume_moments[0][0];
        // scaling factor
        double scaling_factor = 1/volume_moment;
        // edge lengths
        std::array<double, 4> edge_lengths = geometry_info[key].vol_frac_1d;
        // normals
        std::vector<std::vector<std::vector<double>>> normal_derivatives = geometry_info[key].normal_derivatives;
        std::vector<double> normal = normal_derivatives[0][0];
        // iterate through cell edges (left, up, right, down)
        for (int edge = 0; edge < 4; edge++){
          
          // find neighboring cell through this edge
          std::array<int, 2> neighbor_idx = geometry.NeighborCell(i, j, edge);
          // if full edge, then apply appropriate part of 5 pt stencil
          if (edge_lengths[edge] == cell_size){
            SafeMatrixAssign(matrix_id, geometry.IJToGlobal(neighbor_idx[0], neighbor_idx[1], depth_), 
                             scaling_factor);
            SafeMatrixAssign(matrix_id, matrix_id, -scaling_factor);
          }
          // no edge, no flux
          else if (edge_lengths[edge] == 0){
            SafeMatrixAssign(matrix_id, geometry.IJToGlobal(neighbor_idx[0], neighbor_idx[1], depth_), 
                             scaling_factor);
            SafeMatrixAssign(matrix_id, matrix_id, -scaling_factor);
          }
          // partial edge, linearly interpolate 
          else {
            // get aperature
            double aperature = geometry_info[key].boundary_moments[0][0];
            // this cell & neighbor
            SafeMatrixAssign(matrix_id, geometry.IJToGlobal(neighbor_idx[0], neighbor_idx[1], depth_), 
                             aperature*(1 + aperature)/2*scaling_factor);
            SafeMatrixAssign(matrix_id, matrix_id, -aperature*(1 + aperature)/2*scaling_factor);
            // interpolate with inside value
            std::array<std::array<int, 2>, 2> inter_pair = geometry.InterpolationPair(i, j, normal[0], normal[1], edge);
            SafeMatrixAssign(matrix_id, geometry.IJToGlobal(inter_pair[0][0], inter_pair[0][1], depth_), 
                             -aperature*(1 - aperature)/2*scaling_factor);
            SafeMatrixAssign(matrix_id, geometry.IJToGlobal(inter_pair[1][0], inter_pair[1][1], depth_), 
                             aperature*(1 - aperature)/2*scaling_factor);
          }

        }
        // Boundary Flux - Neumann Condition
        double neumann_condition = input->NeumannCondition(cell_center[0], cell_center[1]);
        rhs_[matrix_id] += -scaling_factor*neumann_condition;

        // Assign Right Hand Side
        rhs_[matrix_id] += input->RightHandSide(cell_center[0], cell_center[1]);
      }

      // If interior, five point stencil
      else if (covered_id == 1){
        SafeMatrixAssign(matrix_id, geometry.IJToGlobal(i + 1, j, depth_), -1/std::pow(cell_size, 2));
        SafeMatrixAssign(matrix_id, geometry.IJToGlobal(i - 1, j, depth_), -1/std::pow(cell_size, 2));
        SafeMatrixAssign(matrix_id, geometry.IJToGlobal(i, j + 1, depth_), -1/std::pow(cell_size, 2));
        SafeMatrixAssign(matrix_id, geometry.IJToGlobal(i, j - 1, depth_), -1/std::pow(cell_size, 2));
        SafeMatrixAssign(matrix_id, matrix_id, 4/std::pow(cell_size, 2));
        
        // Assign Right Hand Side
        rhs_[matrix_id] += input->RightHandSide(cell_center[0], cell_center[1]);
      }
      // If exterior set to 0
      else if (covered_id == 0){
        // SafeMatrixAssign(matrix_id, matrix_id, 1);
      }
    }
  }
};


void Laplacian::SafeMatrixAssign(int i_index, int j_index, double value){
  if ((i_index >= 0) && (i_index < num_x_*num_x_) && (j_index >= 0) && (j_index < num_x_*num_x_)){
    matrix_.coeffRef(i_index, j_index) += value;
  }
};


Eigen::VectorXd Laplacian::solve(){
  // Solve system using SuperLU
  Eigen::VectorXd solution(num_x_*num_x_);
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(matrix_); 
  solution = solver.solve(rhs_); 
  return solution;
};

} // namespace solvers

} // namespace boundary