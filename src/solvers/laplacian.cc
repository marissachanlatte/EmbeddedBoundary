#include "solvers/laplacian.h"

#include <iostream>

namespace boundary {

namespace solvers {

Laplacian::Laplacian(boundary::geometry::Boundary geometry){
  num_x_ = geometry.NumX();
  num_y_ = geometry.NumY();
  matrix_.resize(num_x_*num_y_, num_x_*num_y_);
  rhs_ = Eigen::VectorXd::Zero(num_x_*num_y_);
  BuildMatrix(geometry);
};


void Laplacian::BuildMatrix(boundary::geometry::Boundary geometry){
  double cell_size = geometry.CellSize();
  // Get cell map
  std::map<int, int> cell_map = geometry.CellMap();
  // Get geometry information
  std::map<std::array<double, 2>, geometry::geo_info> geometry_info = geometry.BoundaryCells();
  // Iterate through all cells
  for (int i = 0; i < num_y_; i++){ // y-index
    for (int j = 0; j < num_x_; j++){ // x-index
      int global_id = geometry.IJToGlobal(i, j);
      int covered_id = cell_map[global_id];
      // If boundary cell
      if (covered_id == 2){
        std::cout << "boundary cell " << i << " " << j << std::endl;
        // get cell center
        std::array<double, 2> cell_center = {geometry.XMin() + i*cell_size + cell_size/2, 
                                              geometry.YMin() + j*cell_size + cell_size/2};
        // get volume moment (note, actual volume, not fraction)
        double volume_moment = geometry_info[cell_center].volume_moments[0][0];
        // scaling factor
        double scaling_factor = 1/volume_moment;
        std::cout << "scaling factor" << scaling_factor << std::endl;
        // edge lengths
        std::array<double, 4> edge_lengths = geometry_info[cell_center].vol_frac_1d;
        // normals
        std::vector<std::vector<std::vector<double>>> normal_derivatives = geometry_info[cell_center].normal_derivatives;
        std::vector<double> normal = normal_derivatives[0][0];
        // iterate through cell edges (left, up, right, down)
        for (int edge = 0; edge < 4; edge++){
          // find neighboring cell through this edge
          std::array<int, 2> neighbor_idx = geometry.NeighborCell(i, j, edge);
          // if full edge, then apply appropriate part of 5 pt stencil
          if (edge_lengths[edge] == cell_size){
            SafeMatrixAssign(global_id, geometry.IJToGlobal(neighbor_idx[0], neighbor_idx[1]), 
                             -1/std::pow(cell_size, 2) * scaling_factor);
            SafeMatrixAssign(global_id, global_id, 1/std::pow(cell_size, 2) * scaling_factor);
          }
          // no edge, no flux
          else if (edge_lengths[edge] == 0){
          }
          // partial edge, linearly interpolate 
          else {
            // get aperature
            double aperature = geometry_info[cell_center].boundary_moments[0][0];
            // this cell & neighbor
            std::cout << "cell 1 " << geometry.IJToGlobal(neighbor_idx[0], neighbor_idx[1]) << std::endl;
            SafeMatrixAssign(global_id, geometry.IJToGlobal(neighbor_idx[0], neighbor_idx[1]), 
                             aperature*(1 + aperature)/2*scaling_factor);
            SafeMatrixAssign(global_id, global_id, -aperature*(1 + aperature)/2*scaling_factor);
            // interpolate with inside value
            std::array<std::array<int, 2>, 2> inter_pair = geometry.InterpolationPair(i, j, normal[0], normal[1], edge);
            SafeMatrixAssign(global_id, geometry.IJToGlobal(inter_pair[0][0], inter_pair[0][1]), 
                             aperature*(1 - aperature)/2*scaling_factor);
            SafeMatrixAssign(global_id, geometry.IJToGlobal(inter_pair[1][0], inter_pair[1][1]), 
                             -aperature*(1 - aperature)/2*scaling_factor);
            std::cout << "cell 2 " << geometry.IJToGlobal(inter_pair[0][0], inter_pair[0][1]) << std::endl;
            std::cout << "cell 3 " << geometry.IJToGlobal(inter_pair[1][0], inter_pair[1][1]) << std::endl;
          }
        }
        // boundary flux - use Neumann boundary conditions which gives a prescribed flux
        rhs_[global_id] += 1*scaling_factor;

        // Set RHS
        rhs_[global_id] = 1;
      }
      // If interior, five point stencil
      else if (covered_id == 1){
        SafeMatrixAssign(global_id, geometry.IJToGlobal(i + 1, j), -1/std::pow(cell_size, 2));
        SafeMatrixAssign(global_id, geometry.IJToGlobal(i - 1, j), -1/std::pow(cell_size, 2));
        SafeMatrixAssign(global_id, geometry.IJToGlobal(i, j + 1), -1/std::pow(cell_size, 2));
        SafeMatrixAssign(global_id, geometry.IJToGlobal(i, j - 1), -1/std::pow(cell_size, 2));
        SafeMatrixAssign(global_id, global_id, 4/std::pow(cell_size, 2));

        // Set RHS
        rhs_[global_id] = 1;
      }
      // If exterior set to 0
      else if (covered_id == 0){
        SafeMatrixAssign(global_id, global_id, 1);
      }
    }
  }
};


void Laplacian::SafeMatrixAssign(int i_index, int j_index, double value){
  if ((i_index >= 0) && (i_index < num_x_*num_y_) && (j_index >= 0) && (j_index < num_x_*num_y_)){
    matrix_.coeffRef(i_index, j_index) += value;
  }
};


Eigen::VectorXd Laplacian::solve(){
  // Solve system using SuperLU
  Eigen::VectorXd solution(num_x_*num_y_);
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(matrix_); 
  std::cout << matrix_ << std::endl;
  solution = solver.solve(rhs_); 
  return solution;
};

} // namespace solvers

} // namespace boundary