#include "solvers/laplacian.h"

#include <iostream>

namespace boundary {

namespace solvers {

Laplacian::Laplacian(boundary::geometry::Boundary geometry){
  num_x_ = int(std::abs(geometry.XMax() - geometry.XMin())/geometry.CellSize());
  num_y_ = int(std::abs(geometry.YMax() - geometry.YMin())/geometry.CellSize());
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
  for (int j = 0; j < num_x_; j++){ // y-index
    for (int i = 0; i < num_y_; i++){ // x-index
      int global_id = IJToGlobal(i, j, num_x_);
      int covered_id = cell_map[global_id];
      // If boundary cell
      if (covered_id == 2){
        // get cell center
        std::array<double, 2> cell_center = {geometry.XMin() + i*cell_size + cell_size/2, 
                                              geometry.YMin() + j*cell_size + cell_size/2};
        // get volume fraction
        double volume_fraction = geometry_info[cell_center].volume_moments[0][0];
        // get aperature
        double aperature = geometry_info[cell_center].boundary_moments[0][0];
        // scaling factor
        double scaling_factor = 1/(std::pow(cell_size, 2)*volume_fraction);
        // right flux
        safeMatrixAssign(global_id, IJToGlobal(i + 1, j, num_x_), aperature*(1 + aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, global_id, -aperature*(1 + aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, IJToGlobal(i + 1, j + 1, num_x_), aperature*(1 - aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, IJToGlobal(i, j + 1, num_x_), -aperature*(1 - aperature)/2*scaling_factor);
        // left flux
        safeMatrixAssign(global_id, IJToGlobal(i - 1, j, num_x_), -aperature*(1 + aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, global_id, aperature*(1 + aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, IJToGlobal(i - 1, j + 1, num_x_), -aperature*(1 - aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, IJToGlobal(i, j + 1, num_x_), aperature*(1 - aperature)/2*scaling_factor);
        // top flux
        safeMatrixAssign(global_id, IJToGlobal(i, j + 1, num_x_), aperature*(1 + aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, global_id, -aperature*(1 + aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, IJToGlobal(i + 1, j + 1, num_x_), aperature*(1 - aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, IJToGlobal(i + 1, j, num_x_), -aperature*(1 - aperature)/2*scaling_factor);
        // bottom flux
        safeMatrixAssign(global_id, IJToGlobal(i, j - 1, num_x_), -aperature*(1 + aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, global_id, aperature*(1 + aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, IJToGlobal(i + 1, j - 1, num_x_), -aperature*(1 - aperature)/2*scaling_factor);
        safeMatrixAssign(global_id, IJToGlobal(i + 1, j, num_x_), aperature*(1 - aperature)/2*scaling_factor);
        // boundary flux - use Neumann boundary conditions which gives a prescribed flux
        safeMatrixAssign(global_id, global_id, 1);

        // Set RHS
        rhs_[global_id] = 1;
      }
      // If interior, five point stencil
      else if (covered_id == 1){
        safeMatrixAssign(global_id, IJToGlobal(i + 1, j, num_x_), -1/std::pow(cell_size, 2));
        safeMatrixAssign(global_id, IJToGlobal(i - 1, j, num_x_), -1/std::pow(cell_size, 2));
        safeMatrixAssign(global_id, IJToGlobal(i, j + 1, num_x_), -1/std::pow(cell_size, 2));
        safeMatrixAssign(global_id, IJToGlobal(i, j - 1, num_x_), -1/std::pow(cell_size, 2));
        safeMatrixAssign(global_id, global_id, 4/std::pow(cell_size, 2));

        // Set RHS
        rhs_[global_id] = 1;
      }
      // If exterior set to 0
      else if (covered_id == 0){
        safeMatrixAssign(global_id, global_id, 1);
      }
    }
  }
};


void Laplacian::safeMatrixAssign(int i_index, int j_index, double value){
  if ((i_index >= 0) && (i_index < num_x_*num_y_) && (j_index >= 0) && (j_index < num_x_*num_y_)){
    matrix_.coeffRef(i_index, j_index) += value;
  }
};


int Laplacian::IJToGlobal(int x_index, int y_index, int num_x){
  if(x_index < 0 || y_index < 0 || x_index >= num_x_ || y_index >= num_y_){
    return -1;
  }
  return num_x*y_index + x_index;
};


Eigen::VectorXd Laplacian::solve(){
  // Solve system using SuperLU
  Eigen::VectorXd solution(num_x_*num_y_);
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(matrix_); 
  solution = solver.solve(rhs_); 
  return solution;
};

} // namespace solvers

} // namespace boundary