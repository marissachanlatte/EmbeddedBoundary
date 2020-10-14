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
  for (int i = 0; i < num_y_; i++){ // y-index
    for (int j = 0; j < num_x_; j++){ // x-index
      int global_id = IJToGlobal(i, j, num_x_);
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
        // edge lengths
        std::array<double, 4> edge_lengths = geometry_info[cell_center].vol_frac_1d;
        // normals
        std::vector<std::vector<std::vector<double>>> normal_derivatives = geometry_info[cell_center].normal_derivatives;
        std::vector<double> normal = normal_derivatives[0][0];
        // iterate through cell edges (left, up, right, down)
        for (int edge = 0; edge < 4; edge++){
          // find neighboring cell through this edge
          std::array<int, 2> neighbor_idx = neighborCell(i, j, edge);
          // if full edge, then apply appropriate part of 5 pt stencil
          if (edge_lengths[edge] == cell_size){
            std::cout << "edge " << edge << "full." << std::endl;
            safeMatrixAssign(global_id, IJToGlobal(neighbor_idx[0], neighbor_idx[1], num_x_), 
                             -1/std::pow(cell_size, 2) * scaling_factor);
            safeMatrixAssign(global_id, global_id, 1/std::pow(cell_size, 2) * scaling_factor);
          }
          // no edge, no flux
          else if (edge_lengths[edge] == 0){
            std::cout << "edge " << edge << "empty." << std::endl;
          }
          // partial edge, linearly interpolate 
          else {
            // get aperature
            double aperature = geometry_info[cell_center].boundary_moments[0][0];
            // this cell & neighbor
            safeMatrixAssign(global_id, IJToGlobal(neighbor_idx[0], neighbor_idx[1], num_x_), 
                             aperature*(1 + aperature)/2*scaling_factor);
            safeMatrixAssign(global_id, global_id, -aperature*(1 + aperature)/2*scaling_factor);
            // interpolate with inside value
            std::array<std::array<int, 2>, 2> inter_pair = interpolationPair(i, j, normal[0], normal[1], edge);
            safeMatrixAssign(global_id, IJToGlobal(inter_pair[0][0], inter_pair[0][1], num_x_), 
                             aperature*(1 - aperature)/2*scaling_factor);
            safeMatrixAssign(global_id, IJToGlobal(inter_pair[1][0], inter_pair[1][1], num_x_), 
                             -aperature*(1 - aperature)/2*scaling_factor);
          }
        }
        // boundary flux - use Neumann boundary conditions which gives a prescribed flux
        rhs_[global_id] += 1*scaling_factor;

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


int Laplacian::IJToGlobal(int i_index, int j_index, int num_x){
  if(i_index < 0 || j_index < 0 || i_index >= num_y_ || j_index >= num_x_){
    return -1;
  }
  return num_x*j_index + i_index;
};


std::array<int, 2> Laplacian::neighborCell(int i_index, int j_index, int edge){
  std::map<int, std::array<int, 2>> neighbor_map = {{0, {i_index - 1, j_index}},
                                                    {1, {i_index, j_index - 1}},
                                                    {2, {i_index + 1, j_index}},
                                                    {3, {i_index, j_index + 1}}};
  return neighbor_map[edge];
}


int Laplacian::sgn(double v){
    // branchless sign function.
    // https://helloacm.com/how-to-implement-the-sgn-function-in-c/
    return (v > 0) - (v < 0);
}


std::array<int, 2> Laplacian::projected_normal(int side_index, double nx, double ny){
    // (side_index & 1) iff care about y direction.
    // thus we can create a projected normal using comparison + bit mask.
    std::array<int, 2> normal = {sgn((    (side_index & 1)) * nx),
                                 sgn((1 ^ (side_index & 1)) * ny)};
    return normal;
}


int Laplacian::parity(int side_index){
    // return parity of 2 bit number
    return (side_index ^ (side_index >> 1)) & 1;
}


std::array<std::array<int, 2>, 2> Laplacian::interpolationPair(int i, int j, double nx, double ny, int side_index){
  // want to move along projected direction of negative normal.
  std::array<int, 2> steps = projected_normal(side_index, -nx, -ny);
  
  // take step along projected normal.
  int i2 = i + steps[0]; 
  int j2 = j + steps[1];

  int direction = parity(side_index) - (not parity(side_index));
  steps[0] = direction * ((side_index & 1) ^ 1);
  steps[1] = direction *  (side_index & 1);

  int i3 = i2 + steps[0];
  int j3 = j2 + steps[1];

  std::array<int, 2> first_pair = {i2, j2};
  std::array<int, 2> second_pair = {i3, j3};
  std::array<std::array<int, 2>, 2> pair = {first_pair, second_pair};
  return pair;
}



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