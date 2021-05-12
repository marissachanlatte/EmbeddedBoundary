#include "solvers/laplace_operator_ho.h"
#include "helpers/math_helpers.h"

#include <iostream>
#include <fstream>

namespace boundary {

namespace solvers {

/// Outputs laplacian, volume fractions, boundary length, cell centers, boundary flag, and cell size as CSV file
LaplaceOperatorHO::LaplaceOperatorHO(boundary::geometry::Boundary geometry){
  depth_ = geometry.MaxSolverDepth();
  geometry_ = geometry;
  cell_map_ = geometry.CellMap();
  geometry_info_ = geometry.BoundaryCells();
  Q_ = geometry_.Q();
  num_q_ = (Q_ + 1)*(Q_ + 2)/2;
}

void LaplaceOperatorHO::ComputeAndPrint(){
  int num_x_ = std::pow(2, depth_);
  // Loop through all cells at given depth
  double cell_size = geometry_.InitialCellSize()/num_x_;
  // Get cell map
  std::map<double, int> cell_map = geometry_.CellMap();
  // Output file
  std::ofstream laplace_out;
  laplace_out.open ("../outputs/laplace_out.txt");
  // Headers
  laplace_out << "Covered ID,Cell Size,CenterX,CenterY,Boundary Length,Boundary NormalX, Boundary NormalY,Volume Fraction,FirstX,FirstY,Laplacian" << std::endl;
  // Iterate through all cells
  for (int i = 0; i < num_x_; i++){ // y-index
    for (int j = 0; j < num_x_; j++){ // x-index
      // Cell Center
      std::vector<double> cell_center{geometry_.XMin() + j*cell_size + cell_size/2, 
                                      geometry_.YMin() + i*cell_size + cell_size/2};
      
      // Boundary Flag
      double key = boundary::helpers::MortonKey(cell_center, 
                                                depth_, 
                                                geometry_.Maxes(), 
                                                geometry_.Mins());
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
      // Interior & Boundary
      if (covered_id > 0){
        double volume_moment = geometry_info_[key].volume_moments[0][0];
        // iterate through cell edges (left, up, right, down)
        // Find Neighbors
        std::vector<double> neighbors = Neighborhood(cell_center, cell_size);
        for (int edge = 0; edge < 4; edge++){         
          // Construct G
          Eigen::RowVectorXf G = ComputeG(edge, key);
          // Construct W
          Eigen::MatrixXf W = ComputeW(cell_center, cell_size, neighbors);
          // Construct M
          // M is the matrix of volume moments of the neighbors normalized
          Eigen::MatrixXf M = ComputeM(cell_center, cell_size, neighbors);
          // Multiply together to get S (Eq. 23)
          Eigen::VectorXf S = G*helpers::PseudoInverse(W*M)*W;
          // TODO: Why is S a vector? how do we convert to a double for laplacian?
          // laplacian += S;
        }
        // TODO: If near boundary, add boundary flux
        // laplacian += NeumannCondition(boundary_midpoint, normal)*boundary_length;
        
        // Scale by cell size
        laplacian *= 1/std::pow(cell_size, 2);

        // Scale by volume moment
        laplacian *= 1/volume_moment;
      }

      // Write to CSV
      laplace_out << covered_id << "," << cell_size << "," << cell_center[0] << "," << cell_center[1] << "," << boundary_length << "," << boundary_normal_x << "," << boundary_normal_y << "," << volume_fraction << "," << first_x << "," << first_y << "," << laplacian << std::endl;
    }
  }
  laplace_out.close();
};


/// Function Phi for Testing
double LaplaceOperatorHO::Phi(std::vector<double> point){
  return 1.0/16*std::pow((std::pow(point[0], 2) + std::pow(point[1], 2)), 2);
};


/// BC for Phi for Testing
double LaplaceOperatorHO::NeumannCondition(std::vector<double> point, 
                                           std::vector<double> normal){
  double r = std::sqrt(std::pow(point[0], 2) + std::pow(point[1], 2));
  return (normal[0]*point[0] + normal[1]*point[1])*1.0/4*std::pow(r, 2);
};

/// Compute matrix M given cell center
Eigen::MatrixXf LaplaceOperatorHO::ComputeM(std::vector<double> cell_center, 
                                            double cell_size,
                                            std::vector<double> neighbors){
  int num_neighbors = neighbors.size();
  // Initialize M
  Eigen::MatrixXf M(num_neighbors, num_q_);
  // Get relevant moments
  for (int i = 0; i < num_neighbors; i++){
    // If boundary cell, get moments for cell i
    std::vector<std::vector<double>>  moments;
    if (cell_map_[neighbors[i]] == 2){ 
      moments = geometry_info_[neighbors[i]].volume_moments;
    }
    int j = 0;
    M(i, j) = 1;
    // Loop through in lexicographical order
    for (int a = 0; a < (Q_ + 1); a++){
      for (int b = 0; b < (Q_ + 1); b++){
        if ((a + b) == 0 || (a + b) > Q_){
          continue;
        }
        j += 1;
        if (cell_map_[neighbors[i]] == 2){
           M(i, j) = moments[b][a]/moments[0][0];
        }
        else {
          std::vector<int> p{a, b};
          M(i, j) = helpers::PMoment(p, cell_size)/std::pow(cell_size, 2);
        }
      }
    }
  }
  return M;
}


Eigen::MatrixXf LaplaceOperatorHO::ComputeW(std::vector<double> cell_center, 
                                            double cell_size,
                                            std::vector<double> neighbors){
  int num_neighbors = neighbors.size();
  Eigen::MatrixXf W(num_neighbors, num_neighbors);
  // Fill in weights with distance between neighbors and cell center
  for (int i = 0; i < num_neighbors; i++){
    std::vector<double> neighbor_center = geometry_info_[neighbors[i]].cell_center;
    int d = std::sqrt(std::pow(cell_center[0] - neighbor_center[1], 2) +
                      std::pow(cell_center[1] - neighbor_center[1], 2));
    W(i, i) = std::pow(d, -5);
  }
  return W;
}


/// Returns a list of Morton Keys for cells in neighborhood 
std::vector<double> LaplaceOperatorHO::Neighborhood(std::vector<double> cell_center, 
                                                    double cell_size){
  int radius = 3;
  std::vector<double> neighbor_list;
  // iterate through all cells in radius
  for (int i = -radius; i <= radius; i ++){ // y-direction
    for (int j = -radius; j <= radius; j++){ // x-direction
      std::vector<double> neighbor_center{cell_center[0] + j*cell_size, 
                                          cell_center[1] + i*cell_size};
      double neighbor_key = helpers::MortonKey(neighbor_center, 
                                               depth_, 
                                               geometry_.Maxes(), 
                                               geometry_.Mins());
      // Check if neighbor cell is inside the domain and add
      if (cell_map_[neighbor_key] > 0){
        neighbor_list.push_back(neighbor_key);
      }
    }
  }
  return neighbor_list;
}


Eigen::RowVectorXf LaplaceOperatorHO::ComputeG(int edge, double cell_id){
  std::vector<std::vector<double>> moments = geometry_info_[cell_id].volume_moments;
  std::array<int, 2> e{0, 0};
  Eigen::RowVectorXf G(num_q_);
  int d = ((edge == 0 || edge == 2) ? 1 : 0);
  e[d] = 1;
  int i = 0; // place in row
  for (int a = 0; a < (Q_ + 1); a++){
    for (int b = 0; b < (Q_ + 1); b++){
      if ((a + b) > Q_){
        continue;
      }
      std::vector<int> q{a, b};
      std::vector<int> p{a - e[0], b - e[1]};
      G(i) = q[d]*moments[p[0]][p[1]];
      i++;
    }
  }
}
                                          

} // namespace solvers

} // namespace boundary
