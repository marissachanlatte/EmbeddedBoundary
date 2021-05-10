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
}

void LaplaceOperatorHO::ComputeAndPrint(){
  int num_x_ = std::pow(2, depth_);
  // Loop through all cells at given depth
  double cell_size = geometry_.InitialCellSize()/num_x_;
  // Get cell map
  std::map<double, int> cell_map = geometry_.CellMap();
  // Get geometry information
  std::map<double, boundary::geometry::geo_info> geometry_info = geometry_.BoundaryCells();
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
      double key = boundary::helpers::MortonKey(cell_center, depth_, geometry_.Maxes(), geometry_.Mins());
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
        double volume_moment = geometry_info[key].volume_moments[0][0];
        // // Set Volume Fraction
        // volume_fraction = volume_moment/std::pow(cell_size, 2);
        // // Set Boundary Length
        // boundary_length = geometry_info[key].boundary_moments[0][0];
        // // Calculate Laplacian
        // // edge lengths
        // std::array<double, 4> edge_lengths = geometry_info[key].vol_frac_1d;
        // // normals
        // std::vector<std::vector<std::vector<double>>> normal_derivatives = geometry_info[key].normal_derivatives;
        // std::vector<double> normal = normal_derivatives[0][0];
        // boundary_normal_x = normal[0];
        // boundary_normal_y = normal[1];
        // // first boundary moments
        // first_x = geometry_info[key].boundary_moments[1][0];
        // first_y = geometry_info[key].boundary_moments[0][1];
        // std::vector<double> boundary_midpoint{first_x/boundary_length + cell_center[0], 
        //                                       first_y/boundary_length + cell_center[1]};
        // iterate through cell edges (left, up, right, down)
        for (int edge = 0; edge < 4; edge++){
          
          // Construct G
          Eigen::MatrixXf G;
          // Construct W
          Eigen::MatrixXf W;
          // Construct M
          // M is the matrix of volume moments of the neighbors normalized
          Eigen::MatrixXf M = ComputeM(cell_center);
          // Multiply together to get S (Eq. 23)
          // Eigen::VectorXf S = G*helpers::PseudoInverse(W*M)*W;
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
double LaplaceOperatorHO::NeumannCondition(std::vector<double> point, std::vector<double> normal){
  double r = std::sqrt(std::pow(point[0], 2) + std::pow(point[1], 2));
  return (normal[0]*point[0] + normal[1]*point[1])*1.0/4*std::pow(r, 2);
};

/// Compute matrix M given cell center
Eigen::MatrixXf LaplaceOperatorHO::ComputeM(std::vector<double> cell_center){
  Eigen::MatrixXf M;
  // Find Neighbors
  // Get relevant moments
  return M;
}

/// Returns a list of Morton Keys for cells in neighborhood 
std::vector<double> LaplaceOperatorHO::Neighborhood(std::vector<double> cell_center, double cell_size){
  int radius = 3;
  std::vector<double> neighbor_list;
  // iterate through all cells in radius
  for (int i = -radius; i <= radius; i ++){ // y-direction
    for (int j = -radius; j <= radius; j++){ // x-direction
      std::vector<double> neighbor_center{cell_center[0] + j*cell_size, cell_center[1] + i*cell_size};
      double neighbor_key = helpers::MortonKey(neighbor_center, depth_, geometry_.Maxes(), geometry_.Mins());
      // Check if neighbor cell is inside the domain and add
      if (cell_map_[neighbor_key] > 0){
        neighbor_list.push_back(neighbor_key);
      }
    }
  }
  return neighbor_list;
}
} // namespace solvers

} // namespace boundary
