#include "solvers/laplace_operator.h"
#include "helpers/math_helpers.h"

#include <iostream>
#include <fstream>

namespace boundary {

namespace solvers {

/// Outputs laplacian, volume fractions, boundary length, cell centers, boundary flag, and cell size as CSV file
LaplaceOperator::LaplaceOperator(boundary::geometry::Boundary geometry){
  int depth_ = geometry.MaxSolverDepth();
  int num_x_ = std::pow(2, depth_);
  // Loop through all cells at given depth
  double cell_size = geometry.InitialCellSize()/num_x_;
  // Get cell map
  std::map<double, int> cell_map = geometry.CellMap();
  // Get geometry information
  std::map<double, geometry::geo_info> geometry_info = geometry.BoundaryCells();
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
      double key = helpers::MortonKey(cell_center, depth_, geometry.Maxes(), geometry.Mins());
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
        laplacian = 1/std::pow(cell_size, 2)*(-4*Phi(cell_center)
                                               + Phi(right_center)
                                               + Phi(left_center)
                                               + Phi(up_center)
                                               + Phi(down_center));
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
        std::vector<double> boundary_midpoint{first_x/boundary_length + cell_center[0], 
                                              first_y/boundary_length + cell_center[1]};
        // iterate through cell edges (left, up, right, down)
        for (int edge = 0; edge < 4; edge++){
          // find index of neighbor
          std::array<int, 2> neighbor_index = geometry.NeighborCell(i, j, edge);
          // find cell center of neighbor
          std::vector<double> neighbor_center = geometry.IJToCenter(neighbor_index[0], neighbor_index[1], depth_);
          // if full edge, then apply appropriate part of 5 pt stencil
          if (edge_lengths[edge] == cell_size){
            // add contributions to laplacian
            laplacian += Phi(neighbor_center) - Phi(cell_center);
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
            laplacian += aperture*((1 + aperture)/2*(Phi(neighbor_center) - Phi(cell_center))
                                  +(1 - aperture)/2*(Phi(second_pair_center) - Phi(first_pair_center)));
          }
        }
        // Boundary Flux
        // Neumann Condition should be applied at midpoint of front, not cell center
        laplacian += NeumannCondition(boundary_midpoint, normal)*boundary_length;
        
        // Scale by volume moment
        laplacian *= 1/volume_moment;
      }

      // Write to CSV
      laplace_out << covered_id << "," << cell_size << "," << cell_center[0] << "," << cell_center[1] << "," << boundary_length << "," << boundary_normal_x << "," << boundary_normal_y << "," << volume_fraction << "," << first_x << "," << first_y << "," << laplacian << std::endl;
    }
  }
  laplace_out.close();
}

/// Function Phi for Testing
double LaplaceOperator::Phi(std::vector<double> point){
  return 1.0/16*std::pow((std::pow(point[0], 2) + std::pow(point[1], 2)), 2);
};


/// BC for Phi for Testing
double LaplaceOperator::NeumannCondition(std::vector<double> point, std::vector<double> normal){
  double r = std::sqrt(std::pow(point[0], 2) + std::pow(point[1], 2));
  return (normal[0]*point[0] + normal[1]*point[1])*1.0/4*std::pow(r, 2);
};

} // namespace solvers

} // namespace boundary
