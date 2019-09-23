#include "geometry/boundary.h"
#include "helpers/geometry_objects.h"
#include "normals/normals.h"

#include <array>
#include <iostream>

namespace boundary {

namespace geometry {

Boundary::Boundary(boundary::inputs::InputBase* input){
  double cell_size = input->CellSize();
  double y_min = input->YMin();
  double y_max = y_min + cell_size;
  int id_count = 0;
  while (y_max <= input->YMax()){
    double x_min = input->XMin();
    double x_max = x_min + cell_size;
    while (x_max <= input->XMax()){
      // corners of cell
      std::array<std::array<double, 2>, 4> corners;
      corners[0] = {x_min, y_min}; // lower left
      corners[1] = {x_min, y_max}; // upper left
      corners[2] = {x_max, y_max}; // upper right
      corners[3] = {x_max, y_min}; // lower right

      bool is_boundary = IsBoundaryCell(corners[0], corners[1],
                                        corners[2], corners[3], input);
      if (is_boundary){
        // make struct with tag and id
        geo_info cell;
        cell.irregular = true;
        cell.id = id_count;
        // update id count
        id_count += 1;
        // resize vector of normal normal_derivatives
        helpers::Point cell_center = helpers::Point(x_max - cell_size, y_max - cell_size);
        int Q = input->QOrder();
        cell.normal_derivatives.resize(Q+1);
        for (int i=0; i < Q+1; ++i){
          cell.normal_derivatives[i].resize(Q+1);

          for (int j = 0; j < Q + 1; ++j){
            cell.normal_derivatives[i][j].resize(2);
          }
        }
        // Compute derivatives of the normal and store in cell struct
        for (int q_mag=Q; q_mag >= 0 ; q_mag--){
          for (int q1 = 0; q1 < q_mag + 1; q1++){
            std::vector<int> q = {q1, q_mag - q1};
            cell.normal_derivatives[q1][q_mag - q1][0] = normals::Normal::NormalDerivative(q,
                                                            1,
                                                            cell_center,
                                                            input);
            cell.normal_derivatives[q1][q_mag - q1][1] = normals::Normal::NormalDerivative(q,
                                                            2,
                                                            cell_center,
                                                            input);
          }
        }
        // // Compute 1d volume fractions and store in cell
        // // Check if four corners are inside or outside boundary
        // std::vector<int> inside{input->Inside(corners[0]),
        //                         input->Inside(corners[1]),
        //                         input->Inside(corners[2]),
        //                         input->Inside(corners[3])};
        // for (int edge=0; edge < 4; edge++){
        //   // iterate four edges to determine which ones intersect boundary
        //   if (inside[edge] == inside[(edge+1)%4]){ // no change from corner to corner, TODO: come up with more robust way to do this
        //     cell.vol_frac_1d[edge] = cell_size*inside[edge];
        //   }
        //   else { // change from corner to corner, indicating boundary cuts through
        //     // check if edge is horizontal or vertical
        //     if ((corners[edge][0] - corners[(edge+1)%4][0]) == 0){ // horizontal
        //       // find the intersection of x=corners[i, 0]
        //       double y = input->BoundaryFunction(corners[edge][0]);
        //       double edge_inside = abs(y - corners[edge][1]);
        //       if (inside[edge]){
        //         cell.vol_frac_1d[edge] = edge_inside;
        //       }
        //       else{
        //         cell.vol_frac_1d[edge] = cell_size - edge_inside;
        //       }
        //     }
        //     else { // vertical
        //       // find intersection of y=corners[i, 1]
        //       double x = input->BoundaryInverse(corners[edge][1]);
        //       double edge_inside = abs(x - corners[edge][0]);
        //       if (inside[edge]){
        //         cell.vol_frac_1d[edge] = edge_inside;
        //       }
        //       else {
        //         cell.vol_frac_1d[edge] = cell_size - edge_inside;
        //       }
        //     }
        //   }
        // }
        // add point to map
        std::array<double, 2> point = {x_min + cell_size/2, y_min + cell_size/2};
        boundary_cells_.insert(
          std::pair<std::array<double, 2>, geo_info>(point, cell));
      }
      // Go to next cell
      x_max += cell_size;
      x_min += cell_size;
    }
    // Go to next row
    y_max += cell_size;
    y_min += cell_size;
  }
};

bool Boundary::IsBoundaryCell(std::array<double, 2> lower_left,
                              std::array<double, 2> lower_right,
                              std::array<double, 2> upper_right,
                              std::array<double, 2> upper_left,
                              boundary::inputs::InputBase* input){
  // put into array with inside values
  std::vector<int> inside{input->Inside(lower_left),
                          input->Inside(upper_left),
                          input->Inside(upper_right),
                          input->Inside(lower_right)};
  // try changing all 2s to 0s and see if they are all the same
  int it_prev = inside[0];
  bool same;
  for (int it = 1; it < 4; it++){
    int current = inside[it];
    if (current == 2){
      current = 0;
    }
    if (it_prev == current){
      if (it == 3){
        // they're all the same, not a boundary cell
        same = true;
      }
      continue;
    }
    else{
      same = false;
      break;
    }
  }
  // if they're not try changing all the 2s to 1s and see if they're all the same
  if (!same){
    it_prev = inside[0];
    for (int it = 1; it < 4; it++){
      int current = inside[it];
      if (current == 2){
        current = 1;
      }
      if (it_prev == current){
        if (it == 3){
          same = true;
        }
        continue;
      }
      else{
        same = false;
        break;
      }
    }
  }
  return !same;
}


std::map<std::array<double, 2>, geo_info> Boundary::BoundaryCells(){
  return boundary_cells_;
}


} // namespace geometry

} // namespace boundary
