#include "geometry/boundary.h"

#include <array>

namespace boundary {

namespace geometry {

Boundary::Boundary(boundary::inputs::InputBase* input){
  double cell_size = input->CellSize();
  double y_min = input->YMin();
  double y_max = y_min + cell_size;
  int id_count = 0;
  while (y_max < input->YMax()){
    double x_min = input->XMin();
    double x_max = x_min + cell_size;
    while (x_max < input->XMax()){
      // corners of cell
      std::array<double, 2> lower_left{x_min, y_min};
      std::array<double, 2> upper_left{x_min, y_max};
      std::array<double, 2> upper_right{x_max, y_max};
      std::array<double, 2> lower_right{x_max, y_min};

      bool is_boundary = IsBoundaryCell(lower_left, lower_right,
                                        upper_right, upper_left, input);
      if (is_boundary){
        // make struct with tag and id
        geo_info cell;
        cell.irregular = true;
        cell.id = id_count;
        // add point to map
        std::array<double, 2> point = {x_min + cell_size/2, y_min + cell_size/2};
        boundary_cells_.insert(
          std::pair<std::array<double, 2>, geo_info>(point, cell));
      }
    }
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

} // namespace geometry

} // namespace boundary
