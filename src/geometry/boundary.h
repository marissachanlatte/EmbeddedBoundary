#ifndef EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H
#define EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H

#include "inputs/input_base.h"

#include <array>
#include <map>

namespace boundary {

namespace geometry {

  struct geo_info {
    bool irregular; // tells if cell is irregular or not
    int id; // id to index into vector with more cell information
  };

  class Boundary{
    public:
      Boundary(boundary::inputs::InputBase* input);
      static bool IsBoundaryCell(std::array<double, 2> lower_left,
                                 std::array<double, 2> lower_right,
                                 std::array<double, 2> upper_right,
                                 std::array<double, 2> upper_left,
                                 boundary::inputs::InputBase* input);
    private:
      std::map<std::array<double, 2>, geo_info> boundary_cells_;
  };
}

}
#endif // EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H
