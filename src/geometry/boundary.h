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

/// A class describing the boundary geometry.
/**
This class stores a map of all boundary cells with necessary geometry information
*/
  class Boundary{
    public:
      /// Constructor
      /**
      Iterates through all cells in geometry and adds boundary cells to map
      */
      Boundary(boundary::inputs::InputBase* input);
      /// Determines if a cell is a boundary cell
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
