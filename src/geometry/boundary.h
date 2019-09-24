#ifndef EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H
#define EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H

#include "inputs/input_base.h"

#include <array>
#include <map>
#include <vector>

namespace boundary {

namespace geometry {

  struct geo_info {
    /// tells if cell is irregular or not
    bool irregular;
    /// id to index into vector with more cell information
    int id;
    /// derivatives of the normal to the boundary
    std::vector<std::vector<std::vector<double>>> normal_derivatives;
    /// 1d volume fraction for boundary cell edges
    std::array<double, 4> vol_frac_1d;
    /// volume moments
    std::vector<std::vector<double>> volume_moments;
    /// boundary moments
    std::vector<std::vector<double>> boundary_moments;
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
      std::map<std::array<double, 2>, geo_info> BoundaryCells();

    private:
      void CalculateMoments_(std::array<double, 2> cell_center);
      double DIntegral_(double beginning,
                        double end,
                        std::array<int, 2> q,
                        int index,
                        double fixed_value);
      double CalcD_(double bd_length,
                    double fixed_value,
                    std::array<double, 2> cell_center,
                    std::array<int, 2> q,
                    int d,
                    std::array<int, 2> which_d);
      std::map<std::array<double, 2>, geo_info> boundary_cells_;
      boundary::inputs::InputBase* input_;
      int Q_;
      double cell_size_;
  };
}

}
#endif // EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H
