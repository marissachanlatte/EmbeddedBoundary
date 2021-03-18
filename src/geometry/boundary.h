#ifndef EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H
#define EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H

#include "inputs/geometries/geometry_input_base.h"

#include <array>
#include <map>
#include <vector>

namespace boundary {

namespace geometry {
/**
This struct contains all geometry information for a given cell
**/

  struct geo_info {
    /// tells if cell is irregular or not
    bool irregular;
    /// derivatives of the normal to the boundary
    std::vector<std::vector<std::vector<double>>> normal_derivatives;
    /// 1d volume fraction for boundary cell edges (left, up, right, down)
    std::array<double, 4> vol_frac_1d;
    /// volume moments
    std::vector<std::vector<double>> volume_moments;
    /// boundary moments
    std::vector<std::vector<double>> boundary_moments;
    /// cell center
    std::vector<double> cell_center;
    /// cell size
    double cell_size;
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
      Boundary(boundary::inputs::GeometryInputBase* input);
      /// Determines if a cell is a boundary cell
      static bool IsBoundaryCell(std::array<double, 2> lower_left,
                                 std::array<double, 2> lower_right,
                                 std::array<double, 2> upper_right,
                                 std::array<double, 2> upper_left,
                                 boundary::inputs::GeometryInputBase* input);
      std::map<int, geo_info> BoundaryCells();
      /// Tells whether a cell is 0 - exterior, 1 - interior, or 2 - boundary
      std::map<int, int> CellMap();
      static double WhichValue(std::vector<double> values,
                         double first_bound,
                         double second_bound);
      /// Given a cell ID, returns cell center
      std::vector<double> IDtoCenter(int id);
      /// Given an IJ index, returns global index for a certain depth
      int IJToGlobal(int x_index, int y_index, int depth);
      /// Given a cell and an edge, returns (i, j) index of neighboring cell
      std::array<int, 2> NeighborCell(int i_index, int j_index, int edge);
      /// Given a cell and an edge, returns center of neighboring cell
      std::vector<double> IJToCenter(int i_index, int j_index, int depth);
      /// Given a cell edge and normal returns what pair to interpolate with to find partial edge center
      std::array<std::array<int, 2>, 2> InterpolationPair(int i, int j, double nx, double ny, int side_index);
      int MaxSolverDepth();
      double InitialCellSize();
      std::vector<double> Mins();
      std::vector<double> Maxes();
      double XMax();
      double XMin();
      double YMax();
      double YMin();

    private:

      void SetupMesh_(double cell_size, double y_min, double y_max, double x_min, double x_max);
      void CalculateMoments_(int key, double cell_size);
      void RecursiveCalculateMoments_(int key, double cell_size);
      void PropagateUp_();
      double DIntegral_(double beginning,
                        double end,
                        std::array<int, 2> q,
                        int index,
                        double fixed_value);
      double CalcD_(double bd_length,
                    double fixed_value,
                    std::vector<double> cell_center,
                    std::array<int, 2> q,
                    int d,
                    std::array<int, 2> which_d,
                    double cell_size);
      std::map<int, geo_info> boundary_cells_;
      std::map<int, int> cell_map_;
      std::map<int, std::vector<double>> id_to_center_;
      boundary::inputs::GeometryInputBase* input_;
      int Sgn_(double v);
      std::array<int, 2> ProjectedNormal_(int side_index, double nx, double ny);
      int Parity_(int side_index);
      int Q_;
      double initial_cell_size_;
      std::vector<double> mins_;
      std::vector<double> maxes_;
      int max_depth_;
      int max_solver_depth_;
  };
}

}
#endif // EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H
