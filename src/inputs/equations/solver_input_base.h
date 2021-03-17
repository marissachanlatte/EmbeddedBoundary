#ifndef EMBEDDED_BOUNDARY_INPUT_SOLVER_INPUT_BASE_H
#define EMBEDDED_BOUNDARY_INPUT_SOLVER_INPUT_BASE_H

#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {

/// The input file base class.
/**
All input files are based on this class which
contains all relevant input geometry information.
*/
class SolverInputBase {
public:
  virtual ~SolverInputBase() = default;
  /// Returns the value of the derivative at the boundary
  virtual double NeumannCondition(double x_value, double y_value) = 0;
};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_INPUT_SOLVER_INPUT_BASE_H
