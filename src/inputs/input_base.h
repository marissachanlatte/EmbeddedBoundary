#ifndef EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H
#define EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H

#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {

class InputBase {
public:
  virtual ~InputBase() = default;
  virtual double BoundaryFunction(double x_value) = 0;
  virtual double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) = 0;
  virtual double BoundaryInverse(double y_value) = 0;
  virtual int Inside(std::array<double, 2> point) = 0; // returns 0 for outside, 1 for inside, 2 for on boundary
  virtual double XMin() = 0;
  virtual double XMax() = 0;
  virtual double YMin() = 0;
  virtual double YMax() = 0;
  virtual double CellSize() = 0;
  virtual int QOrder() = 0;

};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H
