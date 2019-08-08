#ifndef EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H
#define EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H

#include "helpers/geometry_objects.h"

#include <vector>

namespace boundary {

namespace inputs {

class InputBase {
public:
  virtual ~InputBase() = default;
  virtual double BoundaryFunction(double x_value) = 0;
  virtual double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) = 0;
  virtual double BoundaryInverse(double y_value) = 0;
  virtual double Minimum() = 0;
  virtual double Maximum() = 0;
  virtual int Depth() = 0;

};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H
