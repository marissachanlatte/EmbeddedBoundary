#ifndef EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H
#define EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H

namespace boundary {

namespace inputs {

class InputBase {
public:

  virtual ~InputBase() = default;
  virtual double BoundaryFunction(double x_value) = 0;
  virtual double BoundaryDerivatives(double x_value[2], int degree[2]) = 0;
  virtual double BoundaryInverse(double y_value) = 0;
};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H
