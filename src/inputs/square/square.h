#ifndef EMBEDDED_BOUNDARY_SRC_INPUTS_SQUARE_SQUARE_H
#define EMBEDDED_BOUNDARY_SRC_INPUTS_SQUARE_SQUARE_H

#include "inputs/input_base.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {


/// An input file representing the square [-1, 1] x [-1, 1]
class SquareGeometry : public InputBase{
  public:
    ~SquareGeometry() = default;
    std::vector<double> BoundaryFunction(double x_value) override;
    double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) override;
    std::vector<double> BoundaryInverse(double y_value) override;
    int Inside(std::array<double, 2> point) override;
    double XMin() override;
    double XMax() override;
    double YMin() override;
    double YMax() override;
    // double InitialCellSize() override;
    int MaxDepth() override;
    int MaxSolverDepth() override;
    int QOrder() override;

};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SRC_INPUTS_SQUARE_SQUARE_H
