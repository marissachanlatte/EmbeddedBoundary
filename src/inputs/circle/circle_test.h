#ifndef EMBEDDED_BOUNDARY_SRC_INPUTS_CIRCLE_CIRCLE_TEST_H
#define EMBEDDED_BOUNDARY_SRC_INPUTS_CIRCLE_CIRCLE_TEST_H

#include "inputs/input_base.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {


/// An input file representing the Circle 0 = x^2 + y^2 - 1
class CircleTestGeometry : public InputBase{
  public:
    ~CircleTestGeometry() = default;
    std::vector<double> BoundaryFunction(double x_value) override;
    double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) override;
    std::vector<double> BoundaryInverse(double y_value) override;
    double NeumannCondition(double x_value, double y_value) override;
    int Inside(std::array<double, 2> point) override;
    double XMin() override;
    double XMax() override;
    double YMin() override;
    double YMax() override;
    int MaxDepth() override;
    int MaxSolverDepth() override;
    int QOrder() override;
};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SRC_INPUTS_CIRCLE_CIRCLE_TEST_H
