#ifndef EMBEDDED_BOUNDARY_SRC_INPUTS_LINE_LINE_H
#define EMBEDDED_BOUNDARY_SRC_INPUTS_LINE_LINE_H

#include "inputs/input_base.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {


/// An input file representing the line x=y
class LineGeometry : public InputBase{
  public:
    ~LineGeometry() = default;
    double BoundaryFunction(double x_value) override;
    double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) override;
    double BoundaryInverse(double y_value) override;
    int Inside(std::array<double, 2> point) override;
    double XMin() override;
    double XMax() override;
    double YMin() override;
    double YMax() override;
    double CellSize() override;
    int QOrder() override;

};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SRC_INPUTS_LINE_LINE_H
