#ifndef EMBEDDED_BOUNDARY_SRC_INPUTS_LINE_LINE_H
#define EMBEDDED_BOUNDARY_SRC_INPUTS_LINE_LINE_H

#include "inputs/input_base.h"
#include "helpers/geometry_objects.h"

#include <vector>

namespace boundary {

namespace inputs {

class LineGeometry : public InputBase{
  public:
    ~LineGeometry() = default;
    double BoundaryFunction(double x_value) override;
    double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) override;
    double BoundaryInverse(double y_value) override;

};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SRC_INPUTS_LINE_LINE_H
