#ifndef EMBEDDED_BOUNDARY_SRC_INPUTS_GEOMETRIES_ELLIPSE_ELLIPSE_FLIP_H
#define EMBEDDED_BOUNDARY_SRC_INPUTS_GEOMETRIES_ELLIPSE_ELLIPSE_FLIP_H

#include "inputs/geometries/geometry_input_base.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {


/// An input file representing the Ellipse 0 = 2x^2 + y^2 - 1
class EllipseFlipGeometry : public GeometryInputBase{
  public:
    ~EllipseFlipGeometry() = default;
    std::vector<double> BoundaryFunction(double x_value) override;
    double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) override;
    std::vector<double> BoundaryInverse(double y_value) override;
    int Inside(helpers::Point a_point) override;
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

#endif // EMBEDDED_BOUNDARY_SRC_INPUTS_GEOMETRY_ELLIPSE_ELLIPSE-FLIP_H
