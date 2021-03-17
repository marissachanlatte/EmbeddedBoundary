#ifndef EMBEDDED_BOUNDARY_NORMALS_NORMALS_H_
#define EMBEDDED_BOUNDARY_NORMALS_NORMALS_H_

#include "inputs/geometries/geometry_input_base.h"
#include "helpers/geometry_objects.h"

#include<array>

namespace boundary {

namespace normals {

/// A class to compute the normals and derivatives of the normals
class Normal{
  public:
    ~Normal() = default;
    /// Computes the normal to a given boundary at a given point
    static std::array<double, 2> ComputeNormal(helpers::Point a_point,
                           boundary::inputs::GeometryInputBase* input);
    /// Computes the normalized gradient of the boundary at a given point.
    static double NormalizedGradient(helpers::Point a_point,
                              boundary::inputs::GeometryInputBase* input);
    /// Computes the partial gradient of order p_order at a given point.
    static double PartialNormalizedGradient(std::vector<int> p_order,
                                     helpers::Point a_point,
                                     boundary::inputs::GeometryInputBase* input);
    /// Computes the p-th derivative of the normal to the boundary at a given point.
    static double NormalDerivative(std::vector<int> p_order, int dim,
                            helpers::Point a_point,
                            boundary::inputs::GeometryInputBase* input);
};

} // namespace normals

} // namespace boundary


#endif
