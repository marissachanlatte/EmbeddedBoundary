#ifndef EMBEDDED_BOUNDARY_NORMALS_NORMALS_H_
#define EMBEDDED_BOUNDARY_NORMALS_NORMALS_H_

#include "inputs/input_base.h"
#include "helpers/geometry_objects.h"

namespace boundary {

namespace normals {

class Normal{
  public:
    ~Normal() = default;
    double * ComputeNormal(helpers::Point a_point,
                           boundary::inputs::InputBase* input);
    double NormalizedGradient(helpers::Point a_point,
                              boundary::inputs::InputBase* input);
    double PartialNormalizedGradient(std::vector<int> p_order,
                                     helpers::Point a_point,
                                     boundary::inputs::InputBase* input);
    double NormalDerivative(std::vector<int> p_order, int dim,
                            helpers::Point a_point,
                            boundary::inputs::InputBase* input);
};

} // namespace normals

} // namespace boundary


#endif
