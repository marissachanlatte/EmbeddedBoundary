#ifndef EMBEDDED_BOUNDARY_NORMALS_NORMALS_H_
#define EMBEDDED_BOUNDARY_NORMALS_NORMALS_H_

#include "inputs/input_base.h"

namespace boundary {

namespace normals {

class Normal{
  public:
    ~Normal() = default;
    double * ComputeNormal(double x_value[2],
                           boundary::inputs::InputBase* input);
    double NormalizedGradient(double point[2],
                              boundary::inputs::InputBase* input);
    double PartialNormalizedGradient(int P[2], double point[2],
                                     double (*derivative)(double, int));
    double NormalDerivative(int P[2], int dim, double point[2],
                            double (*derivative)(double, int));
};

} // namespace normals

} // namespace boundary


#endif
