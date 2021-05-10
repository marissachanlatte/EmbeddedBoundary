#ifndef EMBEDDED_BOUNDARY_SRC_INPUTS_EQUATIONS_LAPLACE_NEUMANN_H
#define EMBEDDED_BOUNDARY_SRC_INPUTS_EQUATIONS_LAPLACE_NEUMANN_H

#include "inputs/equations/solver_input_base.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {


/// An input file representing the equation f_xx + f_yy = 0 with f'(x, y) = y along the boundary
class LaplaceNeumann : public SolverInputBase{
  public:
    ~LaplaceNeumann() = default;
    double NeumannCondition(double x_value, double y_value) override;
    double RightHandSide(double x_value, double y_value) override;
};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SRC_INPUTS_EQUATIONS_LAPLACE_NEUMANN_H
