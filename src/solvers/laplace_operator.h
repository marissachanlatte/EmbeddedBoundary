#ifndef EMBEDDED_BOUNDARY_SOLVERS_LAPLACE_OPERATOR_H
#define EMBEDDED_BOUNDARY_SOLVERS_LAPLACE_OPERATOR_H

#include "geometry/boundary.h"

#include <vector>

namespace boundary {

namespace solvers {

/// A Class Computing the Laplace Operator

    class LaplaceOperator{
        public:
            LaplaceOperator(boundary::geometry::Boundary);
        private:
            double NeumannCondition(std::vector<double> point, std::vector<double> normal);
            double Phi(std::vector<double> point);
    };
} // namespace solvers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SOLVERS_LAPLACE_OPERATOR_H