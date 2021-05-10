#ifndef EMBEDDED_BOUNDARY_SOLVERS_LAPLACE_OPERATOR_HO_H
#define EMBEDDED_BOUNDARY_SOLVERS_LAPLACE_OPERATOR_HO_H

#include "geometry/boundary.h"

#include <vector>
#include <Eigen/Dense>

namespace boundary {

namespace solvers {

/// A Class Computing the Laplace Operator

    class LaplaceOperatorHO{
        public:
            LaplaceOperatorHO(boundary::geometry::Boundary geometry);
            std::vector<double> Neighborhood(std::vector<double> cell_center, double cell_size);
            void ComputeAndPrint();
        private:
            double NeumannCondition(std::vector<double> point, std::vector<double> normal);
            double Phi(std::vector<double> point);
            Eigen::MatrixXf ComputeM(std::vector<double> cell_center);
            double depth_;
            boundary::geometry::Boundary geometry_;
            std::map<double, int> cell_map_;
    };
} // namespace solvers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SOLVERS_LAPLACE_OPERATOR_HO_H