#ifndef EMBEDDED_BOUNDARY_SOLVERS_LAPLACIAN_H
#define EMBEDDED_BOUNDARY_SOLVERS_LAPLACIAN_H

#include "geometry/boundary.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Dense>

namespace boundary {

namespace solvers {

/// A Class Computing the Laplacian

  class Laplacian{
    public:
      Laplacian();
      Laplacian(boundary::geometry::Boundary geometry);
      ~Laplacian() = default;
      Eigen::VectorXd solve();
    
    private:
        Eigen::SparseMatrix<double> matrix_;
        Eigen::VectorXd rhs_;
        void BuildMatrix(boundary::geometry::Boundary geometry);
        void SafeMatrixAssign(int i_index, int j_index, double value);
        int num_x_;
        int num_y_;
  };

} // namespace solvers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SOLVERS_LAPLACIAN_H