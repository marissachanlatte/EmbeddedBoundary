#ifndef EMBEDDED_BOUNDARY_SOLVERS_LAPLACIAN_H
#define EMBEDDED_BOUNDARY_SOLVERS_LAPLACIAN_H

#include "geometry/boundary.h"

#include <Eigen/Sparse>

namespace boundary {

namespace solvers {

/// A Class Computing the Laplacian

  class Laplacian{
    public:
      Laplacian();
      Laplacian(boundary::geometry::Boundary geometry);
      ~Laplacian() = default;
    
    private:
        Eigen::SparseMatrix<double> matrix_;
        Eigen::SparseMatrix<double> BuildMatrix(boundary::geometry::Boundary geometry);
        int IJToGlobal(int x_index, int y_index, int num_x);
  };

} // namespace solvers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SOLVERS_LAPLACIAN_H