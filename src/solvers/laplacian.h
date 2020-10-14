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
        void safeMatrixAssign(int i_index, int j_index, double value);
        int IJToGlobal(int x_index, int y_index, int num_x);
        std::array<int, 2> neighborCell(int i_index, int j_index, int edge);
        int sgn(double v);
        std::array<int, 2> projected_normal(int side_index, double nx, double ny);
        int parity(int side_index);
        std::array<std::array<int, 2>, 2> interpolationPair(int i, int j, double nx, double ny, int side_index);
        int num_x_;
        int num_y_;

  };

} // namespace solvers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SOLVERS_LAPLACIAN_H