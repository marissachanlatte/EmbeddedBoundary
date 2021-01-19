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
      /// Constructor for laplacian, takes in a geometry.
      Laplacian(boundary::geometry::Boundary geometry);
      /// Default destructor
      ~Laplacian() = default;
      /// Solves linear system
      Eigen::VectorXd solve();
      
    
    private:
        /// Left hand side of linear system
        Eigen::SparseMatrix<double> matrix_;
        /// Right hand side of linear system
        Eigen::VectorXd rhs_;
        /// Builds appropriate matrix given geometry
        void BuildMatrix(boundary::geometry::Boundary geometry);
        /// Assigns value to matrix, skipping if entry doesn't exist
        void SafeMatrixAssign(int i_index, int j_index, double value);
        /// Number of cells in x direction
        int num_x_;
        /// Number of cells in y direction
        int num_y_; 
  };

} // namespace solvers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SOLVERS_LAPLACIAN_H