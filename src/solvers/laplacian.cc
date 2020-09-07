#include "solvers/laplacian.h"

namespace boundary {

namespace solvers {

Laplacian::Laplacian(boundary::geometry::Boundary geometry){
    matrix_ = BuildMatrix(geometry);
};


Eigen::SparseMatrix<double> Laplacian::BuildMatrix(boundary::geometry::Boundary geometry){
    // Figure out size of matrix
    int num_x = int(std::abs(geometry.XMax() - geometry.XMin())/geometry.CellSize());
    int num_y = int(std::abs(geometry.YMax() - geometry.YMin())/geometry.CellSize());
    Eigen::SparseMatrix<double> matrix(num_x*num_y, num_x*num_y);
    double cell_size = geometry.CellSize();

    // Get cell map
    std::map<int, int> cell_map = geometry.CellMap();
    // Iterate through all cells
    for (int j = 0; j < num_x; j++){ // y-index
      for (int i = 0; i < num_y; i++){ // x-index
        int global_id = IJToGlobal(i, j, num_x);
        int covered_id = cell_map[global_id];
        // If boundary cell
        if (covered_id == 1){

        }
        // If interior, five point stencil
        else if (covered_id == 0){
          matrix.coeffRef(global_id, IJToGlobal(i + 1, j, num_x)) += 1/std::pow(cell_size, 2);
          matrix.coeffRef(global_id, IJToGlobal(i - 1, j, num_x)) += 1/std::pow(cell_size, 2);
          matrix.coeffRef(global_id, IJToGlobal(i, j + 1, num_x)) += 1/std::pow(cell_size, 2);
          matrix.coeffRef(global_id, IJToGlobal(i, j - 1, num_x)) += 1/std::pow(cell_size, 2);
          matrix.coeffRef(global_id, global_id) += 4/std::pow(cell_size, 2);
        }
        // If exterior leave as 0
      }
    }
};


int Laplacian::IJToGlobal(int x_index, int y_index, int num_x){
  return num_x*y_index + x_index;
};


} // namespace solvers

} // namespace boundary