#include "normals/normals.h"
#include "inputs/input_base.h"

#include <cmath>

namespace boundary {

namespace normals {

double * Normal::ComputeNormal(helpers::Point a_point,
                               boundary::inputs::InputBase* input){

  std::vector<int> x_unit{1, 0};
  std::vector<int> y_unit{0, 1};

  double dx = input->BoundaryDerivatives(a_point, x_unit);
  double dy = input->BoundaryDerivatives(a_point, y_unit);
  double normalization = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
  static double gradient[2] = {dx/normalization, dy/normalization};
  return gradient;
}

double NormalizedGradient(helpers::Point a_point,
                          boundary::inputs::InputBase* input){

  std::vector<int> x_unit{1, 0};
  std::vector<int> y_unit{0, 1};

  double dx = input->BoundaryDerivatives(a_point, x_unit);
  double dy = input->BoundaryDerivatives(a_point, y_unit);
  return std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
};

// double PartialNormalizedGradient(int P[2], double point[2],
//                                  boundary::inputs::InputBase* input){
//   double partial_l = 0;
//   for (int dim=1; dim < 3; dim++){ // hard coded for DIM == 2
//     int e[2] = {0, 0};
//     e[dim - 1] = 1;
//     int Q[2];
//     std::transform(P, P + 2, e, Q, std::minus<int>());
//     for (int r1 = 0; r1 < Q[0] + 1; r1++){
//       for (int r2 = 0; r2 < Q[1] + 1; r2++){
//         // Loop structure hard coded for DIM = 2
//         int R[2] = {r1, r2};
//
//       }
//     }
//
//   }
//                                  };
// double NormalDerivative(int P[2], int dim, double point[2],
//                         double (*derivative)(double, int));

} // namespace normals

} // namespace boundary
