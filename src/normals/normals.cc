#include "normals/normals.h"
#include "inputs/input_base.h"

#include <cmath>

namespace boundary {

namespace normals {

double * Normal::ComputeNormal(double x_value[2],
                               boundary::inputs::InputBase* input){
  int x_unit[2] = {1, 0};
  int y_unit[2] = {0, 1};
  double dx = input->BoundaryDerivatives(x_value, x_unit);
  double dy = input->BoundaryDerivatives(x_value, y_unit);
  double normalization = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
  static double gradient[2] = {dx/normalization, dy/normalization};
  return gradient;
}

double NormalizedGradient(double point[2], boundary::inputs::InputBase* input){
  int x_unit[2] = {1, 0};
  int y_unit[2] = {0, 1};
  double dx = input->BoundaryDerivatives(point, x_unit);
  double dy = input->BoundaryDerivatives(point, y_unit);
  return std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
};

// double PartialNormalizedGradient(int P[2], double point[2],
//                                  double (*derivative)(double, int));
// double NormalDerivative(int P[2], int dim, double point[2],
//                         double (*derivative)(double, int));

} // namespace normals

} // namespace boundary
