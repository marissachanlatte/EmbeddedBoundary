#include "normals/normals.h"
#include "inputs/input_base.h"
#include "helpers/math_helpers.cc"

#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

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

double Normal::NormalizedGradient(helpers::Point a_point,
                          boundary::inputs::InputBase* input){

  std::vector<int> x_unit{1, 0};
  std::vector<int> y_unit{0, 1};

  double dx = input->BoundaryDerivatives(a_point, x_unit);
  double dy = input->BoundaryDerivatives(a_point, y_unit);
  return std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
};

double Normal::PartialNormalizedGradient(std::vector<int> p_order,
                                 helpers::Point a_point,
                                 boundary::inputs::InputBase* input){
  double partial_l = 0;
  for (int dim = 1; dim < 3; dim++){ // hard coded for DIM == 2
    std::vector<int> unit {0, 0};
    unit[dim - 1] = 1;
    std::vector<int> q_order;
    std::transform(p_order.begin(), p_order.end(), unit.begin(),
                   q_order.begin(), std::minus<int>());
    for (int r1 = 0; r1 < q_order[0] + 1; r1++){
      for (int r2 = 0; r2 < q_order[1] + 1; r2++){
        // Loop structure hard coded for DIM = 2
        std::vector<int> r_order{r1, r2};
        int binom = helpers::MultiIndexBinomial(q_order, r_order);
        // Make new order vector
        std::vector<int> new_order;
        // subtract Q and R
        std::transform(q_order.begin(), q_order.end(), r_order.begin(),
                       new_order.begin(), std::minus<int>());
        // multiply unit vector
        std::vector<int> two_unit;
        int N = 2;
        std::transform(unit.begin(), unit.end(), two_unit.begin(),
                       [N](int i){return i*N;});
        // all together
        std::transform(new_order.begin(), new_order.end(), two_unit.begin(),
                       new_order.begin(), std::plus<int>());
        partial_l += binom
                     * Normal::NormalDerivative(r_order, dim, a_point,input)
                     * input->BoundaryDerivatives(a_point, new_order);
      }
    }
  }
  return partial_l;
};


double Normal::NormalDerivative(std::vector<int> p_order, int dim,
                        helpers::Point a_point,
                        boundary::inputs::InputBase* input){
  return 0;
};

} // namespace normals

} // namespace boundary
