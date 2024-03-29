#include "normals/normals.h"
#include "inputs/geometries/geometry_input_base.h"
#include "helpers/math_helpers.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

namespace boundary {

namespace normals {

std::array<double, 2> Normal::ComputeNormal(helpers::Point a_point,
                               boundary::inputs::GeometryInputBase* input){

  std::vector<int> x_unit{1, 0};
  std::vector<int> y_unit{0, 1};

  double dx = input->BoundaryDerivatives(a_point, x_unit);
  double dy = input->BoundaryDerivatives(a_point, y_unit);
  double normalization = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
  std::array<double, 2> gradient;
  if (normalization == 0){gradient = {0, 0};}
  else {gradient = {dx/normalization, dy/normalization};}
  return gradient;
}

double Normal::NormalizedGradient(helpers::Point a_point,
                          boundary::inputs::GeometryInputBase* input){

  std::vector<int> x_unit{1, 0};
  std::vector<int> y_unit{0, 1};

  double dx = input->BoundaryDerivatives(a_point, x_unit);
  double dy = input->BoundaryDerivatives(a_point, y_unit);
  return std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
};

double Normal::PartialNormalizedGradient(std::vector<int> p_order,
                                 helpers::Point a_point,
                                 boundary::inputs::GeometryInputBase* input){
  double partial_l = 0;
  for (int dim = 1; dim < 3; dim++){ // hard coded for DIM == 2
    std::vector<int> unit {0, 0};
    unit[dim - 1] = 1;
    std::vector<int> q_order;
    std::transform(p_order.begin(), p_order.end(), unit.begin(),
                   std::back_inserter(q_order), std::minus<int>());
    for (int r1 = 0; r1 < q_order[0] + 1; r1++){
      for (int r2 = 0; r2 < q_order[1] + 1; r2++){
        // Loop structure hard coded for DIM = 2
        std::vector<int> r_order{r1, r2};
        int binom = helpers::MultiIndexBinomial(q_order, r_order);
        // Make new order vector
        std::vector<int> new_order;
        // subtract Q and R
        std::transform(q_order.begin(), q_order.end(), r_order.begin(),
                       std::back_inserter(new_order), std::minus<int>());
        // multiply unit vector
        std::vector<int> two_unit;
        std::transform(unit.begin(), unit.end(), std::back_inserter(two_unit),
                       std::bind(std::multiplies<int>(),
                       std::placeholders::_1, 2));
        // all together
        std::vector<int> sum_order;
        std::transform(new_order.begin(), new_order.end(), two_unit.begin(),
                        std::back_inserter(sum_order), std::plus<int>());
        partial_l += binom
                     * Normal::NormalDerivative(r_order, dim, a_point, input)
                     * input->BoundaryDerivatives(a_point, sum_order);
      }
    }
  }
  return partial_l;
};


double Normal::NormalDerivative(std::vector<int> p_order, int dim,
                        helpers::Point a_point,
                        boundary::inputs::GeometryInputBase* input){
  if (p_order[0] == 0 && p_order[1] == 0){
    return Normal::ComputeNormal(a_point, input)[dim - 1];
  }
  std::vector<int> unit {0, 0};
  unit[dim - 1] = 1;
  double deriv_sums = 0;
  for (int q1 = 0; q1 < (p_order[0] + 1); q1++){
    for (int q2 = 0; q2 < (p_order[1] + 1); q2++){
      std::vector<int> q_order{q1, q2};
      if (q_order[0] == p_order[0] && q_order[1] == p_order[1]){
        continue;
      }
      int binom = helpers::MultiIndexBinomial(p_order, q_order);
      // subtract P and Q
      std::vector<int> new_order;
      std::transform(p_order.begin(), p_order.end(), q_order.begin(),
                     std::back_inserter(new_order), std::minus<int>());
      deriv_sums += binom
                    *Normal::PartialNormalizedGradient(new_order, a_point, input)
                    *Normal::NormalDerivative(q_order, dim, a_point, input);
    }
  }

  // add P and unit
  std::vector<int> unit_added;
  std::transform(p_order.begin(), p_order.end(), unit.begin(),
                 std::back_inserter(unit_added), std::plus<int>());
  double partial_n = input->BoundaryDerivatives(a_point, unit_added) - deriv_sums;
  double normalized_gradient = Normal::NormalizedGradient(a_point, input);
  if (normalized_gradient == 0){return 0;}
  else {return partial_n/Normal::NormalizedGradient(a_point, input);}
};

} // namespace normals

} // namespace boundary
