#include "inputs/equations/cos_sin_neumann.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <iostream>


namespace boundary {

namespace inputs {

double CosSinNeumann::NeumannCondition(double x_value, double y_value){
  return 0;
};


double CosSinNeumann::RightHandSide(double x_value, double y_value){
  return -std::cos(M_PI*x_value)*std::sin(M_PI*y_value);;
};

} // namespace inputs

} // namespace boundary
