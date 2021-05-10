#include "inputs/equations/laplace_neumann.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <iostream>


namespace boundary {

namespace inputs {

double LaplaceNeumann::NeumannCondition(double x_value, double y_value){
  return y_value;
};


double LaplaceNeumann::RightHandSide(double x_value, double y_value){
  return 0;
};

} // namespace inputs

} // namespace boundary
