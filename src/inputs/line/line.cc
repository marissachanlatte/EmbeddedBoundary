#include "inputs/line/line.h"

#include <array>

//! Defines the geometry for the boundary line x = y

namespace boundary {

namespace inputs {

double LineGeometry::BoundaryFunction(double x_value){
  return x_value;
};

double LineGeometry::BoundaryDerivatives(double x_value[2], int degree[2]){
  if (degree[0] == 0 && degree[1] == 0){
    return x_value[1] - x_value[0];
  }

  else if (degree[0] + degree[1] == 1){
    return degree[1] - degree[0];
  }

  else {
    return 0;
  }
};

double LineGeometry::BoundaryInverse(double y_value){
  return y_value;
};


} // namespace inputs

} // namespace boundary
