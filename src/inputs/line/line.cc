#include "inputs/line/line.h"

//! Defines the geometry for the boundary line x = y

namespace boundary {

namespace inputs {

double LineGeometry::BoundaryFunction(double x_value){
  return x_value;
};

double LineGeometry::BoundaryDerivatives(double x_value, int degree){
  if (degree == 0){
    return x_value;
  }

  else if (degree == 1){
    return 1;
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
