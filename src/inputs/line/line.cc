#include "inputs/line/line.h"
#include "helpers/geometry_objects.h"

#include <vector>


//! Defines the geometry for the boundary line x = y

namespace boundary {

namespace inputs {

double LineGeometry::BoundaryFunction(double x_value){
  return x_value;
};

double LineGeometry::BoundaryDerivatives(helpers::Point a_point,
                                         std::vector<int> degree){
  if (degree[0] == 0 && degree[1] == 0){
    return a_point.y_val - a_point.x_val;
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
