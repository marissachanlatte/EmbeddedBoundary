#include "inputs/line/line_zero.h"
#include "helpers/geometry_objects.h"

#include <vector>



//! Defines the geometry for the boundary line x = y

namespace boundary {

namespace inputs {

double LineGeometryZero::BoundaryFunction(double x_value){
  return x_value;
};

double LineGeometryZero::BoundaryDerivatives(helpers::Point a_point,
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


double LineGeometryZero::BoundaryInverse(double y_value){
  return y_value;
};

double LineGeometryZero::Minimum(){
  return -1.0; // domain min
};


double LineGeometryZero::Maximum(){
  return 1.0; // domain max
};


int LineGeometryZero::Depth(){
  return 0; // quad tree depth, num cells = 4^DEPTH
};

int LineGeometryZero::QOrder(){
  return 1; // Q order
}
} // namespace inputs

} // namespace boundary
