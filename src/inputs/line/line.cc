#include "inputs/line/line.h"
#include "helpers/geometry_objects.h"

#include <vector>


namespace boundary {

namespace inputs {

std::vector<double> LineGeometry::BoundaryFunction(double x_value){
  std::vector<double> y_values;
  y_values.push_back(x_value);
  return y_values;
};

double LineGeometry::BoundaryDerivatives(helpers::Point a_point,
                                         std::vector<int> degree){
  if (degree[0] == 0 && degree[1] == 0){
    return a_point.value(1) - a_point.value(0);
  }

  else if (degree[0] + degree[1] == 1){
    return degree[1] - degree[0];
  }

  else {
    return 0;
  }
};


std::vector<double> LineGeometry::BoundaryInverse(double y_value){
  std::vector<double> x_values;
  x_values.push_back(y_value);
  return x_values;
};


int LineGeometry::Inside(std::array<double, 2> point){
  if (point[0] < point[1]){
    return 0;
  }
  else if (point[0] > point[1]){
    return 1;
  }
  else {
    return 2;
  }
};


double LineGeometry::XMin(){
  return -1.0; // domain min
};


double LineGeometry::XMax(){
  return 1.0; // domain max
};


double LineGeometry::YMin(){
  return -1.0; // domain min
};


double LineGeometry::YMax(){
  return 1.0; // domain max
};


double LineGeometry::CellSize(){
  return .25;
};

int LineGeometry::QOrder(){
  return 1; // Q order
}
} // namespace inputs

} // namespace boundary
