#include "inputs/circle/circle_test.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <iostream>


namespace boundary {

namespace inputs {

std::vector<double> CircleTestGeometry::BoundaryFunction(double x_value){
  std::vector<double> y_values;
  y_values.push_back(std::sqrt(1 - std::pow(x_value, 2)));
  y_values.push_back(-std::sqrt(1 - std::pow(x_value, 2)));
  return y_values;
};


double CircleTestGeometry::BoundaryDerivatives(helpers::Point a_point,
                                         std::vector<int> degree){
  if (degree[0] == 0 && degree[1] == 0){
    return std::pow(a_point.value(1), 2) + std::pow(a_point.value(0), 2) - 1;
  }

  else if (degree[0] + degree[1] == 1){
    return 2*(degree[0]*a_point.value(0) + degree[1]*a_point.value(1));
  }

  else if ((degree[0] == 2 && degree[1] == 0) || (degree[0] == 0 && degree[1] == 2)){
    return 2;
  }

  else {
    return 0;
  }
};


std::vector<double> CircleTestGeometry::BoundaryInverse(double y_value){
  std::vector<double> x_values;
  x_values.push_back(std::sqrt(1 - std::pow(y_value, 2)));
  x_values.push_back(-std::sqrt(1 - std::pow(y_value, 2)));
  return x_values;
};


int CircleTestGeometry::Inside(std::array<double, 2> point){
  if ((std::pow(point[0], 2) + std::pow(point[1], 2)) > 1){
    return 0;
  }
  else if ((std::pow(point[0], 2) + std::pow(point[1], 2)) < 1){
    return 1;
  }
  else {
    return 2;
  }
};


double CircleTestGeometry::XMin(){
  return -1.0; // domain min
};


double CircleTestGeometry::XMax(){
  return 1.0; // domain max
};


double CircleTestGeometry::YMin(){
  return -1.0; // domain min
};


double CircleTestGeometry::YMax(){
  return 1.0; // domain max
};


double CircleTestGeometry::CellSize(){
  return .25;
};

int CircleTestGeometry::QOrder(){
  return 1; // Q order
}
} // namespace inputs

} // namespace boundary
