#include "inputs/geometries/circle/circle_shift.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <iostream>


namespace boundary {

namespace inputs {

// Shift Circle
double y_shift = 0;
double x_shift = 0;
double radius_sqd = 2;

std::vector<double> CircleShiftGeometry::BoundaryFunction(double x_value){
  std::vector<double> y_values;
  y_values.push_back(std::sqrt(radius_sqd - std::pow(x_value + x_shift, 2)) - y_shift);
  y_values.push_back(-std::sqrt(radius_sqd - std::pow(x_value + x_shift, 2)) - y_shift);
  return y_values;
};

double CircleShiftGeometry::BoundaryDerivatives(helpers::Point a_point,
                                         std::vector<int> degree){
  if (degree[0] == 0 && degree[1] == 0){
    return std::pow(a_point.val(0) + x_shift, 2) + std::pow(a_point.val(1) + y_shift, 2) - radius_sqd;
  }

  else if (degree[0] + degree[1] == 1){
    return 2*(degree[0]*(a_point.val(0) + x_shift) + degree[1]*(a_point.val(1) + y_shift));
  }

  else if ((degree[0] == 2 && degree[1] == 0) || (degree[0] == 0 && degree[1] == 2)){
    return 2;
  }

  else {
    return 0;
  }
};


std::vector<double> CircleShiftGeometry::BoundaryInverse(double y_value){
  std::vector<double> x_values;
  x_values.push_back(std::sqrt(radius_sqd - std::pow(y_value + y_shift, 2)) - x_shift);
  x_values.push_back(-std::sqrt(radius_sqd - std::pow(y_value + y_shift, 2)) - x_shift);
  return x_values;
};


int CircleShiftGeometry::Inside(helpers::Point a_point){
  if ((std::pow(a_point.val(0) + x_shift, 2) + std::pow(a_point.val(1) + y_shift, 2)) > radius_sqd){
    return 0;
  }
  else if ((std::pow(a_point.val(0) + x_shift, 2) + std::pow(a_point.val(1) + y_shift, 2)) < radius_sqd){
    return 1;
  }
  else {
    return 2;
  }
};


double CircleShiftGeometry::XMin(){
  return -x_shift-radius_sqd; // domain min
};


double CircleShiftGeometry::XMax(){
  return -x_shift+radius_sqd; // domain max
};


double CircleShiftGeometry::YMin(){
  return -y_shift-radius_sqd; // domain min
};


double CircleShiftGeometry::YMax(){
  return -y_shift+radius_sqd; // domain max
};


int CircleShiftGeometry::MaxDepth() {
  return depth;
};


int CircleShiftGeometry::MaxSolverDepth() {
  return depth;
};


int CircleShiftGeometry::QOrder(){
  return 1; // Q order
}
} // namespace inputs

} // namespace boundary
