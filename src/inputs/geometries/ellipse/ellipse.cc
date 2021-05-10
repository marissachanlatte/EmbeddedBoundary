#include "inputs/geometries/ellipse/ellipse.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <iostream>


namespace boundary {

namespace inputs {

std::vector<double> EllipseGeometry::BoundaryFunction(double x_value){
  std::vector<double> y_values;
  y_values.push_back(std::sqrt((1 - std::pow(x_value, 2))/2));
  y_values.push_back(-std::sqrt((1 - std::pow(x_value, 2))/2));
  return y_values;
};

double EllipseGeometry::BoundaryDerivatives(helpers::Point a_point,
                                         std::vector<int> degree){
  if (degree[0] == 0 && degree[1] == 0){
    return 2*std::pow(a_point.val(1), 2) + std::pow(a_point.val(0), 2) - 1;
  }

  else if (degree[0] == 1 && degree[1] == 0){ 
    return 2*a_point.val(0);
  }

  else if (degree[0] == 0 && degree[1] == 1){
    return 4*a_point.val(1);
  }

  else if (degree[0] == 2 && degree[1] == 0){
    return 2;
  }

  else if (degree[0] == 0 && degree[1] == 2){
    return 4;
  }

  else {
    return 0;
  }
};


std::vector<double> EllipseGeometry::BoundaryInverse(double y_value){
  std::vector<double> x_values;
  x_values.push_back(std::sqrt(1 - 2*std::pow(y_value, 2)));
  x_values.push_back(-std::sqrt(1 - 2*std::pow(y_value, 2)));
  return x_values;
};


int EllipseGeometry::Inside(helpers::Point a_point){
  if ((std::pow(a_point.val(0), 2) + 2*std::pow(a_point.val(1), 2)) > 1){
    return 0;
  }
  else if ((std::pow(a_point.val(0), 2) + 2*std::pow(a_point.val(1), 2)) < 1){
    return 1;
  }
  else {
    return 2;
  }
};


double EllipseGeometry::XMin(){
  return -1.0; // domain min
};


double EllipseGeometry::XMax(){
  return 1.0; // domain max
};


double EllipseGeometry::YMin(){
  return -1.0; // domain min
};


double EllipseGeometry::YMax(){
  return 1.0; // domain max
};


int EllipseGeometry::MaxDepth() {
  return depth;
};


int EllipseGeometry::MaxSolverDepth() {
  return depth;
};


int EllipseGeometry::QOrder(){
  return 1; // Q order
}
} // namespace inputs

} // namespace boundary
