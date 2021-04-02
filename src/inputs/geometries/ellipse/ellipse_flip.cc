#include "inputs/geometries/ellipse/ellipse_flip.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <iostream>


namespace boundary {

namespace inputs {

std::vector<double> EllipseFlipGeometry::BoundaryFunction(double x_value){
  std::vector<double> y_values;
  y_values.push_back(std::sqrt(1 - 2*std::pow(x_value, 2)));
  y_values.push_back(-std::sqrt(1 - 2*std::pow(x_value, 2)));
  return y_values;
};

double EllipseFlipGeometry::BoundaryDerivatives(helpers::Point a_point,
                                         std::vector<int> degree){
  if (degree[0] == 0 && degree[1] == 0){
    return 2*std::pow(a_point.value(0), 2) + std::pow(a_point.value(1), 2) - 1;
  }

  else if (degree[0] == 0 && degree[1] == 1){ 
    return 2*a_point.value(1);
  }

  else if (degree[0] == 1 && degree[1] == 0){
    return 4*a_point.value(0);
  }

  else if (degree[0] == 0 && degree[1] == 2){
    return 2;
  }

  else if (degree[0] == 2 && degree[1] == 0){
    return 4;
  }

  else {
    return 0;
  }
};


std::vector<double> EllipseFlipGeometry::BoundaryInverse(double y_value){
  std::vector<double> x_values;
  x_values.push_back(std::sqrt((1 - std::pow(y_value, 2))/2));
  x_values.push_back(-std::sqrt((1 - std::pow(y_value, 2))/2));
  return x_values;
};


int EllipseFlipGeometry::Inside(std::array<double, 2> point){
  if ((2*std::pow(point[0], 2) + std::pow(point[1], 2)) > 1){
    return 0;
  }
  else if ((2*std::pow(point[0], 2) + std::pow(point[1], 2)) < 1){
    return 1;
  }
  else {
    return 2;
  }
};


double EllipseFlipGeometry::XMin(){
  return -1.0; // domain min
};


double EllipseFlipGeometry::XMax(){
  return 1.0; // domain max
};


double EllipseFlipGeometry::YMin(){
  return -1.0; // domain min
};


double EllipseFlipGeometry::YMax(){
  return 1.0; // domain max
};


int EllipseFlipGeometry::MaxDepth() {
  return 3;
};


int EllipseFlipGeometry::MaxSolverDepth() {
  return 3;
};


int EllipseFlipGeometry::QOrder(){
  return 1; // Q order
}
} // namespace inputs

} // namespace boundary
