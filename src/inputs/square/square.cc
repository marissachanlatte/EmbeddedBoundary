#include "inputs/square/square.h"
#include "helpers/geometry_objects.h"

#include <vector>


namespace boundary {

namespace inputs {

std::vector<double> SquareGeometry::BoundaryFunction(double x_value){
  std::vector<double> y_values;
  y_values.push_back(0.0);
  return y_values;
};


double SquareGeometry::BoundaryDerivatives(helpers::Point a_point,
                                         std::vector<int> degree){
  return 0;
};


std::vector<double> SquareGeometry::BoundaryInverse(double y_value){
  std::vector<double> x_values;
  x_values.push_back(0.0);
  return x_values;
};


int SquareGeometry::Inside(std::array<double, 2> point){
  if (point[0] < 1.0 && point[0] > -1.0 && point[1] < 1.0 && point[1] > -1.0){
    return 1;
  }
  else if ((point[0] == 1.0 && point[1] <= 1.0 && point[1] >= -1.0) ||
           (point[1] == 1.0 && point[0] <= 1.0 && point[0] >= -1.0) || 
           (point[1] == -1.0 && point[0] <= 1.0 && point[0] >= -1.0) || 
           (point[0] == -1.0 && point[1] <= 1.0 && point[1] >= -1.0)){
    return 2;
  }
  else {
    return 0;
  }
};


double SquareGeometry::XMin(){
  return -2.0; // domain min
};


double SquareGeometry::XMax(){
  return 2.0; // domain max
};


double SquareGeometry::YMin(){
  return -2.0; // domain min
};


double SquareGeometry::YMax(){
  return 2.0; // domain max
};


double SquareGeometry::InitialCellSize(){
  return 1.0;
};

int SquareGeometry::QOrder(){
  return 1; // Q order
}
} // namespace inputs

} // namespace boundary
