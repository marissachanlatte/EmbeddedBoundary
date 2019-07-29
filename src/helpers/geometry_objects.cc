# include "helpers/geometry_objects.h"

namespace boundary {

namespace helpers {

Point::Point(){
  x_val = y_val = 0;
};

Point::Point(double first_dim, double second_dim){
  x_val = first_dim;
  y_val = second_dim;
};

Point Point::operator + (const Point &a_point){
  Point result;
  result.x_val = x_val + a_point.x_val;
  result.y_val = y_val + a_point.y_val;
  return result;
};

} // namespace helpers

} // namespace boundary
