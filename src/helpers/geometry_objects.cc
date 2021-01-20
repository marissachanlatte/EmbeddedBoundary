#include "helpers/geometry_objects.h"

#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

namespace boundary {

namespace helpers {


Point::Point(){
  _coords.push_back(std::numeric_limits<double>::quiet_NaN());
}


Point::Point(std::vector<double> coords){
  _coords = coords;
  _dim = coords.size();
}

double Point::value(int d) const {
  return _coords[d];
}


Point Point::operator + (const Point &a_point){
  std::vector<double> result;
  // Loop over dimension
  for (int d = 0; d < _dim; d++){
    result.push_back(_coords[d] + a_point.value(d));
  }
  return Point(result);
}

} // namespace helpers

} // namespace boundary
