#ifndef EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
#define EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H

#include <vector>
#include <cmath>

namespace boundary {

namespace helpers {

  class Point{
    public:
      Point();
      Point(double first_dim, double second_dim);
      ~Point() = default;
      Point operator + (const Point &a_point);
      double x_val;
      double y_val;
  };

} // namespace helpers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
