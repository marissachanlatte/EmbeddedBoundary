#ifndef EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
#define EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H

#include <vector>
#include <cmath>

namespace boundary {

namespace helpers {

/// A Class describing a point
/** This class represents a d-dimensional point
*/
  class Point{
    public:
      Point();
      /// A constructor using a vector of the coordinates of the point
      Point(std::vector<double> coords);
      ~Point() = default;
      Point operator + (const Point &a_point);
      /// Returns the d value of the point
      double value(int d) const;
    private:
      /// The dimension of the point
      int _dim;
      /// A vector of the coordinates of the point
      std::vector<double> _coords;
  };

} // namespace helpers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
