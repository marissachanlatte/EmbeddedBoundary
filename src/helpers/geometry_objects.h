#ifndef EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
#define EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H

#include <vector>

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

  class QuadTree{
    public:
      QuadTree();
      QuadTree(QuadTree* north_west,
               QuadTree* north_east,
               QuadTree* south_west,
               QuadTree* south_east);
      QuadTree(Point* cell_center, int degree);
      void AssignNormals(double value, int degree);
      void AssignVolume(double value, int degree);
      void AssignBoundary(double value, int degree);

      Point* GetCellCenter();
      double GetNormal(int degree);
      double GetVolume(int degree);
      double GetBoundary(int degree);
    private:
      // Children
      QuadTree* north_west_;
      QuadTree* north_east_;
      QuadTree* south_west_;
      QuadTree* south_east_;

      // Cell Center
      Point* cell_center_;
      int degree_;

      // Normals
      std::vector<double> normals_;

      // Moments
      std::vector<double> volume_moments_;
      std::vector<double> boundary_moments_;

  };
} // namespace helpers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
