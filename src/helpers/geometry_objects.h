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

  class QuadTree{
    public:
      QuadTree();
      QuadTree(QuadTree* north_west,
               QuadTree* north_east,
               QuadTree* south_west,
               QuadTree* south_east,
               int id);
      QuadTree(Point* cell_center, int degree, int id);
      void AssignNormals(std::vector<double> value, std::vector<int> degree);
      void AssignVolume(double value, std::vector<int> degree);
      void AssignBoundary(double value, std::vector<int> degree);

      Point* GetCellCenter();
      std::vector<double> GetNormal(std::vector<int> degree);
      double GetVolume(std::vector<int> degree);
      double GetBoundary(std::vector<int> degree);
      int GetID();

      QuadTree* NorthWest();
      QuadTree* NorthEast();
      QuadTree* SouthWest();
      QuadTree* SouthEast();

    private:
      // Children
      QuadTree* north_west_;
      QuadTree* north_east_;
      QuadTree* south_west_;
      QuadTree* south_east_;

      // Cell Center
      Point* cell_center_;
      int degree_;
      int id_;

      // Normals
      std::vector<std::vector<std::vector<double>>> normals_;

      // Moments
      std::vector<std::vector<double>> volume_moments_;
      std::vector<std::vector<double>> boundary_moments_;

  };

} // namespace helpers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
