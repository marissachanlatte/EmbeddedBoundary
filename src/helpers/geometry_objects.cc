# include "helpers/geometry_objects.h"

#include <vector>

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


QuadTree::QuadTree(QuadTree* north_west,
                   QuadTree* north_east,
                   QuadTree* south_west,
                   QuadTree* south_east){
  north_west_ = north_west;
  north_east_ = north_east;
  south_west_ = south_west;
  south_east_ = south_east;
};


QuadTree::QuadTree(Point* cell_center, int degree){
  cell_center_ = cell_center;
  degree_ = degree;
  normals_.resize(degree);
  volume_moments_.resize(degree);
  boundary_moments_.resize(degree);
};


void QuadTree::AssignNormals(double value, int degree){
  normals_[degree] = value;
};


void QuadTree::AssignVolume(double value, int degree){
  volume_moments_[degree] = value;
};


void QuadTree::AssignBoundary(double value, int degree){
  boundary_moments_[degree] = value;
};


Point* QuadTree::GetCellCenter(){
  return cell_center_;
};


double QuadTree::GetNormal(int degree){
  return normals_[degree];
};


double QuadTree::GetVolume(int degree){
  return volume_moments_[degree];
};


double QuadTree::GetBoundary(int degree){
  return boundary_moments_[degree];
};

} // namespace helpers

} // namespace boundary
