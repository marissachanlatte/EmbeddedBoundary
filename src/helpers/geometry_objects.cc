#include "helpers/geometry_objects.h"

#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

namespace boundary {

namespace helpers {

Point::Point(){
  x_val = y_val = std::numeric_limits<double>::quiet_NaN();
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

QuadTree::QuadTree(){
  Point center = Point();
  cell_center_ = &center;
  empty_ = true;
};


QuadTree::QuadTree(QuadTree* north_west,
                   QuadTree* north_east,
                   QuadTree* south_west,
                   QuadTree* south_east,
                   int id){
  north_west_ = north_west;
  north_east_ = north_east;
  south_west_ = south_west;
  south_east_ = south_east;
  id_ = id;
  Point center = Point();
  cell_center_ = &center;
  empty_ = false;
};


QuadTree::QuadTree(Point* cell_center, int degree, int id){
  cell_center_ = cell_center;
  degree_ = degree;
  // resize normals
  normals_.resize(degree + 1);
  volume_moments_.resize(degree + 1);
  boundary_moments_.resize(degree + 1);
  for (int i = 0; i < (degree + 1); i++){
    normals_[i].resize(degree + 1);
    volume_moments_[i].resize(degree + 1);
    boundary_moments_[i].resize(degree + 1);
    for (int j = 0; j < (degree + 1); j++){
      normals_[i][j].resize(2);
    }
  }
  id_ = id;
  empty_= false;
  QuadTree dummy_tree = QuadTree();
  north_west_ = &dummy_tree;
  north_east_ = &dummy_tree;
  south_west_ = &dummy_tree;
  south_east_ = &dummy_tree;
};


void QuadTree::AssignNormals(std::vector<double> value,
                             std::vector<int> degree){
  normals_[degree[0]][degree[1]] = value;
};


void QuadTree::AssignVolume(double value, std::vector<int> degree){
  volume_moments_[degree[0]][degree[1]] = value;
};


void QuadTree::AssignBoundary(double value, std::vector<int> degree){
  boundary_moments_[degree[0]][degree[1]] = value;
};


Point* QuadTree::GetCellCenter(){
  return cell_center_;
};


std::vector<double> QuadTree::GetNormal(std::vector<int> degree){
  return normals_[degree[0]][degree[1]];
};


double QuadTree::GetVolume(std::vector<int> degree){
  return volume_moments_[degree[0]][degree[1]];
};


double QuadTree::GetBoundary(std::vector<int> degree){
  return boundary_moments_[degree[0]][degree[1]];
};

QuadTree* QuadTree::NorthWest(){
  return north_west_;
};

QuadTree* QuadTree::NorthEast(){
  return north_east_;
};

QuadTree* QuadTree::SouthEast(){
  return south_east_;
};

QuadTree* QuadTree::SouthWest(){
  return south_west_;
};

int QuadTree::GetID(){
  return id_;
};

bool QuadTree::IsEmpty(){
  return empty_;
}
} // namespace helpers

} // namespace boundary
