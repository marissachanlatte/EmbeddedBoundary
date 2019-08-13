#include "geometry/build_mesh.h"
#include "helpers/geometry_objects.h"
#include "inputs/input_base.h"

#include <iostream>
namespace boundary {

namespace geometry {

MeshTree::MeshTree(boundary::inputs::InputBase* input){
  int level = 0;
  double right = input->Maximum();
  double left = input->Minimum();
  double top = input->Maximum();
  double bottom = input->Minimum();

  num_nodes_ = 0;
  depth_ = input->Depth();
  q_order_ = input->QOrder();

  helpers::QuadTree mesh_tree = MeshTree::BuildMesh(level, right, left, top, bottom);
  mesh_ = &mesh_tree;
  helpers::Point* center = mesh_->GetCellCenter();
  double x = center->x_val;
  std::cout << "mesh_" << x << std::endl;
};

helpers::QuadTree MeshTree::BuildMesh(int level, double right, double left,
                                      double top, double bottom){
  helpers::Point center = helpers::Point((right + left)/2, (top + bottom)/2);
  if (level == depth_){
    helpers::QuadTree mini_tree = helpers::QuadTree(&center, q_order_,
                                                    num_nodes_);
    num_nodes_ += 1;
    return mini_tree;
  }

  helpers::QuadTree north_west = MeshTree::BuildMesh(level + 1, center.x_val,
                                                     left, top, center.y_val);
  helpers::QuadTree north_east = MeshTree::BuildMesh(level + 1, right,
                                                     center.x_val, top,
                                                     center.y_val);
  helpers::QuadTree south_west = MeshTree::BuildMesh(level + 1, center.x_val,
                                                     left, center.y_val,
                                                     bottom);
  helpers::QuadTree south_east = MeshTree::BuildMesh(level + 1, right,
                                                     center.x_val, center.y_val,
                                                     bottom);

  helpers::QuadTree full_tree = helpers::QuadTree(&north_west, &north_east,
                                                  &south_west, &south_east,
                                                  num_nodes_);
  num_nodes_ += 1;
  return full_tree;
};


int MeshTree::NumNodes(){
  return num_nodes_;
};


helpers::QuadTree* MeshTree::GetMesh(){
  return mesh_;
};


} // namespace geometry

} // namespace boundary
