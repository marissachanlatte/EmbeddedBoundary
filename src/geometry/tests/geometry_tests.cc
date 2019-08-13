#include "geometry/build_mesh.h"
#include "inputs/line/line.h"
#include "inputs/line/line_zero.h"

#include "gtest/gtest.h"

namespace boundary {

namespace geometry {

TEST(GeometryTest, NumNodes){
  inputs::LineGeometry input;
  MeshTree mesh_tree = MeshTree(&input);
  EXPECT_EQ(341, mesh_tree.NumNodes());
}

TEST(GeometryTest, NumNodesZero){
  inputs::LineGeometryZero input_zero;
  MeshTree mesh_tree_zero = MeshTree(&input_zero);
  EXPECT_EQ(1, mesh_tree_zero.NumNodes());
}
// 
// TEST(GeometryTest, CellCenterZero){
//   inputs::LineGeometryZero input_zero;
//   MeshTree mesh_tree_zero = MeshTree(&input_zero);
//   helpers::Point* center_zero = mesh_tree_zero.GetMesh()->GetCellCenter();
//   double x = center_zero->x_val;
//   std::cout << "x: " << x << std::endl;
//   double y = center_zero->y_val;
//   EXPECT_EQ(0, x);
//   EXPECT_EQ(0, y);
// }
// TEST(GeometryTest, MeshStructure){
//   inputs::LineGeometry input;
//   MeshTree mesh_tree = MeshTree(&input);
//   helpers::QuadTree* full_mesh = mesh_tree.GetMesh();
//   helpers::QuadTree* child1 = full_mesh->NorthWest();
//   EXPECT_FALSE(child1->IsEmpty());
//   helpers::QuadTree* child2 = child1->NorthWest();
//   EXPECT_FALSE(child2->IsEmpty());
//   helpers::QuadTree* child3 = child2->NorthWest();
//   EXPECT_FALSE(child3->IsEmpty());
//   helpers::QuadTree* child4 = child3->NorthWest();
//   helpers::Point* center = child4->GetCellCenter();
//   EXPECT_EQ(-.9375, center->x_val);
//   EXPECT_EQ(.9375, center->y_val);
//
// }

} // namespace geometry

} // namespace boundary
