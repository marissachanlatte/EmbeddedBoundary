#include "inputs/line/line.h"
#include "helpers/geometry_objects.h"
#include "geometry/build_mesh.h"
#include "normals/normals.h"

#include <stack>
#include <vector>
#include <cmath>
#include <iostream>

int main(){
  // Setup Input
  boundary::inputs::LineGeometry input;
  int q_order = input.QOrder();

  // Build Mesh
  boundary::geometry::MeshTree mesh_tree = boundary::geometry::MeshTree(&input);

  // Build Normals
  for (int q_mag = q_order; q_mag > -1; q_mag--){
    for (int q1 = 0; q1 < (q_mag + 1); q1 ++){
      std::vector<int> q{q1, q_mag - q1};
      // Initialize Stack
      std::stack<boundary::helpers::QuadTree*> stack;
      int num_nodes = mesh_tree.NumNodes();
      std::vector<bool> discovered(num_nodes);
      boundary::helpers::QuadTree* mesh = mesh_tree.GetMesh();
      stack.push(mesh);
      while (!stack.empty()){
        boundary::helpers::QuadTree* node = stack.top();
        stack.pop();
        if (!discovered[node->GetID()]){
          // Mark as Discovered
          discovered[node->GetID()] = true;
          std::cout << "test 3" << std::endl;
          // Identify Children
          boundary::helpers::QuadTree* north_west = node->NorthWest();
          boundary::helpers::QuadTree* north_east = node->NorthEast();
          boundary::helpers::QuadTree* south_east = node->SouthEast();
          boundary::helpers::QuadTree* south_west = node->SouthWest();
          std::cout << "test 4" << std::endl;
          // // Dereference
          // boundary::helpers::QuadTree north_west_dref = *north_west;
          // std::cout << "test 5" << std::endl;
          // Find Cell Centers
          boundary::helpers::Point* center_nw = north_west->GetCellCenter();
          std::cout << "test 6" << std::endl;
          // Get Values
          double x_nw = center_nw->x_val;

          std::cout << "x val type" << x_nw << std::endl;
          std::cout << "is nan working?" << std::isnan(x_nw) << std::endl;
          // // If children exist, push them to stack
          // std::cout << "test4" << std::endl;
          // std::cout << north_west->GetCellCenter()->x_val << std::endl;
          // std::cout << std::isnan(north_west->GetCellCenter()->x_val) << std::endl;
          // if (north_west){
          //   std::cout << "test5" << std::endl;
          //   stack.push(north_west);
          // };
          // if (north_east){
          //   stack.push(north_east);
          // };
          // if (south_east){
          //   stack.push(south_east);
          // };
          // if (south_west){
          //   stack.push(south_west);
          // };

          // If no children exist, calculate normals
          if (!north_west && !north_east && !south_east && !south_west){
              boundary::helpers::Point* center = node->GetCellCenter();
              double x_normal = boundary::normals::Normal::NormalDerivative(q, 1, *center,
                                                                  &input);
              double y_normal = boundary::normals::Normal::NormalDerivative(q, 2, *center,
                                                                  &input);
              std::vector<double> normal_vector{x_normal, y_normal};
              node->AssignNormals(normal_vector, q);
          };
        };
      };
    };
  };
  return 0;
}
