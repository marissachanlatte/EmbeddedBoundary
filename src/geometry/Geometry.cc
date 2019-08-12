#include "inputs/line/line.h"
#include "helpers/geometry_objects.h"
#include "geometry/build_mesh.h"
#include "normals/normals.h"

#include <stack>
#include <vector>

namespace boundary {

namespace geometry {

int main(){
  // Setup Input
  inputs::LineGeometry input;
  int q_order = input.QOrder();

  // Build Mesh
  MeshTree mesh_tree = MeshTree(&input);

  // Build Normals
  for (int q_mag = q_order; q_mag > -1; q_mag--){
    for (int q1 = 0; q1 < (q_mag + 1); q1 ++){
      std::vector<int> q{q1, q_mag - q1};
      // Initialize Stack
      std::stack<helpers::QuadTree*> stack;
      int num_nodes = mesh_tree.NumNodes();
      std::vector<bool> discovered(num_nodes);
      helpers::QuadTree* mesh = &mesh_tree.mesh;
      stack.push(mesh);

      while (!stack.empty()){
        helpers::QuadTree* node = stack.top();
        stack.pop();
        if (!discovered[node->GetID()]){
          // Mark as Discovered
          discovered[node->GetID()] = true;

          // Identify Children
          helpers::QuadTree* north_west = node->NorthWest();
          helpers::QuadTree* north_east = node->NorthEast();
          helpers::QuadTree* south_east = node->SouthEast();
          helpers::QuadTree* south_west = node->SouthWest();

          // If children exist, push them to stack
          if (north_west){
            stack.push(north_west);
          };
          if (north_east){
            stack.push(north_east);
          };
          if (south_east){
            stack.push(south_east);
          };
          if (south_west){
            stack.push(south_west);
          };

          // If no children exist, calculate normals
          if (!north_west && !north_east && !south_east && !south_west){
              helpers::Point* center = node->GetCellCenter();
              double x_normal = normals::Normal::NormalDerivative(q, 1, *center,
                                                                  &input);
              double y_normal = normals::Normal::NormalDerivative(q, 2, *center,
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

} // namespace geometry

} // namespace boundary
