#ifndef EMBEDDED_BOUNDARY_GEOMETRY_BUILD_MESH_H
#define EMBEDDED_BOUNDARY_GEOMETRY_BUILD_MESH_H

#include "inputs/input_base.h"
#include "helpers/geometry_objects.h"

namespace boundary {

namespace geometry {

class MeshTree{
public:
  MeshTree(inputs::InputBase* input);
  ~MeshTree() = default;
  helpers::QuadTree BuildMesh(int level, double right, double left, double top,
                     double bottom);
  helpers::QuadTree mesh;
  int NumNodes();
private:
  int depth_;
  int q_order_;
  int num_nodes_;

};

} // namespace geometry

} // namespace boundary

#endif
