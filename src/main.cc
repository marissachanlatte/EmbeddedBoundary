#include "inputs/line/line.h"
#include "inputs/circle/circle.h"
#include "helpers/geometry_objects.h"
#include "normals/normals.h"
#include "geometry/boundary.h"

#include <stack>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

int main(){
  boundary::inputs::CircleGeometry input;
  boundary::geometry::Boundary line_boundary = boundary::geometry::Boundary(&input);
  std::map<std::array<double, 2>, boundary::geometry::geo_info> boundary_cells = line_boundary.BoundaryCells();
  std::ofstream output;
  output.open ("../outputs/output_boundary.txt");
  for (std::map<std::array<double, 2>, boundary::geometry::geo_info>::iterator it=boundary_cells.begin();
       it != boundary_cells.end(); it++){
         output << it->first[0] << " " << it->first[1] << " " << it->second.volume_moments[0][0] << std::endl;
       }
  output.close();
  return 0;
}
