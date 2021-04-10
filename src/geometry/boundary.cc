#include "geometry/boundary.h"
#include "helpers/geometry_objects.h"
#include "helpers/math_helpers.h"
#include "normals/normals.h"

#include <array>
#include "math.h"
#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>

namespace boundary {

namespace geometry {

Boundary::Boundary(boundary::inputs::GeometryInputBase* input){
  input_ = input;
  // Set global variables
  Q_ = input->QOrder();
  mins_.push_back(input->XMin());
  mins_.push_back(input->YMin());
  maxes_.push_back(input->XMax());
  maxes_.push_back(input->YMax());
  initial_cell_size_ = maxes_[0] - mins_[0];
  max_depth_ = input->MaxDepth();
  max_solver_depth_ = input->MaxSolverDepth();

  SetupMesh_(initial_cell_size_, mins_[1], maxes_[1], mins_[0], maxes_[0]);
  RecursiveCalculateMoments_(1, initial_cell_size_);
  PropagateUp_();
};

void Boundary::RecursiveCalculateMoments_(double key, double cell_size){
  // Check if child exists. If it does, check for their children... 
  bool has_child = false;
  for (int it = 0; it < 4; it++){
    // Iterating through all four possible child cells
    // TODO: Fix for arbitrary dim
    double new_key = 100*key + (it&1) + (it > 1)*10;
    if (boundary_cells_.count(new_key)){
      has_child = true;
      RecursiveCalculateMoments_(new_key, cell_size/2);
    }
  }
  // If no children, calculate moments
  if (!has_child && boundary_cells_.count(key)) {
    CalculateMoments_(key, cell_size);
  }
}


void Boundary::PropagateUp_(){
  for (int depth = max_depth_; depth >= 0; depth--){
    int n = 2*depth;
    double child_cell_volume = std::pow(initial_cell_size_/std::pow(2, depth + 1), 2);
    // Iterate through all cells at this level
    for (int it=0; it < (1 << n); it++){
      double key = std::pow(10, n);
      for (int jt=0; jt < n; jt++){
        key += ((it>>jt)&1)*std::pow(10, jt);
      }
      // If a boundary cell, sum up all children
      if (boundary_cells_.count(key)){
        for (int it = 0; it < 4; it++){ // Loop through children 
          double child_key = 100*key + (it&1) + (it > 1)*10;
          for (int q_mag = Q_; q_mag >= 0; q_mag--){ // Loop through q order
            for (int i = 0; i < q_mag + 1; i++){  // Volume moments
              // Check if child is on boundary
              if (boundary_cells_.count(child_key)){
                // Propagate volume moments
                if (i < q_mag){boundary_cells_[key].volume_moments[i][q_mag - 1 - i] += boundary_cells_[child_key].volume_moments[i][q_mag - 1 - i];}
                // Propagate boundary moments
                boundary_cells_[key].boundary_moments[i][q_mag - i] += boundary_cells_[child_key].boundary_moments[i][q_mag - i];
              }
              // If child cell is interior add whole volume
              else {
                if (cell_map_[child_key] == 1 && (i < q_mag)){
                  boundary_cells_[key].volume_moments[i][q_mag - 1 - i] += child_cell_volume;
                }
              }
            } 
          }
        }
        // If child below solver max level, delete
        if (depth >= max_solver_depth_){
          for (int it = 0; it < 4; it++){ // Loop through children 
            double child_key = 100*key + (it&1) + (it > 1)*10;
            // Delete child
            boundary_cells_.erase(child_key);
          }
        }
      }
    }
  }  
};


void Boundary::SetupMesh_(double cell_size, double y_min, double y_max, double x_min, 
                         double x_max){
  // Caculate depth
  int depth = int(std::log2((initial_cell_size_/cell_size)));
  // Iterate through all cells
  double y_curr = y_min;
  double y_next = y_min + cell_size;
  while (y_next <= y_max){
    double x_curr = x_min;
    double x_next = x_curr + cell_size;
    while (x_next <= x_max){
      // corners of cell
      std::array<std::array<double, 2>, 4> corners;
      corners[0] = {x_curr, y_curr}; // lower left
      corners[1] = {x_curr, y_next}; // upper left
      corners[2] = {x_next, y_next}; // upper right
      corners[3] = {x_next, y_curr}; // lower right

      // cell center
      std::vector<double> center {x_curr + cell_size/2, y_curr + cell_size/2};

      // calculate Morton key
      double key = helpers::MortonKey(center, depth, maxes_, mins_);

      std::vector<int> inside{input_->Inside(corners[0]),
                              input_->Inside(corners[1]),
                              input_->Inside(corners[2]),
                              input_->Inside(corners[3])};

      // Refine
      // TODO: come up with a real criterion, this one just refines all cells to max
      if (depth < max_depth_){
        SetupMesh_(cell_size/2, y_curr, center[1], x_curr, center[0]); // lower left
        SetupMesh_(cell_size/2, center[1], y_next, x_curr, center[0]); // upper left
        SetupMesh_(cell_size/2, center[1], y_next, center[0], x_next); // upper right
        SetupMesh_(cell_size/2, y_curr, center[1], center[0], x_next); // lower right
      }

      // Determine if cell is boundary cell
      bool is_boundary = false;
      if (depth == max_depth_){
        // Check if four corners are inside or outside boundary
        is_boundary = IsBoundaryCell(corners[0], corners[1],
                                     corners[2], corners[3], input_);
      }

      else {
        // Check if any of the child cells are boundary cells
        for (int it = 0; it < 4; it++){ // Loop through children 
          double child_key = 100*key + (it&1) + (it > 1)*10;
          if (cell_map_[child_key] == 2){
            is_boundary = true;
          }
        }
      }

      if (is_boundary){
        // add to cell map 
        cell_map_.insert(std::pair<double, int>(key, 2));
        // make struct with tag and id
        geo_info cell;
        cell.irregular = true;
        helpers::Point cell_center = helpers::Point(center);
        // resize vectors
        cell.normal_derivatives.resize(Q_+1);
        cell.volume_moments.resize(Q_+1);
        cell.boundary_moments.resize(Q_+1);
        for (int i=0; i < Q_+1; ++i){
          cell.normal_derivatives[i].resize(Q_+1);
          cell.volume_moments[i].resize(Q_+1);
          cell.boundary_moments[i].resize(Q_+1);
          for (int j = 0; j < Q_ + 1; ++j){
            cell.normal_derivatives[i][j].resize(2);
          }
        }
        // Compute derivatives of the normal and store in cell struct
        for (int q_mag=Q_; q_mag >= 0 ; q_mag--){
          for (int q1 = 0; q1 < q_mag + 1; q1++){
            std::vector<int> q = {q1, q_mag - q1};
            cell.normal_derivatives[q1][q_mag - q1][0] = normals::Normal::NormalDerivative(q,
                                                            1, cell_center, input_);
            cell.normal_derivatives[q1][q_mag - q1][1] = normals::Normal::NormalDerivative(q,
                                                            2, cell_center, input_);
          }
        }

        // Compute 1d volume fractions and store in cell
        for (int corner=0; corner < 4; corner++){
          // iterate four edges to determine which ones intersect boundary
          if ((inside[corner] == inside[(corner+1)%4]) || inside[corner] == 2){ // no change from corner to corner, TODO: come up with more robust way to do this
            cell.vol_frac_1d[corner] = cell_size*inside[(corner+1)%4];
          }
          else if (inside[(corner+1)%4] == 2){ // next corner is on the boundary
            cell.vol_frac_1d[corner] = cell_size*inside[corner];
          }
          else { // change from corner to corner, indicating boundary cuts through
            // check if edge is horizontal or vertical
            if ((corners[corner][0] - corners[(corner+1)%4][0]) == 0){ // horizontal
              // find the intersection of x=corners[i, 0]
              std::vector<double> y_values = input_->BoundaryFunction(corners[corner][0]);
              double y = 0;
              if (y_values.size() == 1){
                y = input_->BoundaryFunction(corners[corner][0])[0];
              }
              else if (y_values.size() == 2){
                // choose the value that is in the cell
                // Robustness TODO: fix this in the case that both values are in cell
                y = WhichValue(y_values, corners[corner][1], corners[(corner+1)%4][1]);
              }

              double edge_inside = abs(y - corners[corner][1]);
              if (inside[corner]){
                cell.vol_frac_1d[corner] = edge_inside;
              }
              else{
                cell.vol_frac_1d[corner] = cell_size - edge_inside;
              }
            }
            else { // vertical
              // find intersection of y=corners[i, 1]
              std::vector<double> x_values = input_->BoundaryInverse(corners[corner][1]);
              double x = 0;
              if (x_values.size() == 1){
                x = input_->BoundaryFunction(corners[corner][0])[0];
              }
              else if (x_values.size() == 2){
                // choose the value that is in the cell
                // Robustness TODO: fix this in the case that both values are in cell
                x = WhichValue(x_values, corners[corner][0], corners[(corner+1)%4][0]);
              }
              double edge_inside = abs(x - corners[corner][0]);
              if (inside[corner]){
                cell.vol_frac_1d[corner] = edge_inside;
              }
              else {
                cell.vol_frac_1d[corner] = cell_size - edge_inside;
              }
            }
          }
        }
        
        cell.cell_center = center;
        cell.cell_size = cell_size;
        // add cell to map
        boundary_cells_.insert(std::pair<double, geo_info>(key, cell));
      }

      // Else add to cell map
      else {cell_map_.insert(std::pair<double, int>(key, *std::min_element(inside.begin(), inside.end())));}

      // Add id to cell center map
      id_to_center_.insert(std::pair<double, std::vector<double>>(key, center));

      // Go to next cell
      x_next += cell_size;
      x_curr += cell_size;
    }
    // Go to next row
    y_next += cell_size;
    y_curr += cell_size;
  }
}

double Boundary::WhichValue(std::vector<double> values, double first_bound, double second_bound){
  if (first_bound < second_bound){
    for (std::vector<double>::iterator it = values.begin(); it != values.end(); it++){
      if (*it >= first_bound && *it <= second_bound){
        return *it;
      }
    }
  }

  else if (first_bound > second_bound){
    for (std::vector<double>::iterator it = values.begin(); it != values.end(); it++){
      if (*it <= first_bound && *it >= second_bound){
        return *it;
      }
    }
  }

  else {
    throw std::invalid_argument("No change between bounds.");
  }
  return 1;
};


bool Boundary::IsBoundaryCell(std::array<double, 2> lower_left,
                              std::array<double, 2> upper_left,
                              std::array<double, 2> upper_right,
                              std::array<double, 2> lower_right,
                              boundary::inputs::GeometryInputBase* input){
  // put into array with inside values
  std::vector<int> inside{input->Inside(lower_left),
                          input->Inside(upper_left),
                          input->Inside(upper_right),
                          input->Inside(lower_right)};
  // try changing all 2s to 0s and see if they are all the same
  int it_prev = inside[0];
      if (it_prev == 2){
      it_prev = 0;
    }
  bool same;
  for (int it = 1; it < 4; it++){
    int current = inside[it];
    if (current == 2){
      current = 0;
    }
    if (it_prev == current){
      if (it == 3){
        // they're all the same, not a boundary cell
        same = true;
      }
      continue;
    }
    else{
      same = false;
      break;
    }
  }

  // if they're not try changing all the 2s to 1s and see if they're all the same
  if (!same){
    it_prev = inside[0];
    if (it_prev == 2){
      it_prev = 1;
    }
    for (int it = 1; it < 4; it++){
      int current = inside[it];
      if (current == 2){
        current = 1;
      }
      if (it_prev == current){
        if (it == 3){
          same = true;
        }
        continue;
      }
      else{
        same = false;
        break;
      }
    }
  }

  return !same;
};


std::map<double, geo_info> Boundary::BoundaryCells(){
  return boundary_cells_;
};


double Boundary::DIntegral_(double beginning, double end, std::array<int, 2> q,
                            int index, double fixed_value){
  double next_plus = q[(index + 1) % 2] + 1;
  return std::pow(fixed_value, q[index])*(1/(next_plus))*(std::pow(end, next_plus) - std::pow(beginning, next_plus));
};


double Boundary::CalcD_(double bd_length,
              double fixed_value,
              std::vector<double> cell_center,
              std::array<int, 2> q,
              int d,
              std::array<int, 2> which_d,
              double cell_size){
  int d_op = (d + 1) % 2;
  double d_pm;
  if (bd_length == cell_size){
    d_pm = DIntegral_(cell_center[d_op] - cell_size/2,
                            cell_center[d_op] + cell_size/2,
                            q, d, fixed_value);
  }
  else if (bd_length > 0){
    // TODO: fix this to be more general
    std::array<double, 2> corner = {cell_center[0] + cell_size/2*which_d[0],
                                    cell_center[1] + cell_size/2*which_d[1]};
    std::vector<double> int_values;
    if (d == 0){
      int_values = input_->BoundaryFunction(fixed_value);
    }
    else if (d == 1){
      int_values = input_->BoundaryInverse(fixed_value);
    }
    double intersection;
    if (int_values.size() == 1){
      intersection = int_values[0];
    }
    else if (int_values.size() == 2){
      // choose the value that is in the cell
      // Robustness TODO: fix this in the case that both values are in cell
      intersection = WhichValue(int_values, cell_center[d_op] - cell_size/2, 
                                cell_center[d_op] + cell_size/2);
    }
    if (input_->Inside(corner)){
      d_pm = DIntegral_(intersection, corner[d_op], q, d, fixed_value);
    }
    else {
      d_pm = DIntegral_(corner[d_op] - cell_size, intersection, q, d, fixed_value);
    }
  }
  else {
    d_pm = 0;
  }
  return d_pm;
};


void Boundary::CalculateMoments_(double key, double cell_size){
  std::vector<double> cell_center = id_to_center_[key];
  for (int q_mag = Q_; q_mag >= 0; q_mag--){
    Eigen::VectorXf rho = Eigen::VectorXf::Zero(2*(q_mag + 1));
    Eigen::MatrixXf lhs = Eigen::MatrixXf::Zero(2*(q_mag + 1), 2*q_mag + 1);

    int it = -1;
    for (int q1 = 0; q1 < (q_mag + 1); q1++){
      std::array<int, 2> q = {q1, q_mag - q1};
      for (int d = 0; d < 2; d++){
        it += 1;
        std::array<int, 2> e = {0, 0};
        e[d] = 1;

        std::array<int, 2> plus_vec = {1, 1};
        std::array<int, 2> first_minus_vec = {-1, 1};
        std::array<int, 2> second_minus_vec = {1, -1};

        double ptv_bd_length;
        double d_plus;
        double ntve_bd_length;
        double d_minus;
        // set d_plus and d_minus
        if (d == 0){
          // calculate d_plus
          ptv_bd_length = boundary_cells_[key].vol_frac_1d[2];
          double fixed_ptve_x = cell_center[0] + cell_size/2;
          d_plus = CalcD_(ptv_bd_length, fixed_ptve_x, cell_center, q, d, plus_vec, cell_size);

          // calculate d_minus
          ntve_bd_length = boundary_cells_[key].vol_frac_1d[0];
          double fixed_ntve_x = cell_center[0] - cell_size/2;
          d_minus = CalcD_(ntve_bd_length, fixed_ntve_x, cell_center, q,
                           d, first_minus_vec, cell_size);
        }
        else if (d == 1){
          // calculate d_plus
          ptv_bd_length = boundary_cells_[key].vol_frac_1d[1];
          double fixed_ptve_y = cell_center[1] + cell_size/2;
          d_plus = CalcD_(ptv_bd_length, fixed_ptve_y, cell_center, q,
                                 d, plus_vec, cell_size);

          // calculate d_minus
          ntve_bd_length = boundary_cells_[key].vol_frac_1d[3];
          double fixed_ntve_y = cell_center[1] - cell_size/2;
          d_minus = CalcD_(ntve_bd_length, fixed_ntve_y, cell_center, q,
                           d, second_minus_vec, cell_size);
        }
        int s_sum = 0;
        int S = Q_ - q_mag;
        for (int s1 = 0; s1 < S+1; s1++){
          for (int s2 = 0; s2 < S+1; s2++){
            if (s1 + s2 >= S || s1 + s2 < 1){
              continue;
            }
            s_sum += (boundary_cells_[key].normal_derivatives[s1][s2][d]
                      *boundary_cells_[key].boundary_moments[q[0] + s1][q[1] + s2]);
          }
        }
        std::array<int, 2> q_minus_e = {q[0] - e[0], q[1] - e[1]};
        if (q_minus_e[0] >= 0 && q_minus_e[1] >= 0){
          lhs(it, q_minus_e[0]) = q[d];
        }
        // M_B Unknowns (stored based on value of q[0
        lhs(it, q_mag + q[0]) = -boundary_cells_[key].normal_derivatives[0][0][d];
        rho(it) = d_plus - d_minus + s_sum;
         
      }
    }

    // Solve
    Eigen::VectorXf v_and_b = lhs.colPivHouseholderQr().solve(rho);

    // Unpack V and B
    for (int i = 0; i < q_mag; i++){
      boundary_cells_[key].volume_moments[i][q_mag - 1 - i] = v_and_b(i);
    }
    for (int i = 0; i < q_mag + 1; i++){
      boundary_cells_[key].boundary_moments[i][q_mag - i] = v_and_b(q_mag + i);
    }
  }
};


// std::vector<double> Boundary::IDToCenter(double id){
//   return id_to_center_[id];
// };


int Boundary::IJToGlobal(int i_index, int j_index, int depth){
  int num_x = std::pow(2, depth);
  if(i_index < 0 || j_index < 0 || i_index >= num_x || j_index >= num_x){
    return -1;
  }
  return num_x*i_index + j_index;
};


std::array<int, 2> Boundary::NeighborCell(int i_index, int j_index, int edge){
  std::map<int, std::array<int, 2>> neighbor_map = {{0, {i_index, j_index - 1}},
                                                    {1, {i_index + 1, j_index}},
                                                    {2, {i_index, j_index + 1}},
                                                    {3, {i_index - 1, j_index}}};
  return neighbor_map[edge];
}


std::vector<double> Boundary::IJToCenter(int i_index, int j_index, int depth){
  int num_x = std::pow(2, depth);
  double cell_size = (maxes_[0] - mins_[0])/num_x;
  std::vector<double> cell_center{mins_[0] + i_index*cell_size + cell_size/2, 
                                  mins_[1] + j_index*cell_size + cell_size/2};
  return cell_center;
}


int Boundary::Sgn_(double v){
    // branchless sign function.
    // https://helloacm.com/how-to-implement-the-sgn-function-in-c/
    return (v > 0) - (v < 0);
}


std::array<int, 2> Boundary::ProjectedNormal_(int side_index, double nx, double ny){
    // (side_index & 1) iff care about y direction.
    // thus we can create a projected normal using comparison + bit mask.
    std::array<int, 2> normal = {Sgn_((    (side_index & 1)) * nx),
                                 Sgn_((1 ^ (side_index & 1)) * ny)};
    return normal;
}


int Boundary::Parity_(int side_index){
    // return parity of 2 bit number
    return (side_index ^ (side_index >> 1)) & 1;
}


std::array<std::array<int, 2>, 2> Boundary::InterpolationPair(int i, int j, double nx, double ny, int side_index){
  // want to move along projected direction of negative normal.
  std::array<int, 2> steps = ProjectedNormal_(side_index, -nx, -ny);
  
  // take step along projected normal.
  int i2 = i + steps[1]; 
  int j2 = j + steps[0];

  int direction = Parity_(side_index) - (not Parity_(side_index));
  steps[0] = direction * ((side_index & 1) ^ 1);
  steps[1] = direction *  (side_index & 1);

  int i3 = i2 + steps[1];
  int j3 = j2 + steps[0];

  std::array<int, 2> first_pair = {i2, j2};
  std::array<int, 2> second_pair = {i3, j3};
  std::array<std::array<int, 2>, 2> pair = {first_pair, second_pair};
  return pair;
}


std::map<double, int> Boundary::CellMap(){
  return cell_map_;
};


double Boundary::InitialCellSize(){
  return initial_cell_size_;
};


double Boundary::XMin(){
  return mins_[0];
};


double Boundary::XMax(){
  return maxes_[0];
};


double Boundary::YMin(){
  return mins_[1];
};


double Boundary::YMax(){
  return maxes_[1];
};


std::vector<double> Boundary::Mins(){
  return mins_;
};


std::vector<double> Boundary::Maxes(){
  return maxes_;
};


int Boundary::MaxSolverDepth(){
  return max_solver_depth_;
}

} // namespace geometry

} // namespace boundary
