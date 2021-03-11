#include "helpers/math_helpers.h"
#include "helpers/geometry_objects.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

namespace boundary{

namespace helpers{


/**
A function to caculate the factorial of an integer
*/
int Factorial(int alpha){
  // Since alpha will typically be small, using a lookup table for speed
  if (alpha < 10)
    return factorials[alpha];
  else
    return alpha * Factorial(alpha - 1);
}

/**
A function to calculate the multi-index factorial value of a vector of integers
*/
int MultiIndexFactorial(std::vector<int> alpha){
  int total = 1;
  for(std::vector<int>::iterator it = alpha.begin(); it != alpha.end(); ++it){
    total *= Factorial(*it);
  }
  return total;
}

/**
A function to calculate the multi-index binomial coefficient of two vectors of integers
*/
int MultiIndexBinomial(std::vector<int> alpha,
                                   std::vector<int> beta){
  std::vector<int> difference;
  std::transform(alpha.begin(), alpha.end(), beta.begin(),
                 std::back_inserter(difference), std::minus<int>());
  return (MultiIndexFactorial(alpha)/
          (MultiIndexFactorial(beta)*MultiIndexFactorial(difference)));
}


/** 
A function that normalizes a point to between 0 and 1 
*/
std::vector<double> NormalizeVector(std::vector<double> coords, std::vector<double> maxes,
                                          std::vector<double> mins){
  std::vector<double> normalized_vector;
  int dim = coords.size();
  for (int d = 0; d < dim; d++){
    normalized_vector.push_back((coords[d] - mins[d])/(maxes[d] - mins[d]));
  }
  return normalized_vector;
}


/** 
 A function to map points in [0, 1]^n to the integers
 */
std::vector<int> IntegerMap(std::vector<double> scaled_coords, int depth){
  // Vector of indices
  std::vector<int> indices;
  // Cell size
  double cell_size = std::pow(2, -depth);
  for (auto& it : scaled_coords){
    indices.push_back(int(std::floor(it/cell_size)));
  }
  return indices;
}


/**
A function to generate keys from Z Morton Order
*/
int MortonKey(std::vector<double> coords, int depth, std::vector<double> maxes,
                     std::vector<double> mins){
  int dim = coords.size();
  // scale coordinates to 0 - 1
  std::vector<double> scaled = NormalizeVector(coords, maxes, mins);
  // Use mesh spacing to map to integers
  std::vector<int> mesh_mapped = IntegerMap(scaled, depth);
  // Initialize key
  int key = 0;
  // Keep track of number of digits of key
  int digit = 0;
  // Loop through depths
  for (int k = 0; k < depth; k ++){
    // Loop through dimensions
    for (int d = (dim - 1); d >= 0; d--){
      // Get kth bit
      int bit = ((mesh_mapped[d] & (1 << k)) >> k);
      // add to key in digit place
      key += bit * std::pow(10, digit);
      // increment digit
      digit += 1;
    }
  }
  // Add leading 1
  key += std::pow(10, digit);
  return key;
}

}

}
