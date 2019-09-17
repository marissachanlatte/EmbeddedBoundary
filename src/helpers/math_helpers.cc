#include "helpers/geometry_objects.h"
#include <vector>
#include <algorithm>
#include <functional>

namespace boundary{

namespace helpers{

static const int factorials[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};

/**
A function to caculate the factorial of an integer
*/
static int Factorial(int alpha){
  // Since alpha will typically be small, using a lookup table for speed
  if (alpha < 10)
    return factorials[alpha];
  else
    return alpha * Factorial(alpha - 1);
}

/**
A function to calculate the multi-index factorial value of a vector of integers
*/
static int MultiIndexFactorial(std::vector<int> alpha){
  int total = 1;
  for(std::vector<int>::iterator it = alpha.begin(); it != alpha.end(); ++it){
    total *= Factorial(*it);
  }
  return total;
}

/**
A function to calculate the multi-index binomial coefficient of two vectors of integers
*/
static int MultiIndexBinomial(std::vector<int> alpha,
                                   std::vector<int> beta){
  std::vector<int> difference;
  std::transform(alpha.begin(), alpha.end(), beta.begin(),
                 std::back_inserter(difference), std::minus<int>());
  return (MultiIndexFactorial(alpha)/
          (MultiIndexFactorial(beta)*MultiIndexFactorial(difference)));
}

}

}
