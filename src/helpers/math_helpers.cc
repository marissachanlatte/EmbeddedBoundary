#include "helpers/math_helpers.h"

#include <stdint.h>

namespace boundary{

namespace helpers{

static const int factorials[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};

int MathHelper::Factorial(int alpha){
  // Since alpha will typically be small, using a lookup table for speed
  if (alpha < 10)
    return factorials[alpha];
  else
    return alpha * Factorial(alpha - 1);
}


// int MathHelpers::MultiIndexFactorial(std::vector<int> alpha){
//
// }
//
//
// int MathHelpers::MultiIndexBinomial(std::vector<int> alpha,
//                                     std::vector<int> beta){
// }

}

}
