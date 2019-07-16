#ifndef EMBEDDED_BOUNDARY_HELPERS_MATH_HELPERS_H
#define EMBEDDED_BOUNDARY_HELPERS_MATH_HELPERS_H

#include <vector>

namespace boundary {

namespace helpers {

  class MathHelper{
    public:
      ~MathHelper() = default;
      int Factorial(int alpha);
      int MultiIndexFactorial(std::vector<int> alpha);
      int MultiIndexBinomial(std::vector<int> alpha,
                                std::vector<int> beta);
  };

} // namespace helpers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_HELPERS_MATH_HELPERS_H
