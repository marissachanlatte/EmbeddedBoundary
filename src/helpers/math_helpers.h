#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

#include <vector>
#include <Eigen/Dense>

namespace boundary{

namespace helpers{

const int factorials[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};

/**
A function to caculate the factorial of an integer
*/
int Factorial(int alpha);

/**
A function to calculate the multi-index factorial value of a vector of integers
*/
int MultiIndexFactorial(std::vector<int> alpha);

/**
A function to calculate the multi-index binomial coefficient of two vectors of integers
*/
int MultiIndexBinomial(std::vector<int> alpha,
                                   std::vector<int> beta);

/** 
A function that normalizes a point to between 0 and 1 
*/
std::vector<double> NormalizeVector(std::vector<double> coords, std::vector<double> maxes,
                                          std::vector<double> mins);


/** 
 A function to map points in [0, 1]^n to the integers
 */
std::vector<int> IntegerMap(std::vector<double> scaled_coords, int depth);


/**
A function to generate keys from Z Morton Order
*/
double MortonKey(std::vector<double> coords, int depth, std::vector<double> maxes,
                     std::vector<double> mins);


/** 
A function to take the Moore-Penrose Pseudo-Inverse as recommended by Eigen developers
Adapated from https://gist.github.com/pshriwise/67c2ae78e5db3831da38390a8b2a209f
*/
Eigen::MatrixXf PseudoInverse(const Eigen::MatrixXf &a);

/**
A function to caculate the pth moment over a square volume analytically
*/
double PMoment(std::vector<int> p, double cell_size);

} // namespace helpers

} // namespace boundary

#endif // MATH_HELPERS_H