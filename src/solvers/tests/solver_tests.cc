#include "solvers/laplacian.h"
#include "solvers/laplace_operator_ho.h"
#include "inputs/geometries/circle/circle_test.h"
#include "helpers/math_helpers.h"
#include "geometry/boundary.h"

#include "gtest/gtest.h"
#include <iostream>

TEST(HOOperatorTest, Neighborhood){
    // Input
    boundary::inputs::CircleTestGeometry input;
    // Make Geometry
    boundary::geometry::Boundary boundary = boundary::geometry::Boundary(&input);
    boundary::solvers::LaplaceOperatorHO laplace_operator = boundary::solvers::LaplaceOperatorHO(boundary);
    std::vector<double> center_coords{-0.625, -0.625};
    std::vector<double> neighbor_list = laplace_operator.Neighborhood(center_coords, .25);
    std::vector<double> first_neighbor{-0.375, -0.625};
    double first_key = boundary::helpers::MortonKey(first_neighbor, 3, boundary.Maxes(), boundary.Mins());
    EXPECT_TRUE(std::count(neighbor_list.begin(), neighbor_list.end(), first_key));
    std::vector<double> second_neighbor{-0.625, 0.125};
    double second_key = boundary::helpers::MortonKey(second_neighbor, 3, boundary.Maxes(), boundary.Mins());
    EXPECT_TRUE(std::count(neighbor_list.begin(), neighbor_list.end(), second_key));
    std::vector<double> third_neighbor{-0.875, -0.875};
    double third_key = boundary::helpers::MortonKey(third_neighbor, 3, boundary.Maxes(), boundary.Mins());
    EXPECT_FALSE(std::count(neighbor_list.begin(), neighbor_list.end(), third_key));


}


