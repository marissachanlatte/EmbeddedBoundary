 #include "inputs/geometries/circle/circle_test.h"
 #include "helpers/geometry_objects.h"

 #include "gtest/gtest.h"
 #include <vector>

 class CircleInputTest : public ::testing::Test{
   protected:
     CircleInputTest();
 };

 TEST(CircleInputTest, Boundary){
   boundary::inputs::CircleTestGeometry circle;
   EXPECT_EQ(0, circle.BoundaryFunction(1)[0]);
   EXPECT_EQ(std::sqrt(.75), circle.BoundaryFunction(.5)[0]);
   EXPECT_EQ(-std::sqrt(.75), circle.BoundaryFunction(.5)[1]);
 }

 TEST(CircleInputTest, Inverse){
   boundary::inputs::CircleTestGeometry circle;
   EXPECT_EQ(0, circle.BoundaryInverse(1)[0]);
   EXPECT_EQ(std::sqrt(.75), circle.BoundaryInverse(.5)[0]);
   EXPECT_EQ(-std::sqrt(.75), circle.BoundaryInverse(.5)[1]);
 }
