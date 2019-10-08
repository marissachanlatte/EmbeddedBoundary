 #include "inputs/circle/circle.h"
 #include "helpers/geometry_objects.h"

 #include "gtest/gtest.h"
 #include <vector>

 class CircleInputTest : public ::testing::Test{
   protected:
     CircleInputTest();
 };

 TEST(CircleInputTest, Boundary){
   boundary::inputs::CircleGeometry circle;
   EXPECT_EQ(0, circle.BoundaryFunction(1)[0]);
   EXPECT_EQ(std::sqrt(.75), circle.BoundaryFunction(.5)[0]);
   EXPECT_EQ(-std::sqrt(.75), circle.BoundaryFunction(.5)[1]);
 }

 // TEST(CircleInputTest, Derivative){
 //   boundary::inputs::CircleGeometry circle;
 //
 //   boundary::helpers::Point point1 = boundary::helpers::Point(1, 1);
 //   boundary::helpers::Point point2 = boundary::helpers::Point(3, 3);
 //   boundary::helpers::Point point3 = boundary::helpers::Point(4, 4);
 //   boundary::helpers::Point point4 = boundary::helpers::Point(6, 6);
 //
 //   std::vector<int> order1{1, 0};
 //   std::vector<int> order2{0, 1};
 //   std::vector<int> order3{0, 0};
 //   std::vector<int> order4{4, 5};
 //
 //   EXPECT_EQ(-1, circle.BoundaryDerivatives(point1, order1));
 //   EXPECT_EQ(1, circle.BoundaryDerivatives(point2, order2));
 //   EXPECT_EQ(0, circle.BoundaryDerivatives(point3, order3));
 //   EXPECT_EQ(0, circle.BoundaryDerivatives(point4, order4));
 // }
 //
 TEST(CircleInputTest, Inverse){
   boundary::inputs::CircleGeometry circle;
   EXPECT_EQ(0, circle.BoundaryInverse(1)[0]);
   EXPECT_EQ(std::sqrt(.75), circle.BoundaryInverse(.5)[0]);
   EXPECT_EQ(-std::sqrt(.75), circle.BoundaryInverse(.5)[1]);
 }
