#include "../Normals.h"

#include "gtest/gtest.h"

class NormalTest : public ::testing::Test{
  protected:
    NormalTest();
};

TEST(NormalTest, dummy){
  Normal test_normal;
  EXPECT_EQ(0, test_normal.makeNormal());
}
