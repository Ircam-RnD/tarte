#include <duffing.h>

#include <gtest/gtest.h>

TEST(DuffingTest, Instanciates)
{
    tarte::Duffing<float> model(44100);
    EXPECT_EQ(2, 2);
}