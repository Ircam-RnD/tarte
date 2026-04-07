#include <utility/biquad.h>

#include <gtest/gtest.h>

TEST(BiquadTest, Instanciates)
{
    tarte::Biquad filter(44100.0f);
}

TEST(BiquadTest, Processes)
{
    tarte::Biquad filter(44100.0f);
    double out = filter.Process(1.0f);
}