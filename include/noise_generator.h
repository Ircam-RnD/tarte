#pragma once

#include <array>
#include <bit>
#include <cstdint>
#include <random>

namespace tarte {

enum NoiseColor { White, Pink };

class NoiseGenerator {
private:
    static constexpr uint32_t kNumStreams = 16;

    NoiseColor color_;

    uint32_t counter_ = 0;

    std::array<float, kNumStreams> state_{};
    float running_sum_ = 0.0f;

    std::mt19937 rng_{std::random_device{}()};
    std::uniform_real_distribution<float> dist_{-1.0f, 1.0f};

    float GenerateWhite() { return dist_(rng_); }

public:
    explicit NoiseGenerator(NoiseColor color = NoiseColor::Pink): color_(color)
    {
        if (color_ == NoiseColor::Pink) {
            for (auto& s: state_) {
                s = GenerateWhite();
                running_sum_ += s;
            }
        }
    }

    float Process()
    {
        if (color_ == NoiseColor::White)
            return GenerateWhite();

        ++counter_;

        uint32_t changed = std::countr_zero(counter_);

        if (changed < kNumStreams) {
            float old_value = state_[changed];
            float new_value = GenerateWhite();

            state_[changed] = new_value;
            running_sum_ += new_value - old_value;
        }

        return running_sum_ / static_cast<float>(kNumStreams);
    }
};

} // namespace tarte