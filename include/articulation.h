#pragma once

#include "vowels.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace tarte {

class Articulation {
private:
    std::array<float, 6> areas_, positions_;
    std::array<float, 4> etas_; // Etas are the "warping exponent"

public:
    Articulation(vowels::Vowel const& v = vowels::a) { SetFromVowel(v); };

    void SetFromVowel(vowels::Vowel const& v)
    {
        positions_ = v.positions;
        areas_ = v.areas;
        etas_ = v.etas;
    };

    void Interpolate2Vowels(vowels::Vowel const& v1, vowels::Vowel const& v2, float alpha)
    {
        // Vowel = (1 - alpha) * v1 + alpha * v2
        for (size_t i = 0; i < 6; ++i) {
            positions_[i] = v1.positions[i] + alpha * (v2.positions[i] - v1.positions[i]);
            areas_[i] = v1.areas[i] + alpha * (v2.areas[i] - v1.areas[i]);
            if (i < 4) {
                etas_[i] = v1.etas[i] + alpha * (v2.etas[i] - v1.etas[i]); // Only 4 exponents
            }
        }
    };

    void InterpolateNVowels(vowels::Vowel const* vs, float const* alphas, size_t const n)
    {
        // Weighted interpolation between N vowels
        float sum_alphas = 0;
        for (size_t i = 0; i < n; ++i) {
            sum_alphas += alphas[i];
        }
        for (size_t i = 0; 6; ++i) {
            positions_[i] = 0;
            areas_[i] = 0;
            if (i < 4) {
                etas_[i] = 0;
            }
            for (size_t j = 0; j < n; ++j) {
                positions_[i] += alphas[j] / sum_alphas * vs[j].positions[i];
                areas_[i] += alphas[j] / sum_alphas * vs[j].areas[i];
                if (i < 4) {
                    etas_[i] += alphas[j] / sum_alphas * vs[j].etas[i];
                }
            }
        }
    };

    // Find N closest vowels to target F0 and F1 frequencies and compute interpolation coefficients
    // Fills vectors closest_vowels and alphas thzat can be then passed to InterpolateNVowels
    void FindClosestVowels(float targetF0,
                           float targetF1,
                           size_t N,
                           std::vector<vowels::Vowel>& closest_vowels,
                           std::vector<float>& alphas)
    {
        if (N == 0)
            return;

        // Create vector of pairs: (distance, vowel pointer)
        std::vector<std::pair<float, vowels::Vowel const*>> vowel_distances;

        // Calculate distance for each vowel
        for (const auto& vowel: vowels::vowels) {
            if (vowel.formant_frequencies.size() >= 2) {
                float f0 = vowel.formant_frequencies[0];
                float f1 = vowel.formant_frequencies[1];

                // Euclidean distance in F0-F1 space
                float distance = std::sqrt((targetF0 - f0) * (targetF0 - f0) + (targetF1 - f1) * (targetF1 - f1));

                vowel_distances.push_back({distance, &vowel});
            }
        }

        // Sort by distance (closest first)
        std::sort(vowel_distances.begin(), vowel_distances.end());

        // Take the N closest vowels
        size_t numVowels = std::min(N, vowel_distances.size());
        closest_vowels.clear();
        closest_vowels.reserve(numVowels);

        // Calculate weights based on inverse distance
        std::vector<float> weights;
        weights.reserve(numVowels);

        for (size_t i = 0; i < numVowels; ++i) {
            closest_vowels.push_back(*vowel_distances[i].second);

            float distance = vowel_distances[i].first;
            float weight = 1.0f / (distance + 1e-6f);
            weights.push_back(weight);
        }

        // Normalize weights to sum to 1
        float totalWeight = 0.0f;
        for (float w: weights) {
            totalWeight += w;
        }

        alphas.clear();
        alphas.reserve(numVowels);
        for (float w: weights) {
            alphas.push_back(w / totalWeight);
        }
    }

    // Convenience function that directly sets the articulation from F0/F1
    // Works by interpolating between N nearest vowels in the F0/F1 space
    void SetFromFormants(float targetF0, float targetF1, size_t N = 3)
    {
        std::vector<vowels::Vowel> closets_vowels;
        std::vector<float> alphas;

        FindClosestVowels(targetF0, targetF1, N, closets_vowels, alphas);

        if (!closets_vowels.empty()) {
            InterpolateNVowels(closets_vowels.data(), alphas.data(), closets_vowels.size());
        }
    }

    template<typename T>
    void getAreas(T const* evaluation_positions, T* out, std::size_t const size)
    {
        double xscale = positions_[5] / evaluation_positions[size - 1];
        for (size_t i = 0; i < size; ++i) {
            if (evaluation_positions[i] * xscale < positions_[0]) {
                out[i] = areas_[0];
            } else if (evaluation_positions[i] * xscale < positions_[1]) {
                out[i] = (areas_[1] + areas_[0]) / 2 +
                         (areas_[1] - areas_[0]) / 2 *
                             std::cos(M_PI * std::pow((positions_[1] - evaluation_positions[i] * xscale) /
                                                          (positions_[1] - positions_[0]),
                                                      etas_[0]));
            } else if (evaluation_positions[i] * xscale < positions_[2]) {
                out[i] = (areas_[2] + areas_[1]) / 2 +
                         (areas_[2] - areas_[1]) / 2 *
                             std::cos(M_PI * std::pow((positions_[2] - evaluation_positions[i] * xscale) /
                                                          (positions_[2] - positions_[1]),
                                                      etas_[1]));
            } else if (evaluation_positions[i] * xscale < positions_[3]) {
                out[i] = (areas_[3] + areas_[2]) / 2 +
                         (areas_[3] - areas_[2]) / 2 *
                             std::cos(M_PI * std::pow((positions_[3] - evaluation_positions[i] * xscale) /
                                                          (positions_[3] - positions_[2]),
                                                      etas_[2]));
            } else if (evaluation_positions[i] * xscale < positions_[4]) {
                out[i] = (areas_[4] + areas_[3]) / 2 +
                         (areas_[4] - areas_[3]) / 2 *
                             std::cos(M_PI * std::pow((positions_[4] - evaluation_positions[i] * xscale) /
                                                          (positions_[4] - positions_[3]),
                                                      etas_[3]));
            } else {
                out[i] = areas_[5];
            }
            // Bound the area
            out[i] = std::max(T(1e-8), out[i]);
        }
    };
};
} // namespace tarte
