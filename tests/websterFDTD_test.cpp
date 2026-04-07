#include <cmath>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

#include "articulation.h"
#include "websterFDTD.h" // adjust path as needed

using namespace tarte;

// ──────────────────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────────────────

static constexpr float kSR = 44100.0f;
static constexpr double kSRd = 44100.0;
static constexpr float kLen = 0.17f; // 17 cm, typical vocal-tract length

/// Run N silent (zero-input) samples and return the last radiated pressure.
template<typename F>
F runSilent(WebsterFDTD<F>& wt, int n = 512)
{
    F last = 0;
    for (int i = 0; i < n; ++i) {
        wt.Process(F{0});
        last = wt.ReadRadiatedPressure();
    }
    return last;
}

/// Check that a value is finite (not NaN, not Inf).
template<typename F>
void expectFinite(F v, const char* label = "")
{
    EXPECT_TRUE(std::isfinite(v)) << label << " is not finite: " << v;
}

// ──────────────────────────────────────────────────────────────────────────────
// 1. Construction & basic setup
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, ConstructsWithDefaultLength)
{
    ASSERT_NO_THROW(WebsterFDTD<float> wt(kSR));
}

TEST(WebsterFDTD, ConstructsWithExplicitLength)
{
    ASSERT_NO_THROW(WebsterFDTD<float> wt(kSR, kLen));
}

TEST(WebsterFDTD, DoubleTemplateCompiles)
{
    ASSERT_NO_THROW(WebsterFDTD<double> wt(kSRd, 0.17));
}

TEST(WebsterFDTD, DefaultPhysicalParameters)
{
    WebsterFDTD<float> wt(kSR);
    EXPECT_FLOAT_EQ(wt.get_c0(), 340.0f);
    EXPECT_FLOAT_EQ(wt.get_rho0(), 1.2f);
    EXPECT_FLOAT_EQ(wt.get_l0(), 17e-2f);
}

TEST(WebsterFDTD, DefaultWallParameters)
{
    WebsterFDTD<float> wt(kSR);
    EXPECT_FLOAT_EQ(wt.get_wall_area_mass(), 15.0f);
    EXPECT_FLOAT_EQ(wt.get_wall_area_damping(), 16000.0f);
}

// ──────────────────────────────────────────────────────────────────────────────
// 2. Setters – physical parameters
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, SetSpeedOfSound)
{
    WebsterFDTD<float> wt(kSR);
    wt.set_c0(350.0f);
    EXPECT_FLOAT_EQ(wt.get_c0(), 350.0f);
}

TEST(WebsterFDTD, SetLength)
{
    WebsterFDTD<float> wt(kSR);
    wt.set_l0(0.20f);
    EXPECT_FLOAT_EQ(wt.get_l0(), 0.20f);
}

TEST(WebsterFDTD, SetDensity)
{
    WebsterFDTD<float> wt(kSR);
    wt.set_rho0(1.3f);
    EXPECT_FLOAT_EQ(wt.get_rho0(), 1.3f);
}

// ──────────────────────────────────────────────────────────────────────────────
// 3. Constant-section geometry
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, SetConstantSectionPositive)
{
    WebsterFDTD<float> wt(kSR, kLen);
    ASSERT_NO_THROW(wt.SetConstantSection(1e-4f)); // 1 cm² — a narrow tube
}

TEST(WebsterFDTD, SetConstantSectionLarger)
{
    WebsterFDTD<float> wt(kSR, kLen);
    ASSERT_NO_THROW(wt.SetConstantSection(5e-4f));
}

// After setting a constant section, a silent input should produce finite output.
TEST(WebsterFDTD, ConstantSectionSilentRun)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    expectFinite(runSilent(wt), "radiated pressure (constant section, silent)");
}

// ──────────────────────────────────────────────────────────────────────────────
// 4. Process() – numerical stability
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, ProcessSilenceYieldsFiniteOutput)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    for (int i = 0; i < 2048; ++i) {
        wt.Process(0.0f);
        expectFinite(wt.ReadRadiatedPressure(), "radiated pressure");
        expectFinite(wt.ReadInputPressure(), "input pressure");
    }
}

TEST(WebsterFDTD, ProcessImpulseYieldsFiniteOutput)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    wt.Process(1e-5f); // single impulse
    for (int i = 0; i < 2048; ++i) {
        wt.Process(0.0f);
        expectFinite(wt.ReadRadiatedPressure(), "radiated pressure post-impulse");
    }
}

TEST(WebsterFDTD, ProcessSilenceYieldsNearZeroOutput)
{
    // A resting tube driven by zero input should settle to (near) zero.
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    float p_ = runSilent(wt, 8192);
    EXPECT_NEAR(p_, 0.0f, 1e-6f);
}

// ──────────────────────────────────────────────────────────────────────────────
// 5. Linearity: superposition
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, SuperpositionHolds)
{
    // p_(a+b) ≈ p_(a) + p_(b) for a linear system
    const double A = 1e-5f, B = 3e-5f;
    const int N = 128;

    // Simulation for input A
    WebsterFDTD<double> wtA(kSR, kLen);
    wtA.SetConstantSection(2e-4f);
    std::vector<double> pA(N);
    for (int i = 0; i < N; ++i) {
        wtA.Process(i == 0 ? A : 0.0f);
        pA[i] = wtA.ReadRadiatedPressure();
    }

    // Simulation for input B
    WebsterFDTD<double> wtB(kSR, kLen);
    wtB.SetConstantSection(2e-4f);
    std::vector<double> pB(N);
    for (int i = 0; i < N; ++i) {
        wtB.Process(i == 0 ? B : 0.0f);
        pB[i] = wtB.ReadRadiatedPressure();
    }

    // Simulation for combined input A+B
    WebsterFDTD<double> wtAB(kSR, kLen);
    wtAB.SetConstantSection(2e-4f);
    std::vector<double> pAB(N);
    for (int i = 0; i < N; ++i) {
        wtAB.Process(i == 0 ? A + B : 0.0f);
        pAB[i] = wtAB.ReadRadiatedPressure();
    }

    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(pAB[i], pA[i] + pB[i], 1e-8f) << "Superposition fails at sample " << i;
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 6. IO linear dependency coefficients
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, IOCoefficientsAreFinite)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    auto [alpha, beta] = wt.GetIOLinearDependencyCoefficients();
    expectFinite(alpha, "alpha");
    expectFinite(beta, "beta");
}

// ──────────────────────────────────────────────────────────────────────────────
// 7. Yielding walls
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, YieldingWallsCanBeToggled)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    ASSERT_NO_THROW(wt.set_yielding_walls(true));
    ASSERT_NO_THROW(wt.set_yielding_walls(false));
}

TEST(WebsterFDTD, YieldingWallsNumericallyStable)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    wt.set_yielding_walls(true);
    wt.Process(1e-5f);
    for (int i = 0; i < 2048; ++i) {
        wt.Process(0.0f);
        expectFinite(wt.ReadRadiatedPressure(), "yielding-wall radiated pressure");
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 8. Time-varying geometry
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, TimeVaryingGeometryToggle)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    ASSERT_NO_THROW(wt.set_time_varying_geometry(true));
    ASSERT_NO_THROW(wt.set_time_varying_geometry(false));
}

// ──────────────────────────────────────────────────────────────────────────────
// 9. LPF geometry filtering
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, SetLpFrequencyDoesNotThrow)
{
    WebsterFDTD<float> wt(kSR, kLen);
    ASSERT_NO_THROW(wt.set_lp_frequency(0, 1000.0f));
}

TEST(WebsterFDTD, SetLpQDoesNotThrow)
{
    WebsterFDTD<float> wt(kSR, kLen);
    ASSERT_NO_THROW(wt.set_lp_Q(0, 0.707f));
}

TEST(WebsterFDTD, SetAllLpFrequencies)
{
    WebsterFDTD<float> wt(kSR, kLen);
    ASSERT_NO_THROW(wt.set_lp_frequencies(500.0f));
}

TEST(WebsterFDTD, SetAllLpQs)
{
    WebsterFDTD<float> wt(kSR, kLen);
    ASSERT_NO_THROW(wt.set_lp_Qs(0.5f));
}

// ──────────────────────────────────────────────────────────────────────────────
// 10. DspSetup re-initialization
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, DspSetupWithSameSampleRate)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    ASSERT_NO_THROW(wt.DspSetup(kSR));
}

TEST(WebsterFDTD, DspSetupWithNewSampleRate)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    ASSERT_NO_THROW(wt.DspSetup(48000.0f));
    // Physical parameters must remain unchanged after re-setup
    EXPECT_FLOAT_EQ(wt.get_c0(), 340.0f);
}

TEST(WebsterFDTD, DspSetupResetsState)
{
    // After an impulse, re-setup should clear internal state.
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    wt.Process(1.0f);
    float pBefore = wt.ReadRadiatedPressure();

    wt.DspSetup(kSR);
    wt.SetConstantSection(2e-4f);
    float pAfter = wt.ReadRadiatedPressure();

    // Post-reset pressure should be closer to zero than the mid-ring value.
    EXPECT_LT(std::abs(pAfter), std::abs(pBefore) + 1e-9f);
}

// ──────────────────────────────────────────────────────────────────────────────
// 11. Articulation-based geometry
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, SetTargetFromDefaultArticulation)
{
    WebsterFDTD<float> wt(kSR, kLen);
    Articulation art; // default-constructed; assumes Articulation has a default ctor
    ASSERT_NO_THROW(wt.SetTargetGeometryFromArticulation(art));
}

TEST(WebsterFDTD, ArticulationGeometryNumericallyStable)
{
    WebsterFDTD<float> wt(kSR, kLen);
    Articulation art;
    wt.SetTargetGeometryFromArticulation(art);
    wt.set_time_varying_geometry(true);
    wt.Process(1e-5f);
    for (int i = 0; i < 1024; ++i) {
        wt.Process(0.0f);
        expectFinite(wt.ReadRadiatedPressure(), "articulation radiated pressure");
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 12. Output pressure – sign / causality sanity
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, PositiveFlowEventuallyProducesNonZeroOutput)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);

    // Drive with a positive flow impulse.
    wt.Process(1e-4f);

    // The radiated pressure should become nonzero within the tube's propagation
    // delay (roughly l0 / c0 * sr samples ≈ 22 samples at 44.1 kHz / 17 cm).
    bool nonzeroSeen = false;
    for (int i = 0; i < 512; ++i) {
        wt.Process(0.0f);
        if (std::abs(wt.ReadRadiatedPressure()) > 1e-10f) {
            nonzeroSeen = true;
            break;
        }
    }
    EXPECT_TRUE(nonzeroSeen) << "Radiated pressure never became nonzero after an impulse";
}

TEST(WebsterFDTD, InputPressureAndRadiatedPressureAreBothFinite)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    wt.Process(1e-5f);
    expectFinite(wt.ReadInputPressure(), "input pressure");
    expectFinite(wt.ReadRadiatedPressure(), "radiated pressure");
}

// ──────────────────────────────────────────────────────────────────────────────
// 13. Effect of changing sound speed on output
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, ChangingSoundSpeedAffectsOutput)
{
    // Two otherwise identical tubes, different c0 → different impulse responses.
    WebsterFDTD<float> wt1(kSR, kLen), wt2(kSR, kLen);
    wt1.SetConstantSection(2e-4f);
    wt2.SetConstantSection(2e-4f);
    wt2.set_c0(360.0f);

    wt1.Process(1e-5f);
    wt2.Process(1e-5f);

    std::vector<float> p1, p2;
    for (int i = 0; i < 256; ++i) {
        wt1.Process(0.0f);
        wt2.Process(0.0f);
        p1.push_back(wt1.ReadRadiatedPressure());
        p2.push_back(wt2.ReadRadiatedPressure());
    }

    float diffNorm = 0.0f;
    for (int i = 0; i < 256; ++i)
        diffNorm += std::abs(p1[i] - p2[i]);
    EXPECT_GT(diffNorm, 1e-10f) << "Different c0 produced identical impulse responses";
}

// ──────────────────────────────────────────────────────────────────────────────
// 14. Edge cases
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD, VeryLargeFlowInputRemainsFinite)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    // Intentionally large (but physical) flow magnitude.
    wt.Process(1.0f);
    for (int i = 0; i < 128; ++i) {
        wt.Process(0.0f);
        expectFinite(wt.ReadRadiatedPressure(), "large impulse radiated pressure");
    }
}

TEST(WebsterFDTD, NegativeFlowInputRemainsFinite)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    wt.Process(-1e-5f);
    for (int i = 0; i < 128; ++i) {
        wt.Process(0.0f);
        expectFinite(wt.ReadRadiatedPressure(), "negative impulse radiated pressure");
    }
}

TEST(WebsterFDTD, AlternatingSigns_NumericallyStable)
{
    WebsterFDTD<float> wt(kSR, kLen);
    wt.SetConstantSection(2e-4f);
    for (int i = 0; i < 1024; ++i) {
        wt.Process((i % 2 == 0) ? 1e-5f : -1e-5f);
        expectFinite(wt.ReadRadiatedPressure(), "alternating-sign pressure");
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// 15. Double-precision variant
// ──────────────────────────────────────────────────────────────────────────────

TEST(WebsterFDTD_Double, BasicRunIsStable)
{
    WebsterFDTD<double> wt(kSRd, 0.17);
    wt.SetConstantSection(2e-4);
    wt.Process(1e-5);
    for (int i = 0; i < 512; ++i) {
        wt.Process(0.0);
        EXPECT_TRUE(std::isfinite(wt.ReadRadiatedPressure()));
    }
}

TEST(WebsterFDTD_Double, DefaultPhysicalParameters)
{
    WebsterFDTD<double> wt(kSRd);
    EXPECT_DOUBLE_EQ(wt.get_c0(), 340.0);
    EXPECT_DOUBLE_EQ(wt.get_rho0(), 1.2);
    EXPECT_DOUBLE_EQ(wt.get_l0(), 17e-2);
}