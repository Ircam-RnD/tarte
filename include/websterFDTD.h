#pragma once

#include <articulation.h>
#include <utility/biquad.h>

#include <Eigen/Dense>
#include <algorithm>

namespace tarte {

// Maximum number of discretization elements supported.
// All arrays are pre-allocated to this size to avoid dynamic allocation
// during real-time processing and to make DspSetup / Process safe to
// call from different threads without heap-reallocation races.
static constexpr int kMaxN = 40;

template<typename ftype>
class WebsterFDTD {
private:
    // Size N+1  →  kMaxN + 1
    using ArrayNp1 = Eigen::Array<ftype, kMaxN + 1, 1>;
    // Size N    →  kMaxN
    using ArrayN = Eigen::Array<ftype, kMaxN, 1>;
    // Size N-1  →  kMaxN - 1
    using ArrayNm1 = Eigen::Array<ftype, kMaxN - 1, 1>;

    // Physical parameters
    ftype c0_{340}, rho0_{1.2}, l0_{17e-2};    // Acoustic
    ftype lip_radius{0}, L_rad_{0}, R_rad_{0}; // Radiation
    bool yielding_walls{false};
    ftype wall_area_mass_{15}, wall_area_damping_{16000};  // Yielding walls, per-area values
    ArrayN wall_mass_, wall_dissipation_, wall_stiffness_; // Yielding walls

    // Articulation
    ArrayNp1 S_direct_, S_target_, S_direct_last_;
    ArrayNm1 S_dual_;
    ArrayN S_primal_, S_primal_last_, d_S_primal_;
    ArrayN G_lips_, G_glottis_;

    // Discretization parameters
    ftype dt_{0}, sr_{0}, h_{0};
    ArrayNp1 x_direct_;
    ArrayNm1 x_dual_;
    ArrayN x_primal_;
    int N_{0}; // Active number of elements; always <= kMaxN

    // State variables
    ArrayN rho_now_, rho_next_;                               // Acoustic
    ArrayNm1 vel_;                                            // Acoustic velocity
    ArrayN wall_displacement_, wall_vel_now_, wall_vel_next_; // Walls
    ftype radiation_flow{0};                                  // Radiation

    // LPF  (N_lpf_ <= kMaxN + 1)
    std::array<Biquad, kMaxN + 1> lp_filters_;
    float lpf_frequency_{10.0f};
    int N_lpf_{kMaxN + 1};

    // FDTD coefficient arrays
    ArrayN intermediary_;
    ArrayN d_plus_v_, A_, B_, D_, E_, F_, G_, A_rad_, B_rad_;
    ArrayNm1 C_top_, C_low_;

    // Private helpers
    void SetNStability();
    void ComputeDiscreteGreometry();
    bool time_varying_geometry_{false};
    void UpdateWallParameters();
    void UpdateRadiationParameters();
    void UpdateCoefficients();

    void filterSdirectTarget();
    void initializeLPFStates();
    void set_N_lpf(int num);

public:
    WebsterFDTD(ftype sampleRate, ftype length = ftype(17e-2));

    // Recomputes everything (no heap allocation occurs here)
    void DspSetup(ftype sampleRate);

    // Geometry setters
    template<typename intype>
    void SetTargetGeometry(intype const* in, std::size_t const size)
    {
        // Only write into the active N+1 segment
        std::size_t safe_size = std::min(size, std::size_t(N_ + 1));
        for (std::size_t i = 0; i < safe_size; ++i) {
            S_target_[i] = static_cast<ftype>(std::max(float(1e-8), float(in[i])));
        }
    }
    void SetTargetGeometryFromArticulation(Articulation articulation);
    void SetConstantSection(ftype section);

    // Geometry smoothing
    void set_lp_frequency(int index, ftype freq);
    void set_lp_Q(int index, ftype Q);
    void set_lp_frequencies(ftype freq);
    void set_lp_Qs(ftype Q);

    // DSP
    void Process(ftype inputFlow);

    std::tuple<ftype, ftype> GetIOLinearDependencyCoefficients();

    // Listeners
    inline ftype ReadInputPressure() { return c0_ * c0_ * rho_now_(0); }
    inline ftype ReadRadiatedPressure() { return c0_ * c0_ * rho_now_(N_ - 1); }

    // Getters
    inline std::size_t get_N() { return static_cast<std::size_t>(N_); }
    inline ftype get_c0() { return c0_; }
    inline ftype get_rho0() { return rho0_; }
    inline ftype get_l0() { return l0_; }
    inline ftype get_wall_area_mass() { return wall_area_mass_; }
    inline ftype get_wall_area_damping() { return wall_area_damping_; }
    inline float get_lpf_frequency() { return lpf_frequency_; }

    // Setters
    inline void set_yielding_walls(bool isYielding)
    {
        yielding_walls = isYielding;
        UpdateCoefficients();
    }
    inline void set_time_varying_geometry(bool isVarying) { time_varying_geometry_ = isVarying; }
    void set_c0(ftype sound_velocity)
    {
        c0_ = sound_velocity;
        DspSetup(sr_);
    }
    void set_l0(ftype length)
    {
        l0_ = length;
        DspSetup(sr_);
    }
    void set_rho0(ftype rest_density)
    {
        rho0_ = rest_density;
        UpdateRadiationParameters();
        UpdateCoefficients();
    }
};

} // namespace tarte
