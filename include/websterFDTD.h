#pragma once

#include <articulation.h>
#include <utility/biquad.h>

#include <Eigen/Dense>

namespace tarte {

template<typename ftype>
class WebsterFDTD {
private:
    // Physical parameters
    ftype c0_{340}, rho0_{1.2}, l0_{17e-2};    // Acoustic
    ftype lip_radius{0}, L_rad_{0}, R_rad_{0}; // Radiation
    bool yielding_walls{false};
    ftype wall_area_mass_{15}, wall_area_damping_{16000};                      // Yielding walls, per area values
    Eigen::Array<ftype, -1, 1> wall_mass_, wall_dissipation_, wall_stiffness_; // Yielding walls

    // Articulation
    Eigen::Array<ftype, -1, 1> S_direct_, S_dual_, S_primal_, G_lips_, G_glottis_;
    Eigen::Array<ftype, -1, 1> S_direct_last_, S_primal_last_;
    Eigen::Array<ftype, -1, 1> S_target_;
    Eigen::Array<ftype, -1, 1> d_S_primal_;
    float m_lambdaS{10};

    // Discretization parameters
    ftype dt_{0}, sr_{0}, h_{0};
    Eigen::Array<ftype, -1, 1> x_primal_, x_dual_, x_direct_;
    int N_;

    // State variables
    Eigen::Array<ftype, -1, 1> rho_now_, rho_next_, vel_;                         // Acoustic
    Eigen::Array<ftype, -1, 1> wall_displacement_, wall_vel_now_, wall_vel_next_; // Walls
    ftype radiation_flow;                                                         // Radiation

    // LPF
    std::vector<Biquad> lp_filters_;
    int N_lpf_{4}; // Default number of LPF instances
    void set_N_lpf(int num);

    // Update coefficients
    Eigen::Array<ftype, -1, 1> intermediary_;
    Eigen::Array<ftype, -1, 1> d_plus_v_, A_, B_, C_top_, C_low_, D_, E_, F_, G_, A_rad_, B_rad_;

    // Set the number of discretization elements to be at the stability condition
    void SetNStability();
    // Computes the discrete geometry vectors
    void ComputeDiscreteGreometry();
    bool time_varying_geometry_{false};
    // Recomputes discrete wall parameters from current geometry
    void UpdateWallParameters();
    // Recomputes radiation parameters from current geometry
    void UpdateRadiationParameters();
    // Recomputes intermediary variables from current geometry
    void UpdateCoefficients();

    void filterSdirectTarget();
    void initializeLPFStates();

public:
    WebsterFDTD(float sampleRate, float length = 1.0f);

    // Recomputes everything
    void DspSetup(float sampleRate);

    // Setters for target geometry
    void SetTargetGeometryFromArticulation(Articulation articulation);
    void SetConstantSection(ftype section);

    // Geometry filtering functions
    void set_lp_frequency(int index, float freq);
    void set_lp_Q(int index, float freq);
    void set_lp_frequencies(float freq);
    void set_lp_Qs(float freq);

    // Process a sample
    void Process(ftype inputFlow);

    // Linear dependency of the next input pressure (conjugated output) on the
    // next input flow (input)
    std::tuple<ftype, ftype> GetIOLinearDependencyCoefficients();

    // "Listening" functions
    inline ftype ReadInputPressure() { return c0_ * c0_ * rho_now_(0); }
    inline ftype ReadRadiatedPressure() { return c0_ * c0_ * rho_now_(N_ - 1); }

    // Getters
    inline ftype get_c0() { return c0_; }
    inline ftype get_rho0() { return rho0_; }
    inline ftype get_l0() { return l0_; }
    inline ftype get_wall_area_mass() { return wall_area_mass_; }
    inline ftype get_wall_area_damping() { return wall_area_damping_; }

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
