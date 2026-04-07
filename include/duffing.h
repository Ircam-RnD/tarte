#pragma once

#include <algorithm>
#include <cmath>
#include <tuple>

#include <export.h>

// TODO : change integration to exact for the linear part
namespace tarte {

template<typename T>
class Duffing {
private:
    // Numerical epsilon value
    constexpr static T kNumEps{1e-14};

    // Linear part: system parameters
    T mass_, stiffness_, dissipation_;
    // Higer level modal parameters
    T amplitude_, pulsation_, decay_time_;

    // Nonlinear parameters
    T eta_{0};

    // Time-scheme parameters
    T sr_;
    T dt_;
    bool ctrl_term_{true};
    T lambda_ctrl_{0};

    // System state
    T q_last_, q_now_, q_next_;
    T r_now_, r_last_;

    // Intermediate vectors
    T rhs_, lhs_;

    // Nonlinear function evaluation
    T g_, f_nl_;
    T E_nl_;

    // Drift variable
    T epsilon_{0}, max_E_nl_{0};

    void ComputePhysicalParameters();

    void ComputeVAndVprime(T q);

    void ComputeV(T q);

public:
    Duffing(float sample_rate);

    void ReinitDsp(float sample_rate);

    void Process(T input_force);

    // "Listening" functions
    T ReadDisplacement() { return q_now_; };
    T ReadVelocity() { return (q_next_ - q_last_) / (2 * dt_); };

    // Getters and setters
    // Physical parameters
    void set_physical_parameters(T mass, T stiffness, T dissipation, T eta_nl = 0);
    // Higher level modal parameters
    void set_linear_parameters(T amplitude, T pulsation, T decay_time);
    void set_amplitude(T amplitude);
    void set_frequency(T frequency);
    void set_decay_time(T decay_time);
    void set_nonlinearity(T eta) { eta_ = eta; };
    // SAV drift control term (should be useless for duffing if parameters are not modulated in time)
    void set_lambda_ctrl(T lambda_ctrl) { lambda_ctrl_ = std::clamp(lambda_ctrl, T(0), T(10000)); };
};

} // namespace tarte