#include "duffing.h"

namespace tarte {

template<typename T>
Duffing<T>::Duffing(float sampleRate)
{
    ReinitDsp(sampleRate);
};

template<typename T>
void Duffing<T>::ReinitDsp(float sampleRate)
{
    sr_ = sampleRate;
    dt_ = 1 / (float(sr_));

    // Reinit state
    q_next_ = q_now_ = q_last_ = 0;

    r_now_ = 0;

    rhs_ = lhs_ = q_now_;

    // Reinit nonlinear variables
    g_ = q_now_;
    f_nl_ = q_now_;
    E_nl_ = 0;

    // Reinit system matrices
    mass_ = stiffness_ = dissipation_ = 0;

    amplitude_ = mass_;
    pulsation_ = 2 * M_PI * 100 * amplitude_;
    decay_time_ = mass_;
};

template<typename T>
void Duffing<T>::Process(T input_force)
{
    // Nonlinear part
    ComputeVAndVprime(q_now_);
    g_ = f_nl_ / (sqrt(2 * E_nl_) + kNumEps);
    if (E_nl_ > max_E_nl_) {
        max_E_nl_ = E_nl_;
    }

    if (ctrl_term_) {
        ComputeV((q_now_ + q_last_) / 2);
        epsilon_ = r_now_ - sqrt(2 * E_nl_);
        g_ -= lambda_ctrl_ * epsilon_ * dt_ * ((T(0) < (q_now_ - q_last_)) - ((q_now_ - q_last_) < T(0))) /
              (abs(q_now_ - q_last_) + kNumEps);
    }

    // Linear part
    rhs_ = (-stiffness_ + 2 * mass_ / (dt_ * dt_)) * (q_now_) -
           (mass_ / (dt_ * dt_) - (dissipation_) / (2 * dt_)) * (q_last_) + input_force;

    // Nonlinear part
    rhs_ += 0.25 * g_ * g_ * q_last_ - g_ * r_now_;

    // Solve using Shermann-Morrisson
    auto A0_inv = dt_ * dt_ / (mass_ + dt_ * (dissipation_) / 2);
    q_next_ = A0_inv * rhs_ - (A0_inv * g_ * A0_inv * g_ * rhs_) / (4 + A0_inv * g_ * g_);

    r_last_ = r_now_;
    r_now_ = r_now_ + 0.5 * g_ * (q_next_ - q_last_);

    q_last_ = q_now_;
    q_now_ = q_next_;
};

// Intermediary and internal solver functions
template<typename T>
void Duffing<T>::ComputePhysicalParameters()
{
    dissipation_ = 1 / pulsation_;
    mass_ = 1 / (2 * 6.9) * decay_time_ * dissipation_;
    stiffness_ = mass_ * pulsation_ * pulsation_;
};

template<typename T>
void Duffing<T>::ComputeVAndVprime(T q)
{
    E_nl_ = stiffness_ * eta_ / 4 * pow(q, 4);
    f_nl_ = stiffness_ * eta_ * pow(q, 3);
};

template<typename T>
void Duffing<T>::ComputeV(T q)
{
    E_nl_ = stiffness_ * eta_ / 4 * pow(q, 4);
};

// Setters and getters
template<typename T>
void Duffing<T>::set_physical_parameters(T mass, T stiffness, T dissipation, T eta_nl)
{
    dissipation_ = dissipation;
    mass_ = mass;
    stiffness_ = stiffness;
    eta_ = eta_nl;
};

template<typename T>
void Duffing<T>::set_linear_parameters(T amplitude, T pulsation, T decay_time)
{
    amplitude_ = amplitude;
    pulsation_ = pulsation;
    decay_time_ = decay_time;

    ComputePhysicalParameters();
};

template<typename T>
void Duffing<T>::set_amplitude(T amplitude)
{
    amplitude_ = amplitude;
    ComputePhysicalParameters();
};

template<typename T>
void Duffing<T>::set_frequency(T Freq)
{
    pulsation_ = 2 * M_PI * Freq;
    ComputePhysicalParameters();
};

template<typename T>
void Duffing<T>::set_decay_time(T decay_time)
{
    decay_time_ = decay_time;
    ComputePhysicalParameters();
};

template class Duffing<double>;
template class Duffing<float>;

} // namespace tarte