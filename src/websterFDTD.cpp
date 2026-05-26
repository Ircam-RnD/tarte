#include "websterFDTD.h"

#include <iostream>

namespace tarte {

template<typename ftype, int kMaxN>
WebsterFDTD<ftype, kMaxN>::WebsterFDTD(ftype sampleRate, ftype length, Articulation* art)
{
    l0_ = length;
    DspSetup(sampleRate, art);
}

// DspSetup – no heap allocation; all arrays are fixed-size
template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::DspSetup(ftype sampleRate, Articulation* art)
{
    sr_ = sampleRate;
    dt_ = 1 / sr_;

    c02_ = c0_ * c0_; // squared sound velocity

    // Determine N (clamped to kMaxN)
    SetNStability();

    vel_coeff_ = c02_ / rho0_ / h_ * dt_;
    // Spatial grids (write only into active segment)
    x_direct_.setZero();
    x_dual_.setZero();
    x_direct_.head(N_ + 1) = Eigen::Array<ftype, -1, 1>::LinSpaced(N_ + 1, -h_ / 2, l0_ + h_ / 2);
    // x_primal_.head(N_) = Eigen::Array<ftype, -1, 1>::LinSpaced(N_, 0, l0_);
    // x_dual_ is a view into x_direct_[1 .. N-1]; we replicate that here:
    x_dual_.head(N_ - 1) = x_direct_.segment(1, N_ - 1);

    // Sections
    S_target_.setZero();
    S_direct_.setZero();
    S_primal_.setZero();
    S_direct_last_.setZero();
    S_primal_last_.setZero();
    d_S_primal_.setZero();
    gamma_primal_.setZero();
    if (art) {
        SetTargetGeometryFromArticulation(*art, true);
    } else {
        S_target_.head(N_ + 1).setOnes();
        S_direct_.head(N_ + 1) = S_target_.head(N_ + 1);
        S_direct_last_.head(N_ + 1) = S_direct_.head(N_ + 1);
        ComputeDiscreteGreometry();
    }
    S_primal_last_.head(N_) = S_primal_.head(N_);

    // Radiation
    UpdateRadiationParameters();

    // Intermediary arrays reset
    d_plus_v_.setZero();
    intermediary_.setZero();
    A_.setZero();
    B_.setZero();
    C_top_.setZero();
    C_low_.setZero();
    D_.setZero();
    E_.setZero();
    A_rad_.setZero();
    B_rad_.setZero();

    // State reset
    flip_ = false;
    rho_now_ac().setZero();
    rho_next_ac().setZero();
    vel_.setZero();
    wall_displacement_.setZero();
    wall_momentum_now_ac().setZero();
    wall_momentum_next_ac().setZero();
    radiation_flow = 0;

    UpdateCoefficients();

    // LFPs
    N_lpf_ = N_ + 1; // one filter per direct-grid point
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i] = Biquad(sr_, kLowPass, lpf_frequency_, 0.0f, 0.5f);
        if (i < N_ + 1) {
            lp_filters_[i].InitializeState(static_cast<double>(S_target_[i]));
        }
    }
    set_lp_Qs(0.5);
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::SetNStability()
{
    h_ = c0_ / sr_;
    N_ = static_cast<int>(std::floor(l0_ / h_));

    if (N_ > kMaxN) {
        N_ = kMaxN;
    }
    if (N_ < 2) {
        N_ = 2;
    }

    h_ = l0_ / N_;
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::SetTargetGeometryFromArticulation(Articulation articulation, bool force_direct)
{
    articulation.getAreas(x_direct_.data(), S_target_.data(), N_ + 1);
    if (!time_varying_geometry_ or force_direct) {
        S_direct_.head(N_ + 1) = S_target_.head(N_ + 1);
        ComputeDiscreteGreometry();
        UpdateRadiationParameters();
        UpdateCoefficients();
    }
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::SetConstantSection(ftype section)
{
    S_target_.head(N_ + 1).setConstant(section);
    if (!time_varying_geometry_) {
        S_direct_.head(N_ + 1) = S_target_.head(N_ + 1);
        ComputeDiscreteGreometry();
        UpdateRadiationParameters();
        UpdateCoefficients();
    }
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::ComputeDiscreteGreometry()
{
    S_direct_(0) = S_direct_(1);
    S_direct_(N_) = S_direct_(N_ - 1);

    // S_dual_[0 .. N-2]  = S_direct_[1 .. N-1]
    S_dual_.head(N_ - 1) = S_direct_.segment(1, N_ - 1);

    // S_primal_[0 .. N-1] = 0.5 * (S_direct_[0..N-1] + S_direct_[1..N])
    S_primal_.head(N_) = 0.5 * (S_direct_.head(N_) + S_direct_.segment(1, N_));
    gamma_primal_.head(N_) =
        2 * S_primal_.head(N_).sqrt() * static_cast<ftype>(std::sqrt(M_PI)); // Assume circular shape
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::UpdateRadiationParameters()
{
    ftype Aout = S_primal_[N_ - 1];
    lip_radius = std::sqrt(Aout / static_cast<ftype>(M_PI));
    ftype Z0 = rho0_ * c0_ / Aout;
    L_rad_ = 8 * lip_radius / (3 * static_cast<ftype>(M_PI) * c0_) * Z0;
    R_rad_ = 128 / (9 * static_cast<ftype>(M_PI) * static_cast<ftype>(M_PI)) * Z0;
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::UpdateCoefficients()
{
    // Convenience aliases onto active segments (avoid repeating .head(N_) everywhere)
    auto Sp = S_primal_.head(N_);
    auto Sd = S_dual_.head(N_ - 1);
    auto gamma = gamma_primal_.head(N_);
    auto mw = wall_area_mass_;
    auto rw = wall_area_damping_;
    auto kw = wall_area_stiffness_;

    ftype rhoc2 = rho0_ * c02_;

    if (yielding_walls) {
        intermediary_.head(N_) = 2 / dt_ + rw / mw + kw / (2 * mw) * dt_;

        A_.head(N_) = 1 / dt_ + gamma * rhoc2 / (2 * mw * intermediary_.head(N_));
        A_(N_ - 1) += rhoc2 / (Sp(N_ - 1) * 2 * h_ * R_rad_) + rhoc2 * dt_ / (Sp(N_ - 1) * 4 * h_ * L_rad_);

        D_.head(N_) = -2 * rho0_ * gamma / (dt_ * mw * intermediary_.head(N_) * Sp);
        E_.head(N_) = rho0_ * gamma / (intermediary_.head(N_) * Sp) * kw / mw;

        A_rad_.head(N_) = intermediary_.head(N_) * 0.5;
        B_rad_.head(N_) = -A_rad_.head(N_) + 2 / dt_;

    } else {
        A_.head(N_) = 1 / dt_;
        A_(N_ - 1) += rhoc2 / (Sp(N_ - 1) * 2 * h_ * R_rad_) + rhoc2 * dt_ / (Sp(N_ - 1) * 4 * h_ * L_rad_);
    }

    B_.head(N_) = -A_.head(N_) + 2 / dt_;
    C_top_.head(N_ - 1) = -1 / Sp.head(N_ - 1) * rho0_ / h_ * Sd;
    C_low_.head(N_ - 1) = 1 / Sp.tail(N_ - 1) * rho0_ / h_ * Sd;

    F_ = -rho0_ / (h_ * Sp(N_ - 1)); // lips boundary
    G_ = rho0_ / (h_ * Sp(0));       // Glottis boudary
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::Process(ftype inputFlow)
{
    // These views should not have any cost
    auto rho_now = rho_now_ac().head(N_);
    auto rho_next = rho_next_ac().head(N_);
    auto vel = vel_.head(N_ - 1);

    auto A = A_.head(N_);
    auto B = B_.head(N_);
    auto C_top = C_top_.head(N_ - 1);
    auto C_low = C_low_.head(N_ - 1);
    auto D = D_.head(N_);
    auto E = E_.head(N_);
    auto A_rad = A_rad_.head(N_);
    auto B_rad = B_rad_.head(N_);

    auto dv = d_plus_v_.head(N_);
    auto dSp = d_S_primal_.head(N_);
    auto Sp = S_primal_.head(N_);

    auto mw = wall_area_mass_;
    auto kw = wall_area_stiffness_;
    auto wdisp = wall_displacement_.head(N_);
    auto wp_now = wall_momentum_now_ac().head(N_);
    auto wp_next = wall_momentum_next_ac().head(N_);

    if (yielding_walls) {
        dv.head(N_ - 1) = C_top * vel;
        dv(N_ - 1) = 0;
        dv.tail(N_ - 1) += C_low * vel;

        rho_next = (1 / A) * (B * rho_now + dv + D * wp_now + E * wdisp - rho0_ * (dSp / Sp));

        rho_next(0) += G_ * inputFlow / A(0);
        rho_next(N_ - 1) += F_ * radiation_flow / A(N_ - 1);

        vel = vel - vel_coeff_ * (rho_next.tail(N_ - 1) - rho_next.head(N_ - 1));

        wp_next = (1 / A_rad) * (B_rad * wp_now - kw * wdisp + (c02_ * ftype(0.5)) * (rho_now + rho_next));

        wdisp += dt_ * ftype(0.5) / wall_area_mass_ * (wp_now + wp_next);
        radiation_flow += dt_ * c02_ / L_rad_ * ftype(0.5) * (rho_next(N_ - 1) + rho_now(N_ - 1));

    } else {
        dv.head(N_ - 1) = C_top * vel;
        dv(N_ - 1) = 0;
        dv.tail(N_ - 1) += C_low * vel;

        rho_next = (1 / A) * (B * rho_now + dv - rho0_ * (1 / Sp * dSp));
        rho_next(0) += G_ * inputFlow / A(0);
        rho_next(N_ - 1) += F_ * radiation_flow / A(N_ - 1);

        vel = vel - vel_coeff_ * (rho_next.tail(N_ - 1) - rho_next.head(N_ - 1));

        radiation_flow += dt_ * c02_ / L_rad_ * ftype(0.5) * (rho_next(N_ - 1) + rho_now(N_ - 1));
    }

    // Swap buffers to advance state
    flip_ = !flip_;

    if (time_varying_geometry_) {
        // ~19 ms total
        for (int i = 0; i < N_ + 1 && i < N_lpf_; ++i) {
            S_direct_[i] = static_cast<ftype>(lp_filters_[i].Process(static_cast<double>(S_target_[i])));
        } // ~9ms

        ComputeDiscreteGreometry();  // ~1ms
        UpdateRadiationParameters(); // negligible
        UpdateCoefficients();        // ~ 7 ms

        S_direct_last_.head(N_ + 1) = S_direct_.head(N_ + 1);
        S_primal_last_.head(N_) = S_primal_.head(N_);
        d_S_primal_.head(N_) = (S_primal_.head(N_) - S_primal_last_.head(N_)) / dt_;
    }
}

template<typename ftype, int kMaxN>
std::tuple<ftype, ftype> WebsterFDTD<ftype, kMaxN>::GetIOLinearDependencyCoefficients()
{
    return {c02_ * ftype(0.5) *
                (rho_now_ac()(0) + 1 / A_(0) *
                                       (B_(0) * rho_now_ac()(0) - rho0_ / h_ * S_dual_(0) / S_primal_(0) * vel_(0) -
                                        rho0_ * (1 / S_primal_(0) * d_S_primal_(0)))),
            ftype(0.5) * c02_ * (1 / A_(0)) * G_};
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::set_N_lpf(int num)
{
    N_lpf_ = std::clamp(num, 1, kMaxN + 1);
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i] = Biquad(sr_, kLowPass, 10.0f, 0.0f, 0.7f);
    }
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::set_lp_frequency(int index, ftype freq)
{
    if (index >= 0 && index < N_lpf_) {
        lp_filters_[index].set_freq(freq);
    }
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::set_lp_Q(int index, ftype Q)
{
    if (index >= 0 && index < N_lpf_) {
        lp_filters_[index].set_Q(Q);
    }
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::set_lp_frequencies(ftype freq)
{
    lpf_frequency_ = std::clamp(freq, ftype(0.1), ftype(100.0));
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i].set_freq(lpf_frequency_);
    }
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::set_lp_Qs(ftype Q)
{
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i].set_Q(Q);
    }
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::filterSdirectTarget()
{
    for (int i = 0; i < N_ + 1 && i < N_lpf_; ++i) {
        S_target_[i] = static_cast<ftype>(lp_filters_[i].Process(static_cast<double>(S_target_[i])));
    }
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::initializeLPFStates()
{
    for (int i = 0; i < N_lpf_ && i < N_ + 1; ++i) {
        lp_filters_[i].InitializeState(static_cast<double>(S_target_[i]));
    }
}

// Default to kMaxN = 50, enough for voice at 96000 Hz
template class WebsterFDTD<float>;
template class WebsterFDTD<double>;

// also allow for kMaxN = 300, sould be enough for every musical resonator.
template class WebsterFDTD<float, 300>;
template class WebsterFDTD<double, 300>;

} // namespace tarte
