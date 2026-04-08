#include "websterFDTD.h"

#include <iostream>

namespace tarte {

template<typename ftype>
WebsterFDTD<ftype>::WebsterFDTD(ftype sampleRate, ftype length)
{
    l0_ = length;
    DspSetup(sampleRate);
}

// DspSetup – no heap allocation; all arrays are fixed-size
template<typename ftype>
void WebsterFDTD<ftype>::DspSetup(ftype sampleRate)
{
    sr_ = sampleRate;
    dt_ = 1 / sr_;

    // Determine N (clamped to kMaxN)
    SetNStability();

    // Spatial grids (write only into active segment)
    x_direct_.head(N_ + 1) = Eigen::Array<ftype, -1, 1>::LinSpaced(N_ + 1, -h_ / 2, l0_ + h_ / 2);
    x_primal_.head(N_) = Eigen::Array<ftype, -1, 1>::LinSpaced(N_, 0, l0_);
    // x_dual_ is a view into x_direct_[1 .. N-1]; we replicate that here:
    x_dual_.head(N_ - 1) = x_direct_.segment(1, N_ - 1);

    // Sections
    S_target_.head(N_ + 1).setOnes();
    S_direct_.head(N_ + 1) = S_target_.head(N_ + 1);
    S_direct_last_.head(N_ + 1) = S_direct_.head(N_ + 1);

    ComputeDiscreteGreometry();
    S_primal_last_.head(N_) = S_primal_.head(N_);

    d_S_primal_.head(N_).setZero();

    // Boundaries
    G_lips_.head(N_).setZero();
    G_lips_(N_ - 1) = 1;
    G_glottis_.head(N_).setZero();
    G_glottis_(0) = 1;

    // Wall and radiation
    UpdateWallParameters();
    UpdateRadiationParameters();

    // Intermediary arrays reset
    d_plus_v_.head(N_).setZero();
    intermediary_.head(N_).setZero();
    A_.head(N_).setZero();
    B_.head(N_).setZero();
    C_top_.head(N_ - 1).setZero();
    C_low_.head(N_ - 1).setZero();
    D_.head(N_).setZero();
    E_.head(N_).setZero();
    F_.head(N_).setZero();
    G_.head(N_).setZero();
    A_rad_.head(N_).setZero();
    B_rad_.head(N_).setZero();

    // State reset
    rho_now_.head(N_).setZero();
    rho_next_.head(N_).setZero();
    vel_.head(N_ - 1).setZero();
    wall_displacement_.head(N_).setZero();
    wall_vel_now_.head(N_).setZero();
    wall_vel_next_.head(N_).setZero();
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

template<typename ftype>
void WebsterFDTD<ftype>::SetNStability()
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

template<typename ftype>
void WebsterFDTD<ftype>::SetTargetGeometryFromArticulation(Articulation articulation)
{
    articulation.getAreas(x_direct_.data(), S_target_.data(), N_ + 1);
    if (!time_varying_geometry_) {
        S_direct_.head(N_ + 1) = S_target_.head(N_ + 1);
        ComputeDiscreteGreometry();
        UpdateCoefficients();
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::SetConstantSection(ftype section)
{
    S_target_.head(N_ + 1).setConstant(section);
    if (!time_varying_geometry_) {
        S_direct_.head(N_ + 1) = S_target_.head(N_ + 1);
        ComputeDiscreteGreometry();
        UpdateRadiationParameters();
        UpdateCoefficients();
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::ComputeDiscreteGreometry()
{
    S_direct_(0) = S_direct_(1);
    S_direct_(N_) = S_direct_(N_ - 1);

    // S_dual_[0 .. N-2]  = S_direct_[1 .. N-1]
    S_dual_.head(N_ - 1) = S_direct_.segment(1, N_ - 1);

    // S_primal_[0 .. N-1] = 0.5 * (S_direct_[0..N-1] + S_direct_[1..N])
    S_primal_.head(N_) = 0.5 * (S_direct_.head(N_) + S_direct_.segment(1, N_));
}

template<typename ftype>
void WebsterFDTD<ftype>::UpdateWallParameters()
{
    wall_mass_.head(N_) = wall_area_mass_ * 2 * S_primal_.head(N_).sqrt() * static_cast<ftype>(std::sqrt(M_PI));
    wall_stiffness_.head(N_) = wall_mass_.head(N_) * (2 * M_PI * 70) * (2 * M_PI * 70);
    wall_dissipation_.head(N_) =
        wall_area_damping_ * 2 * S_primal_.head(N_).sqrt() * static_cast<ftype>(std::sqrt(M_PI));
}

template<typename ftype>
void WebsterFDTD<ftype>::UpdateRadiationParameters()
{
    ftype Aout = S_primal_[N_ - 1];
    lip_radius = std::sqrt(Aout / static_cast<ftype>(M_PI));
    ftype Z0 = rho0_ * c0_ / Aout;
    L_rad_ = 8 * lip_radius / (3 * static_cast<ftype>(M_PI) * c0_) * Z0;
    R_rad_ = 128 / (9 * static_cast<ftype>(M_PI) * static_cast<ftype>(M_PI)) * Z0;
}

template<typename ftype>
void WebsterFDTD<ftype>::UpdateCoefficients()
{
    // Convenience aliases onto active segments (avoid repeating .head(N_) everywhere)
    auto Sp = S_primal_.head(N_);
    auto Sd = S_dual_.head(N_ - 1);
    auto Gl = G_lips_.head(N_);
    auto Gg = G_glottis_.head(N_);
    auto wm = wall_mass_.head(N_);
    auto wd = wall_dissipation_.head(N_);
    auto ws = wall_stiffness_.head(N_);

    ftype rhoc2 = rho0_ * c0_ * c0_;

    if (yielding_walls) {
        intermediary_.head(N_) = 2 / dt_ + wd / wm + ws / (2 * wm) * dt_;

        A_.head(N_) = 1 / dt_ + 4 * static_cast<ftype>(M_PI) * rhoc2 / (2 * wm * intermediary_.head(N_)) +
                      rhoc2 * Gl / (Sp * 2 * h_ * R_rad_) + rhoc2 * Gl * dt_ / (Sp * 4 * h_ * L_rad_);

        B_.head(N_) = 1 / dt_ - 4 * static_cast<ftype>(M_PI) * rhoc2 / (2 * wm * intermediary_.head(N_)) -
                      rhoc2 * Gl / (Sp * 2 * h_ * R_rad_) - rhoc2 * Gl * dt_ / (Sp * 4 * h_ * L_rad_);

        C_top_.head(N_ - 1) = -1 / Sp.head(N_ - 1) * rho0_ / h_ * Sd;
        C_low_.head(N_ - 1) = 1 / Sp.tail(N_ - 1) * rho0_ / h_ * Sd;

        D_.head(N_) = -2 * rho0_ / (dt_ * intermediary_.head(N_) * Sp);
        E_.head(N_) = rho0_ / (intermediary_.head(N_) * Sp) * ws / wm;
        F_.head(N_) = -rho0_ * Gl / (h_ * Sp);
        G_.head(N_) = rho0_ * Gg / (h_ * Sp);

        A_rad_.head(N_) = 1 / dt_ + wd / (2 * wm) + ws / (4 * wm) * dt_;
        B_rad_.head(N_) = 1 / dt_ - wd / (2 * wm) - ws / (4 * wm) * dt_;

    } else {
        intermediary_.head(N_).setConstant(2 / dt_);

        A_.head(N_) = 1 / dt_ + rhoc2 * Gl / (Sp * 2 * h_ * R_rad_) + rhoc2 * Gl * dt_ / (Sp * 4 * h_ * L_rad_);

        B_.head(N_) = 1 / dt_ - rhoc2 * Gl / (Sp * 2 * h_ * R_rad_) - rhoc2 * Gl * dt_ / (Sp * 4 * h_ * L_rad_);

        C_top_.head(N_ - 1) = -1 / Sp.head(N_ - 1) * rho0_ / h_ * Sd;
        C_low_.head(N_ - 1) = 1 / Sp.tail(N_ - 1) * rho0_ / h_ * Sd;

        D_.head(N_) = -2 * rho0_ / (dt_ * intermediary_.head(N_) * Sp);
        E_.head(N_).setZero();
        F_.head(N_) = -rho0_ * Gl / (h_ * Sp);
        G_.head(N_) = rho0_ * Gg / (h_ * Sp);

        A_rad_.head(N_).setConstant(1 / dt_);
        B_rad_.head(N_).setConstant(1 / dt_);
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::Process(ftype inputFlow)
{
    // These views should not have any cost
    auto rho_now = rho_now_.head(N_);
    auto rho_next = rho_next_.head(N_);
    auto vel = vel_.head(N_ - 1);

    auto A = A_.head(N_);
    auto B = B_.head(N_);
    auto C_top = C_top_.head(N_ - 1);
    auto C_low = C_low_.head(N_ - 1);
    auto D = D_.head(N_);
    auto E = E_.head(N_);
    auto F = F_.head(N_);
    auto Gv = G_.head(N_);
    auto A_rad = A_rad_.head(N_);
    auto B_rad = B_rad_.head(N_);

    auto dv = d_plus_v_.head(N_);
    auto dSp = d_S_primal_.head(N_);
    auto Sp = S_primal_.head(N_);

    auto wm = wall_mass_.head(N_);
    auto ws = wall_stiffness_.head(N_);
    auto wdisp = wall_displacement_.head(N_);
    auto wvn = wall_vel_now_.head(N_);
    auto wvx = wall_vel_next_.head(N_);

    if (yielding_walls) {
        dv.setZero();
        dv.head(N_ - 1) += C_top * vel;
        dv.tail(N_ - 1) += C_low * vel;

        rho_next = (1 / A) * (B * rho_now + dv + D * wvn + E * wdisp + F * radiation_flow + Gv * inputFlow -
                              rho0_ * (1 / Sp * dSp));

        vel = vel - c0_ * c0_ / rho0_ / h_ * dt_ * (rho_next.tail(N_ - 1) - rho_next.head(N_ - 1));

        wvx =
            (1 / A_rad) * (B_rad * wvn - (ws / wm) * wdisp + (Sp * c0_ * c0_ / wm * ftype(0.5)) * (rho_now + rho_next));

        wdisp += dt_ * ftype(0.5) * (wvn + wvx);
        radiation_flow += dt_ * c0_ * c0_ / L_rad_ * ftype(0.5) * (rho_next(N_ - 1) + rho_now(N_ - 1));

    } else {
        dv.setZero();
        dv.head(N_ - 1) += C_top * vel;
        dv.tail(N_ - 1) += C_low * vel;

        rho_next = (1 / A) * (B * rho_now + dv + F * radiation_flow + Gv * inputFlow - rho0_ * (1 / Sp * dSp));

        vel = vel - c0_ * c0_ / rho0_ / h_ * dt_ * (rho_next.tail(N_ - 1) - rho_next.head(N_ - 1));

        radiation_flow += dt_ * c0_ * c0_ / L_rad_ * ftype(0.5) * (rho_next(N_ - 1) + rho_now(N_ - 1));
    }

    rho_now_.head(N_) = rho_next_.head(N_);
    wall_vel_now_.head(N_) = wall_vel_next_.head(N_);

    S_direct_last_.head(N_ + 1) = S_direct_.head(N_ + 1);
    S_primal_last_.head(N_) = S_primal_.head(N_);

    if (time_varying_geometry_) {
        for (int i = 0; i < N_ + 1 && i < N_lpf_; ++i) {
            S_direct_[i] = static_cast<ftype>(lp_filters_[i].Process(static_cast<double>(S_target_[i])));
        }

        ComputeDiscreteGreometry();
        UpdateWallParameters();
        UpdateRadiationParameters();
        UpdateCoefficients();

        d_S_primal_.head(N_) = (S_primal_.head(N_) - S_primal_last_.head(N_)) / dt_;
    }
}

template<typename ftype>
std::tuple<ftype, ftype> WebsterFDTD<ftype>::GetIOLinearDependencyCoefficients()
{
    return {c0_ * c0_ * ftype(0.5) *
                (rho_now_(0) + 1 / A_(0) *
                                   (B_(0) * rho_now_(0) - rho0_ / h_ * S_dual_(0) / S_primal_(0) * vel_(0) -
                                    rho0_ * (1 / S_primal_(0) * d_S_primal_(0)))),
            ftype(0.5) * c0_ * c0_ * (1 / A_(0)) * G_(0)};
}

template<typename ftype>
void WebsterFDTD<ftype>::set_N_lpf(int num)
{
    N_lpf_ = std::clamp(num, 1, kMaxN + 1);
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i] = Biquad(sr_, kLowPass, 10.0f, 0.0f, 0.7f);
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::set_lp_frequency(int index, ftype freq)
{
    if (index >= 0 && index < N_lpf_) {
        lp_filters_[index].set_freq(freq);
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::set_lp_Q(int index, ftype Q)
{
    if (index >= 0 && index < N_lpf_) {
        lp_filters_[index].set_Q(Q);
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::set_lp_frequencies(ftype freq)
{
    lpf_frequency_ = std::clamp(freq, ftype(0.1), ftype(100.0));
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i].set_freq(lpf_frequency_);
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::set_lp_Qs(ftype Q)
{
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i].set_Q(Q);
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::filterSdirectTarget()
{
    for (int i = 0; i < N_ + 1 && i < N_lpf_; ++i) {
        S_target_[i] = static_cast<ftype>(lp_filters_[i].Process(static_cast<double>(S_target_[i])));
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::initializeLPFStates()
{
    for (int i = 0; i < N_lpf_ && i < N_ + 1; ++i) {
        lp_filters_[i].InitializeState(static_cast<double>(S_target_[i]));
    }
}

template class WebsterFDTD<float>;
template class WebsterFDTD<double>;

} // namespace tarte
