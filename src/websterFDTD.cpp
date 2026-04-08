#include "websterFDTD.h"

#include <iostream>

namespace tarte {

template<typename ftype>
WebsterFDTD<ftype>::WebsterFDTD(ftype sampleRate, ftype length)
{
    l0_ = length;
    DspSetup(sampleRate);
}

template<typename ftype>
void WebsterFDTD<ftype>::DspSetup(ftype sampleRate)
{
    sr_ = sampleRate;
    dt_ = 1 / sr_;
    SetNStability();
    x_direct_ = Eigen::Array<ftype, -1, 1>::LinSpaced(N_ + 1, -h_ / 2, l0_ + h_ / 2);
    x_dual_ = x_direct_.segment(1, N_ - 1);
    x_primal_ = Eigen::Array<ftype, -1, 1>::LinSpaced(N_, 0, l0_);

    S_target_.resize(N_ + 1);
    S_target_.setOnes();
    S_direct_.resize(N_ + 1);
    S_dual_.resize(N_ - 1);
    S_primal_.resize(N_);
    S_direct_ = S_target_;
    S_direct_last_ = S_direct_;
    ComputeDiscreteGreometry();
    S_primal_last_ = S_primal_;
    d_S_primal_.resize(N_);
    d_S_primal_.setZero();

    G_lips_.resize(N_);
    G_lips_.setZero();
    G_lips_(N_ - 1) = 1;
    G_glottis_.resize(N_);
    G_glottis_.setZero();
    G_glottis_(0) = 1;

    wall_mass_.resize(N_);
    wall_stiffness_.resize(N_);
    wall_dissipation_.resize(N_);
    UpdateWallParameters();
    UpdateRadiationParameters();

    d_plus_v_.resize(N_);
    d_plus_v_.setZero();
    intermediary_.resize(N_);
    intermediary_.setZero();
    A_.resize(N_);
    A_.setZero();
    B_.resize(N_);
    B_.setZero();
    C_top_.resize(N_ - 1);
    C_top_.setZero();
    C_low_.resize(N_ - 1);
    C_low_.setZero();
    D_.resize(N_);
    D_.setZero();
    E_.resize(N_);
    E_.setZero();
    F_.resize(N_);
    F_.setZero();
    G_.resize(N_);
    G_.setZero();
    A_rad_.resize(N_);
    A_rad_.setZero();
    B_rad_.resize(N_);
    B_rad_.setZero();

    rho_now_.resize(N_);
    rho_now_.setZero();
    rho_next_.resize(N_);
    rho_next_.setZero();
    vel_.resize(N_ - 1);
    vel_.setZero();
    wall_displacement_.resize(N_);
    wall_displacement_.setZero();
    wall_vel_now_.resize(N_);
    wall_vel_now_.setZero();
    wall_vel_next_.resize(N_);
    wall_vel_next_.setZero();
    radiation_flow = 0;

    UpdateCoefficients();

    // Initialize LPF filters
    N_lpf_ = N_ + 1;
    lp_filters_.resize(N_lpf_);
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i] = Biquad(sr_, kLowPass, 10.0f, 0.0f, 0.5f);
        // Initialize filter state to corresponding SdirectTarget value
        if (i < S_target_.size()) {
            lp_filters_[i].InitializeState(static_cast<double>(S_target_[i]));
        }
    }
    set_lp_frequencies(10.0);
    set_lp_Qs(0.5);
}

template<typename ftype>
void WebsterFDTD<ftype>::SetNStability()
{
    h_ = c0_ / sr_;
    N_ = std::floor(l0_ / h_);
    h_ = l0_ / N_;
}

template<typename ftype, typename intype>
void WebsterFDTD<ftype>::SetTargetGeometry(intypt const* in, std::size_t const size)
{
    size_t safe_size = std::min(size, std::size_t(S_target_.size()));
    for (size_t i = 0; i < safe_size; ++i) {
        S_target_[i] = std::max(float(1e-8), float(in[i]));
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::SetTargetGeometryFromArticulation(Articulation articulation)
{
    articulation.getAreas(x_direct_.data(), S_target_.data(), N_ + 1);
    if (!time_varying_geometry_) {
        S_direct_ = S_target_;
        ComputeDiscreteGreometry();
        UpdateCoefficients();
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::SetConstantSection(ftype section)
{
    S_target_ = Eigen::Array<ftype, -1, 1>::Ones(N_ + 1) * section;
    if (!time_varying_geometry_) {
        S_direct_ = S_target_;
        ComputeDiscreteGreometry();
        UpdateRadiationParameters();
        UpdateCoefficients();
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::ComputeDiscreteGreometry()
{
    // Sd = 0.25 * (Sdirect.segment(0, N-1) + 2*Sdirect.segment(1, N-1) +
    // Sdirect.segment(2, N-1));
    S_direct_(0) = S_direct_(1);
    // Sdirect(N-1) = 1e-10;
    S_direct_(N_) = S_direct_(N_ - 1);
    S_dual_ = S_direct_.segment(1, N_ - 1);
    S_primal_ = 0.5 * (S_direct_.segment(0, N_) + S_direct_.segment(1, N_));
    // Sp.segment(1, N-2) = 0.25 * (Sp.segment(0, N-2) + 2*Sp.segment(1, N-2) +
    // Sp.segment(2, N-2));
}

template<typename ftype>
void WebsterFDTD<ftype>::UpdateWallParameters()
{
    wall_mass_ = wall_area_mass_ * 2 * sqrt(S_primal_ * M_PI);
    wall_stiffness_ = wall_mass_ * (2 * M_PI * 70) * (2 * M_PI * 70);
    wall_dissipation_ = wall_area_damping_ * 2 * sqrt(S_primal_ * M_PI);
}

template<typename ftype>
void WebsterFDTD<ftype>::UpdateRadiationParameters()
{
    ftype Aout = S_primal_[N_ - 1];
    lip_radius = sqrt(Aout / M_PI);
    ftype Z0 = rho0_ * c0_ / Aout;
    L_rad_ = 8 * lip_radius / (3 * M_PI * c0_) * Z0;
    R_rad_ = 128 / (9 * M_PI * M_PI) * Z0;
}

template<typename ftype>
void WebsterFDTD<ftype>::UpdateCoefficients()
{
    if (yielding_walls) {
        ftype rhoc2 = rho0_ * c0_ * c0_;
        intermediary_ = 2 / dt_ + wall_dissipation_ / wall_mass_ + wall_stiffness_ / (2 * wall_mass_) * dt_;

        A_ = (1 / dt_ + 4 * M_PI * rhoc2 / (2 * wall_mass_ * intermediary_) +
              rhoc2 * G_lips_ / (S_primal_ * 2 * h_ * R_rad_) + rhoc2 * G_lips_ * dt_ / (S_primal_ * 4 * h_ * L_rad_));
        B_ = (1 / dt_ - 4 * M_PI * rhoc2 / (2 * wall_mass_ * intermediary_) -
              rhoc2 * G_lips_ / (S_primal_ * 2 * h_ * R_rad_) - rhoc2 * G_lips_ * dt_ / (S_primal_ * 4 * h_ * L_rad_));
        C_top_ = -1 / S_primal_.head(N_ - 1) * rho0_ / h_ * S_dual_;
        C_low_ = 1 / S_primal_.tail(N_ - 1) * rho0_ / h_ * S_dual_;

        D_ = -2 * rho0_ / (dt_ * intermediary_ * S_primal_);
        E_ = rho0_ / (intermediary_ * S_primal_) * wall_stiffness_ / (wall_mass_);
        F_ = -rho0_ * G_lips_ / (h_ * S_primal_);
        G_ = rho0_ * G_glottis_ / (h_ * S_primal_);

        A_rad_ = 1 / dt_ + wall_dissipation_ / (2 * wall_mass_) + wall_stiffness_ / (4 * wall_mass_) * dt_;
        B_rad_ = 1 / dt_ - wall_dissipation_ / (2 * wall_mass_) - wall_stiffness_ / (4 * wall_mass_) * dt_;
    } else {
        ftype rhoc2 = rho0_ * c0_ * c0_;
        intermediary_ = 2 / dt_;

        A_ = (1 / dt_ + rhoc2 * G_lips_ / (S_primal_ * 2 * h_ * R_rad_) +
              rhoc2 * G_lips_ * dt_ / (S_primal_ * 4 * h_ * L_rad_));
        B_ = (1 / dt_ - rhoc2 * G_lips_ / (S_primal_ * 2 * h_ * R_rad_) -
              rhoc2 * G_lips_ * dt_ / (S_primal_ * 4 * h_ * L_rad_));
        C_top_ = -1 / S_primal_.head(N_ - 1) * rho0_ / h_ * S_dual_;
        C_low_ = 1 / S_primal_.tail(N_ - 1) * rho0_ / h_ * S_dual_;

        D_ = -2 * rho0_ / (dt_ * intermediary_ * S_primal_);
        E_.setZero();
        F_ = -rho0_ * G_lips_ / (h_ * S_primal_);
        G_ = rho0_ * G_glottis_ / (h_ * S_primal_);

        A_rad_ = 1 / dt_;
        B_rad_ = 1 / dt_;
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::Process(ftype inputFlow)
{
    if (yielding_walls) {
        d_plus_v_.setZero();
        d_plus_v_.head(N_ - 1) = C_top_ * vel_;
        d_plus_v_.tail(N_ - 1) += C_low_ * vel_;
        rho_next_ = (1 / A_) * (B_ * rho_now_ + d_plus_v_ + D_ * wall_vel_now_ + E_ * wall_displacement_ +
                                F_ * radiation_flow + G_ * inputFlow - rho0_ * (1 / S_primal_ * d_S_primal_));
        vel_ = (1 / (1 / dt_)) *
               ((1 / dt_) * vel_ - c0_ * c0_ / (rho0_) / h_ * (rho_next_.tail(N_ - 1) - rho_next_.head(N_ - 1)));
        wall_vel_next_ =
            (1 / A_rad_) * (B_rad_ * wall_vel_now_ - (wall_stiffness_ / wall_mass_) * (wall_displacement_) +
                            (S_primal_ * c0_ * c0_ / wall_mass_ * 0.5) * (rho_now_ + rho_next_));
        wall_displacement_ = wall_displacement_ + dt_ * 0.5 * (wall_vel_now_ + wall_vel_next_);
        radiation_flow = radiation_flow + dt_ * c0_ * c0_ / L_rad_ * 0.5 * (rho_next_(N_ - 1) + rho_now_(N_ - 1));
    } else {
        d_plus_v_.setZero();
        d_plus_v_.head(N_ - 1) = C_top_ * vel_;
        d_plus_v_.tail(N_ - 1) += C_low_ * vel_;
        rho_next_ = (1 / A_) * (B_ * rho_now_ + d_plus_v_ + F_ * radiation_flow + G_ * inputFlow -
                                rho0_ * (1 / S_primal_ * d_S_primal_));
        vel_ = (1 / (1 / dt_)) *
               ((1 / dt_) * vel_ - c0_ * c0_ / (rho0_) / h_ * (rho_next_.tail(N_ - 1) - rho_next_.head(N_ - 1)));
        radiation_flow = radiation_flow + dt_ * c0_ * c0_ / L_rad_ * 0.5 * (rho_next_(N_ - 1) + rho_now_(N_ - 1));
    }
    rho_now_ = rho_next_;
    wall_vel_now_ = wall_vel_next_;

    S_direct_last_ = S_direct_;
    S_primal_last_ = S_primal_;

    if (time_varying_geometry_) {
        // Apply filtering to each element of SdirectTarget
        for (int i = 0; i < S_target_.size() && i < N_lpf_; i++) {
            S_direct_[i] = static_cast<ftype>(lp_filters_[i].Process(static_cast<double>(S_target_[i])));
            // Sdirect[i] = SdirectTarget[i];
        }

        ComputeDiscreteGreometry();
        UpdateWallParameters();
        UpdateRadiationParameters();
        UpdateCoefficients();

        d_S_primal_ = (S_primal_ - S_primal_last_) / dt_;
    }
}

template<typename ftype>
std::tuple<ftype, ftype> WebsterFDTD<ftype>::GetIOLinearDependencyCoefficients()
{
    return {c0_ * c0_ * 0.5 *
                (rho_now_(0) + 1 / A_(0) *
                                   (B_(0) * rho_now_(0) - rho0_ / (h_)*S_dual_(0) / S_primal_(0) * vel_(0) -
                                    rho0_ * (1 / S_primal_(0) * d_S_primal_(0)))),
            0.5 * c0_ * c0_ * (1 / A_(0)) * G_(0)};
}

// LPF filter methods implementation
template<typename ftype>
void WebsterFDTD<ftype>::set_N_lpf(int num)
{
    N_lpf_ = std::max(1, num);
    lp_filters_.resize(N_lpf_);
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
    for (int i = 0; i < N_lpf_; ++i) {
        lp_filters_[i].set_freq(freq);
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
    for (int i = 0; i < S_target_.size() && i < N_lpf_; ++i) {
        S_target_[i] = static_cast<ftype>(lp_filters_[i].Process(static_cast<double>(S_target_[i])));
    }
}

template<typename ftype>
void WebsterFDTD<ftype>::initializeLPFStates()
{
    for (int i = 0; i < N_lpf_ && i < S_target_.size(); ++i) {
        lp_filters_[i].InitializeState(static_cast<double>(S_target_[i]));
    }
}

template class WebsterFDTD<float>;
template class WebsterFDTD<double>;

} // namespace tarte