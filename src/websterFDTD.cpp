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
    S_dual_.head(N_ - 1) = 0.5 * (S_direct_.head(N_ - 1) + S_direct_.segment(1, N_ - 1));
    S_primal_(0) = S_direct_(0);
    // The following is to ensure stability while keeping direct control of the
    // geometry of the primal grid at the boundaries
    if (S_primal_(0) < 0.5 * S_dual_(0)) {
        S_dual_(0) = 2 * S_primal_(0);
    }
    S_primal_(N_ - 1) = S_direct_(N_ - 1);
    if (S_primal_(N_ - 1) < 0.5 * S_dual_(N_ - 2)) {
        S_dual_(N_ - 2) = 2 * S_primal_(N_ - 1);
    }
    S_primal_.segment(1, N_ - 2) = 0.5 * (S_dual_.head(N_ - 2) + S_dual_.segment(1, N_ - 2));
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
    auto Sp = S_primal_.head(N_);
    auto Sd = S_dual_.head(N_ - 1);
    auto gamma = gamma_primal_.head(N_);
    auto mw = wall_area_mass_;
    auto rw = wall_area_damping_;
    auto kw = wall_area_stiffness_;

    ftype rhoc2 = rho0_ * c02_;

    if (yielding_walls) {
        intermediary_.head(N_) = 2 / dt_ + rw / mw + kw / (2 * mw) * dt_;

        A_.head(N_) = 1 / dt_ + gamma * rhoc2 / (2 * mw * Sp * intermediary_.head(N_));
        A_(N_ - 1) += rhoc2 / (Sp(N_ - 1) * 2 * h_ * R_rad_) + rhoc2 * dt_ / (Sp(N_ - 1) * 4 * h_ * L_rad_);

        D_.head(N_) = -2 * rho0_ * gamma / (dt_ * mw * intermediary_.head(N_) * Sp);
        E_.head(N_) = rho0_ * gamma * kw / (intermediary_.head(N_) * Sp * mw);

        A_rad_.head(N_) = intermediary_.head(N_) * 0.5;
        B_rad_.head(N_) = 2 / dt_ - A_rad_.head(N_);

    } else {
        A_.head(N_) = 1 / dt_;
        A_(N_ - 1) += rhoc2 / (Sp(N_ - 1) * 2 * h_ * R_rad_) + rhoc2 * dt_ / (Sp(N_ - 1) * 4 * h_ * L_rad_);
    }

    B_.head(N_) = 2 / dt_ - A_.head(N_);
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
                (rho_now_ac()(0) +
                 (1 / A_(0)) * (B_(0) * rho_now_ac()(0) - S_dual_(0) / S_primal_(0) * rho0_ / h_ * vel_(0) +
                                D_(0) * wall_momentum_now_ac()(0) + E_(0) * wall_displacement_(0) -
                                rho0_ * (d_S_primal_(0) / S_primal_(0)))),
            ftype(0.5) * c02_ * (1 / A_(0)) * G_};
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::ComputePowers()
{
    kinetic_energy_fluid_[!flip_] = 0;
    potential_energy_fluid_[!flip_] = 0;
    kinetic_energy_walls_[!flip_] = 0;
    potential_energy_walls_[!flip_] = 0;
    kinetic_energy_radiation_[!flip_] = 0;

    P_stored_fluid_ = (kinetic_energy_fluid_[!flip_] - kinetic_energy_fluid_[flip_] + potential_energy_fluid_[!flip_] -
                       potential_energy_fluid_[flip_]) /
                      dt_;
    P_stored_walls_ = (kinetic_energy_walls_[!flip_] - kinetic_energy_walls_[flip_] + potential_energy_walls_[!flip_] -
                       potential_energy_walls_[flip_]) /
                      dt_;
    P_stored_radiation_ = (kinetic_energy_radiation_[!flip_] - kinetic_energy_radiation_[flip_]) / dt_;
    P_stored_tot_ = P_stored_fluid_ + P_stored_walls_ + P_stored_radiation_;

    P_diss_walls_ = 0;
    P_diss_radiation_ = 0;
    P_diss_tot_ = P_diss_walls_ + P_diss_radiation_;
    P_in_ = 0;
    P_tot_ = P_in_ + P_diss_tot_ + P_stored_tot_;
}

template<typename ftype, int kMaxN>
void WebsterFDTD<ftype, kMaxN>::BuildLaplaceStateSpace(Eigen::MatrixXd& matinternal,
                                                       Eigen::VectorXd& Sxu,
                                                       Eigen::RowVectorXd& matoutZ,
                                                       Eigen::RowVectorXd& matoutTFFlow,
                                                       Eigen::RowVectorXd& matoutTFPressure) const
{
    using Cplx = std::complex<double>;

    const int Nx = N_;     // # pressure (Sp) states
    const int Nv = N_ - 1; // # velocity (Sd) states

    Eigen::ArrayXd Sp = S_primal_.head(Nx).template cast<double>();
    Eigen::ArrayXd Sd = S_dual_.head(Nv).template cast<double>();
    const double h = static_cast<double>(h_);
    const double rho0 = static_cast<double>(rho0_);
    const double c0 = static_cast<double>(c0_);
    const double Lrad = static_cast<double>(L_rad_);
    const double Rrad = static_cast<double>(R_rad_);

    // Finite-difference matrix (same as Dplus/Dmin in the Python version)
    Eigen::MatrixXd Dplus = Eigen::MatrixXd::Zero(Nx, Nv);
    for (int i = 0; i < Nv; ++i) {
        Dplus(i, i) += 1.0 / h;
        Dplus(i + 1, i) += -1.0 / h;
    }
    Eigen::MatrixXd Dmin = -Dplus.transpose();                       // Nv x Nx
    Eigen::MatrixXd S0 = -Dmin * Sp.inverse().matrix().asDiagonal(); // Nv x Nx

    Eigen::ArrayXd Gl = Eigen::ArrayXd::Zero(Nx);
    Gl(Nx - 1) = 1.0; // outlet
    Eigen::ArrayXd Gg = Eigen::ArrayXd::Zero(Nx);
    Gg(0) = 1.0; // glottis

    auto fillCommonTF = [&](Eigen::RowVectorXd& mTF, int total) {
        mTF = Eigen::RowVectorXd::Zero(total);
        mTF(Nv + Nx - 1) = rho0 / (Rrad * h * Sp(Nx - 1)); // resistive branch
        mTF(total - 1) = 1.0 / Lrad;                       // inductive branch
    };

    if (yielding_walls) {
        // Layout: [Nv | Nx | Nx | Nx | 1]
        const int total = Nv + 3 * Nx + 1;
        Eigen::ArrayXd perim = gamma_primal_.head(Nx).template cast<double>();
        const double mw = static_cast<double>(wall_area_mass_);
        const double kw = static_cast<double>(wall_area_stiffness_);
        const double bw = static_cast<double>(wall_area_damping_);

        Eigen::ArrayXd Hdiag(total);
        Hdiag.segment(0, Nv) = h * Sd * rho0;
        Hdiag.segment(Nv, Nx) = h * Sp * (c0 * c0) / rho0;
        Hdiag.segment(Nv + Nx, Nx) = h * perim / mw;
        Hdiag.segment(Nv + 2 * Nx, Nx) = h * perim * kw;
        Hdiag(total - 1) = Lrad;

        Eigen::ArrayXd SRrad = -(rho0 * rho0) / (h * Rrad) * (Gl / Sp.square());
        Eigen::ArrayXd SLrad = -rho0 / Lrad * (Gl / Sp);
        Eigen::ArrayXd SRw = -bw / perim;
        Eigen::ArrayXd S0w = -1.0 / perim;
        Eigen::ArrayXd Sfw = -rho0 / Sp;

        Eigen::MatrixXd Sxx = Eigen::MatrixXd::Zero(total, total);
        Sxx.block(0, Nv, Nv, Nx) = S0;
        Sxx.block(Nv, 0, Nx, Nv) = -S0.transpose();
        Sxx.block(Nv, Nv, Nx, Nx) = SRrad.matrix().asDiagonal();
        Sxx.block(Nv, Nv + Nx, Nx, Nx) = Sfw.matrix().asDiagonal();
        Sxx.block(Nv, Nv + 3 * Nx, Nx, 1) = SLrad.matrix();
        Sxx.block(Nv + Nx, Nv, Nx, Nx) = (-Sfw.matrix()).asDiagonal();
        Sxx.block(Nv + Nx, Nv + Nx, Nx, Nx) = SRw.matrix().asDiagonal();
        Sxx.block(Nv + Nx, Nv + 2 * Nx, Nx, Nx) = S0w.matrix().asDiagonal();
        Sxx.block(Nv + 2 * Nx, Nv + Nx, Nx, Nx) = (-S0w.matrix()).asDiagonal();
        Sxx.block(total - 1, Nv, 1, Nx) = -SLrad.matrix().transpose();
        Sxx /= h;

        Eigen::ArrayXd SxuArr = Eigen::ArrayXd::Zero(total);
        SxuArr.segment(Nv, Nx) = rho0 / Sp * Gg;
        SxuArr /= h;

        matinternal = (Sxx * Hdiag.matrix().asDiagonal());
        Sxu = SxuArr.matrix();
        matoutZ = (SxuArr.matrix().transpose() * Hdiag.matrix().asDiagonal());

        Eigen::RowVectorXd mTF;
        fillCommonTF(mTF, total);
        matoutTFFlow = (mTF * Hdiag.matrix().asDiagonal());

        matoutTFPressure.resize(total);
        matoutTFPressure.setZero();
        matoutTFPressure(Nv + Nx - 1) = rho0_ / (Sp[Nx - 1] * h) * Hdiag(Nv + Nx - 1);
    } else {
        // Layout: [Nv | Nx | 1]
        const int total = Nv + Nx + 1;

        Eigen::ArrayXd Hdiag(total);
        Hdiag.segment(0, Nv) = h * Sd * rho0;
        Hdiag.segment(Nv, Nx) = h * Sp * (c0 * c0) / rho0;
        Hdiag(total - 1) = Lrad;

        Eigen::ArrayXd S1 = -(rho0 * rho0) / (h * Rrad) * (Gl / Sp.square());
        Eigen::ArrayXd S2 = -rho0 / Lrad * (Gl / Sp);

        Eigen::MatrixXd Sxx = Eigen::MatrixXd::Zero(total, total);
        Sxx.block(0, Nv, Nv, Nx) = S0;
        Sxx.block(Nv, 0, Nx, Nv) = -S0.transpose();
        Sxx.block(Nv, Nv, Nx, Nx) = S1.matrix().asDiagonal();
        Sxx.block(Nv, Nv + Nx, Nx, 1) = S2.matrix();
        Sxx.block(total - 1, Nv, 1, Nx) = -S2.matrix().transpose();
        Sxx /= h;

        Eigen::ArrayXd SxuArr = Eigen::ArrayXd::Zero(total);
        SxuArr.segment(Nv, Nx) = rho0 / Sp * Gg;
        SxuArr /= h;

        matinternal = (Sxx * Hdiag.matrix().asDiagonal());
        Sxu = SxuArr.matrix();
        matoutZ = (SxuArr.matrix().transpose() * Hdiag.matrix().asDiagonal());

        Eigen::RowVectorXd mTF;
        fillCommonTF(mTF, total);
        matoutTFFlow = (mTF * Hdiag.matrix().asDiagonal());

        matoutTFPressure.resize(total);
        matoutTFPressure.setZero();
        matoutTFPressure(Nv + Nx - 1) = rho0_ / (Sp[Nx - 1] * h) * Hdiag(Nv + Nx - 1);
    }
}

template<typename ftype, int kMaxN>
WebsterFDTD<ftype, kMaxN>::PolesResidues WebsterFDTD<ftype, kMaxN>::ComputePolesResidues()
{
    Eigen::MatrixXd matinternal;
    Eigen::VectorXd Sxu;
    Eigen::RowVectorXd matoutZ, matoutTFFlow, matOutTFPressure;
    BuildLaplaceStateSpace(matinternal, Sxu, matoutZ, matoutTFFlow, matOutTFPressure);

    Eigen::EigenSolver<Eigen::MatrixXd> es(matinternal);
    Eigen::VectorXcd poles = es.eigenvalues();
    Eigen::MatrixXcd V = es.eigenvectors();
    Eigen::MatrixXcd Vinv = V.inverse();

    Eigen::VectorXcd impedance_residues = (matoutZ * V).transpose().array() * (Vinv * Sxu).array();
    Eigen::VectorXcd tranfer_function_flow_residues = (matoutTFFlow * V).transpose().array() * (Vinv * Sxu).array();
    Eigen::VectorXcd tranfer_function_pressure_residues =
        (matOutTFPressure * V).transpose().array() * (Vinv * Sxu).array();

    PolesResidues out;
    out.poles.resize(poles.size());
    Eigen::VectorXcd::Map(&out.poles[0], poles.size()) = poles;

    out.impedanceResidues.resize(impedance_residues.size());
    Eigen::VectorXcd::Map(&out.impedanceResidues[0], impedance_residues.size()) = impedance_residues;

    out.tranferFunctionFlowResidues.resize(tranfer_function_flow_residues.size());
    Eigen::VectorXcd::Map(&out.tranferFunctionFlowResidues[0], tranfer_function_flow_residues.size()) =
        tranfer_function_flow_residues;

    out.tranferFunctionPressureResidues.resize(tranfer_function_pressure_residues.size());
    Eigen::VectorXcd::Map(&out.tranferFunctionPressureResidues[0], tranfer_function_pressure_residues.size()) =
        tranfer_function_pressure_residues;
    return out;
}

template<typename ftype, int kMaxN>
std::vector<typename WebsterFDTD<ftype, kMaxN>::FrequencyResponse> WebsterFDTD<ftype, kMaxN>::ComputeFrequencyResponse(
    const std::vector<double>& frequenciesHz)
{
    PolesResidues poles_residues = ComputePolesResidues();

    const int n_poles = static_cast<int>(poles_residues.poles.size());

    std::vector<FrequencyResponse> out;
    out.reserve(frequenciesHz.size());

    for (double freq: frequenciesHz) {
        const std::complex<double> s(0.0, 2.0 * M_PI * freq);

        std::complex<double> Z(0.0, 0.0);
        std::complex<double> TFflow(0.0, 0.0);
        std::complex<double> TFpressure(0.0, 0.0);

        for (int i = 0; i < n_poles; ++i) {
            const std::complex<double> denom = s - poles_residues.poles.at(i);
            Z += poles_residues.impedanceResidues.at(i) / denom;
            TFflow += poles_residues.tranferFunctionFlowResidues.at(i) / denom;
            TFpressure += poles_residues.tranferFunctionPressureResidues.at(i) / denom;
        }

        FrequencyResponse resp;
        resp.impedance = Z;
        resp.transferFunctionFlow = TFflow;
        resp.transferFunctionPressure = TFpressure;
        out.push_back(resp);
    }

    return out;
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
