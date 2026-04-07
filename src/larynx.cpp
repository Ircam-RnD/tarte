#include "larynx.h"

#include <iostream>
namespace tarte {
template<typename ftype>
Larynx<ftype>::Larynx(ftype samplerate, bool yielding_walls)
{
    sr_ = samplerate;
    dt_ = 1 / sr_;

    resonator_ = std::make_shared<WebsterFDTD<ftype>>(sr_);
    resonator_->set_yielding_walls(yielding_walls);

    A0_inv_ = dt_;
    mass_matrix_inv_ = masses_.inverse();

    // clang-format off
  elongation_matrix_ << 1, 0, -1,
                      0, 1, -1,
                      0, 0, 1,
                      1, -1, 0;
    // clang-format on
    stiffness_matrix_ = elongation_matrix_.transpose() * stiffnesses_ * elongation_matrix_;

    dissipation_coefficients_.diagonal()(0) = 2 * xi_ * sqrt(stiffnesses_.diagonal()(0) * masses_.diagonal()(0));
    dissipation_coefficients_.diagonal()(1) = 2 * xi_ * sqrt(stiffnesses_.diagonal()(1) * masses_.diagonal()(1));
    dissipation_coefficients_.diagonal()(2) = xi_ * sqrt(stiffnesses_.diagonal()(2) * masses_.diagonal()(2));
    dissipation_coefficients_.diagonal()(3) = 0;
    dissipation_matrix_ = elongation_matrix_.transpose() * dissipation_coefficients_ * elongation_matrix_;

    p_.setZero();
    q_.setZero();
    r_.setZero();
    Psub_centered_.setZero();
    Psub_.setZero();

    kinetic_energy_.setZero();
    potential_energy_.setZero();
};

template<typename ftype>
void Larynx<ftype>::FillMassesInterpenetrationsAndAreas()
{
    masses_interpenetrations_ = q_(idx_next_, Eigen::all).transpose() - rest_positions_;

    areas_below_masses_ = widths_.cwiseProduct(softplusMatrix(-masses_interpenetrations_,
                                                              kEpsilonSmooth_)); // 2 Comes from symmetric configuration
    smoothed_is_opened_ =
        (-(masses_interpenetrations_ / kEpsilonSmooth_).array().tanh().matrix() + Eigen::Vector<ftype, 3>::Ones()) / 2;
    masses_interpenetrations_derivatives_ = softplusDerivativeMatrix(masses_interpenetrations_, kEpsilonSmooth_);
    masses_interpenetrations_ = softplusMatrix(masses_interpenetrations_, kEpsilonSmooth_);
}

template<typename ftype>
void Larynx<ftype>::ComputeEffectiveAreas()
{
    area_ratio_ = areas_below_masses_(1) / (areas_below_masses_(0) + 1e-14);
    effective_surfaces_Psub_.setZero();
    effective_surfaces_Psup_.setZero();
    if (area_ratio_ > 1) {
        effective_surfaces_Psup_ = widths_.cwiseProduct(lengths_).cwiseProduct(smoothed_is_opened_);
    } else {
        effective_surfaces_Psub_(0) =
            widths_(0) * lengths_(0) * (1 - area_ratio_ * area_ratio_) * (smoothed_is_opened_(0));
        effective_surfaces_Psup_(0) = widths_(0) * lengths_(0) * (area_ratio_ * area_ratio_) * (smoothed_is_opened_(0));
        effective_surfaces_Psup_(1) = widths_(1) * lengths_(1) * (smoothed_is_opened_(1));
    }
}

template<typename ftype>
void Larynx<ftype>::ComputeNonlienarDissipationVector()
{
    area_min_ = areas_below_masses_(Eigen::seq(0, 1)).minCoeff();
    mean_flow_ = area_min_ * sqrt(2 / (kt_ * rho0_) * abs(Psub_(idx_now_) - Psup_)) * sgn(Psub_(idx_now_) - Psup_);
    Rk_ = mean_flow_ / (Psub_(idx_now_) - Psup_ + std::copysign(1e-14, Psub_(idx_now_) - Psup_));
}

template<typename ftype>
void Larynx<ftype>::ComputeSavVector()
{
    elongations_ = elongation_matrix_ * q_(idx_next_, Eigen::all).transpose();

    Enl_ = 0.25 * eta_stiffness_ *
           (stiffnesses_.diagonal().array() * elongations_.array() * elongations_.array() * elongations_.array() *
            elongations_.array())
               .sum();
    Enl_ += 0.5 * (contact_stiffness_ * (masses_interpenetrations_(0) * masses_interpenetrations_(0) +
                                         masses_interpenetrations_(1) * masses_interpenetrations_(1))) +
            contact_stiffness_ * eta_contact_stiffness_ *
                (pow(masses_interpenetrations_(0), alpha_contact_stiffness_ + 1) +
                 pow(masses_interpenetrations_(1), alpha_contact_stiffness_ + 1)) /
                (alpha_contact_stiffness_ + 1); // Contact

    // Fnl_ should be modified to include the derivative of the smoothing function
    // here
    Fnl_ =
        eta_stiffness_ * elongation_matrix_.transpose() *
        (stiffnesses_.diagonal().array() * elongations_.array() * elongations_.array() * elongations_.array()).matrix();
    Fnl_(0) += contact_stiffness_ *
               (eta_contact_stiffness_ * pow(masses_interpenetrations_(0), alpha_contact_stiffness_) +
                masses_interpenetrations_(0)) *
               masses_interpenetrations_derivatives_(0);
    Fnl_(1) += contact_stiffness_ *
               (eta_contact_stiffness_ * pow(masses_interpenetrations_(1), alpha_contact_stiffness_) +
                masses_interpenetrations_(1)) *
               masses_interpenetrations_derivatives_(1);

    g_sav_ = Fnl_ / (sqrt(2 * Enl_) + 1e-14);

    if (control_term_) {
        elongations_ = elongation_matrix_ * 0.5 * (q_(idx_next_, Eigen::all) + q_(idx_now_, Eigen::all)).transpose();

        masses_interpenetrations_ =
            softplusMatrix(0.5 * (q_(idx_next_, Eigen::all) + q_(idx_now_, Eigen::all)).transpose() - rest_positions_,
                           kEpsilonSmooth_);

        Enl_ = 0.25 * eta_stiffness_ *
               (stiffnesses_.diagonal().array() * elongations_.array() * elongations_.array() * elongations_.array() *
                elongations_.array())
                   .sum();
        Enl_ +=
            0.5 * (contact_stiffness_ * (masses_interpenetrations_(0) * masses_interpenetrations_(0) +
                                         masses_interpenetrations_(1) * masses_interpenetrations_(1))) +
            contact_stiffness_ * eta_contact_stiffness_ *
                (pow(masses_interpenetrations_(0), alpha_contact_stiffness_ + 1) + pow(masses_interpenetrations_(1),
                                                                                       alpha_contact_stiffness_ + 1)) /
                (alpha_contact_stiffness_ + 1); // Contact

        epsilon_sav_ = r_(idx_now_) - sqrt(2 * Enl_);
        g_sav_ =
            g_sav_.transpose() - (lambda_sav_ * epsilon_sav_ * dt_ *
                                  (p_(idx_now_, Eigen::all).array() > 0)
                                      .select(Eigen::Array<ftype, 1, 3>::Ones(), -Eigen::Array<ftype, 1, 3>::Ones()) /
                                  (p_(idx_now_, Eigen::all).template lpNorm<1>() + 1e-14))
                                     .matrix();
    }
}

template<typename ftype>
void Larynx<ftype>::Process(ftype Pin)
{
    // Optional computations needed for power balance variables
    sub_glottal_flow =
        -Rk_ * (Psub_(idx_next_) - Psup_) + 0.5 * effective_surfaces_Psub_.transpose() * mass_matrix_inv_ *
                                                (p_(idx_now_, Eigen::all) + p_(idx_next_, Eigen::all)).transpose();
    dissipated_power_flow_ = Rk_ * pow(Psub_(idx_next_) - Psup_, 2);
    auto pmid = 0.5 * (p_(idx_now_, Eigen::all) + p_(idx_next_, Eigen::all)).transpose();
    dissipated_power_folds_ = pmid.transpose() * mass_matrix_inv_ * dissipation_matrix_ * mass_matrix_inv_ * pmid;
    dissipated_power_ = dissipated_power_flow_ + dissipated_power_folds_;

    external_power_sub_ = sub_glottal_flow * Psub_(idx_next_);
    external_power_sup_ = sup_glottal_flow * Psup_;
    external_power_ = external_power_sub_ + external_power_sup_;

    // Step 0: Advance state
    idx_now_ = idx_next_;
    idx_next_ = (idx_now_ + 1) % 2;

    Psub_centered_(idx_next_) = Pin;                                                 // P^{n+1}
    Psub_(idx_next_) = (Psub_centered_(idx_next_) + Psub_centered_(idx_now_)) * 0.5; // P^{n+1/2}
    // Step 1: q_ update
    q_(idx_next_, Eigen::all) =
        q_(idx_now_, Eigen::all).transpose() + dt_ * mass_matrix_inv_ * p_(idx_now_, Eigen::all).transpose();

    // Step 2: Rk_ and g_sav_ explicit computation
    FillMassesInterpenetrationsAndAreas();
    ComputeEffectiveAreas();
    ComputeNonlienarDissipationVector();
    ComputeSavVector();

    // Step 3: get resonator_ feedback coefficients
    std::tie(a_resonator_, b_resonator_) = resonator_->GetIOLinearDependencyCoefficients();

    C0_feedback_ =
        1 / (Rk_ + 1 / b_resonator_) *
        (a_resonator_ / b_resonator_ + Rk_ * Psub_(idx_next_) +
         0.5 * effective_surfaces_Psup_.transpose() * mass_matrix_inv_ * p_(idx_now_, Eigen::all).transpose());
    C1_feedback_ = 1 / (Rk_ + 1 / b_resonator_) * 0.5 * mass_matrix_inv_ * effective_surfaces_Psup_;

    // // Step 4: solve for pnext using Woodburry
    rhs_ = -stiffness_matrix_ * q_(idx_next_, Eigen::all).transpose() +
           (1 / dt_ * Eigen::Matrix<ftype, 3, 3>::Identity() - dissipation_matrix_ * mass_matrix_inv_) *
               p_(idx_now_, Eigen::all).transpose() -
           dt_ / 4 * g_sav_ * (g_sav_.transpose() * mass_matrix_inv_ * p_(idx_now_, Eigen::all).transpose()) -
           g_sav_ * r_[idx_now_] - effective_surfaces_Psub_ * Psub_(idx_next_) -
           effective_surfaces_Psup_ * C0_feedback_;

    u_woodburry_.col(0) = dt_ / 4 * g_sav_;
    u_woodburry_.col(1) = effective_surfaces_Psup_;

    v_woodburry_.row(0) = mass_matrix_inv_ * g_sav_;
    v_woodburry_.row(1) = C1_feedback_;
    woodburry_inv_ = (Eigen::Matrix<ftype, 2, 2>::Identity() + v_woodburry_ * A0_inv_ * u_woodburry_).inverse();

    p_(idx_next_, Eigen::all) =
        A0_inv_ * rhs_ - A0_inv_ * A0_inv_ * u_woodburry_ * woodburry_inv_ * (v_woodburry_ * rhs_);

    // Step 5: r_ update
    r_(idx_next_) = r_(idx_now_) + 0.5 * dt_ * g_sav_.transpose() * mass_matrix_inv_ *
                                       (p_(idx_now_, Eigen::all) + p_(idx_next_, Eigen::all)).transpose();

    // Step 6: Resonator update
    Psup_ = C0_feedback_ + C1_feedback_.transpose() * p_(idx_next_, Eigen::all).transpose();
    sup_glottal_flow = (Psup_ - a_resonator_) / b_resonator_;
    // sup_glottal_flow
    //     = Rk_ * (Psub_(idx_next_) - Psup_)
    //       + 0.5 * effective_surfaces_Psup_.transpose() * mass_matrix_inv_
    //             * (p_(idx_now_, Eigen::all) + p_(idx_next_, Eigen::all)).transpose();
    resonator_->Process(sup_glottal_flow);

    // Optional computations needed for power balance variables
    auto qmid = 0.5 * (q_(idx_now_, Eigen::all) + q_(idx_next_, Eigen::all)).transpose();
    kinetic_energy_(idx_next_) =
        0.5 * p_(idx_now_, Eigen::all) *
        (Eigen::Matrix<ftype, 3, 3>::Identity() - dt_ * dt_ / 4 * mass_matrix_inv_ * stiffness_matrix_ -
         dt_ / 2 * mass_matrix_inv_ * dissipation_matrix_) *
        mass_matrix_inv_ * p_(idx_now_, Eigen::all).transpose();
    potential_energy_(idx_next_) =
        0.5 * qmid.transpose() * stiffness_matrix_ * qmid + 0.5 * r_(idx_now_) * r_(idx_now_);

    stored_power_kinetic_ = (kinetic_energy_(idx_next_) - kinetic_energy_(idx_now_)) / dt_;
    stored_power_potential_ = (potential_energy_(idx_next_) - potential_energy_(idx_now_)) / dt_;
    stored_power_ = stored_power_kinetic_ + stored_power_potential_;
};

template class Larynx<float>;
template class Larynx<double>;

} // namespace tarte