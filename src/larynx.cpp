#include "larynx.h"

#include <iostream>
namespace tarte {
template<typename ftype>
Larynx<ftype>::Larynx(ftype samplerate, bool yielding_walls)
{
    sr_ = samplerate;
    dt_ = 1 / sr_;

    left_vf_ = std::make_shared<BodyCoverVF<ftype>>();
    right_vf_ = std::make_shared<BodyCoverVF<ftype>>();

    resonator_ = std::make_shared<WebsterFDTD<ftype>>(sr_);
    resonator_->set_yielding_walls(yielding_walls);
    resonator_->set_time_varying_geometry(true);

    A0_inv_ = dt_;
    RecomputeMatrices();

    p_.setZero();
    q_.setZero();
    r_.setZero();
    Psub_centered_.setZero();
    Psub_.setZero();

    kinetic_energy_.setZero();
    potential_energy_.setZero();
};

template<typename ftype>
void Larynx<ftype>::RecomputeMatrices(bool update_from_muscles)
{
    if (update_from_muscles) {
        left_vf_->ComputeParametersFromMuscularActivity();
        right_vf_->ComputeParametersFromMuscularActivity();
    }

    mass_matrix_inv_.diagonal().head(3) = left_vf_->masses().inverse().diagonal();
    mass_matrix_inv_.diagonal().tail(3) = right_vf_->masses().inverse().diagonal();

    stiffness_matrix_left_ = BodyCoverVF<ftype>::elongation_matrix_.transpose() * left_vf_->stiffnesses() *
                             BodyCoverVF<ftype>::elongation_matrix_;
    stiffness_matrix_right_ = BodyCoverVF<ftype>::elongation_matrix_.transpose() * right_vf_->stiffnesses() *
                              BodyCoverVF<ftype>::elongation_matrix_;

    dissipation_matrix_left_ = BodyCoverVF<ftype>::elongation_matrix_.transpose() *
                               left_vf_->dissipation_coefficients() * BodyCoverVF<ftype>::elongation_matrix_;
    dissipation_matrix_right_ = BodyCoverVF<ftype>::elongation_matrix_.transpose() *
                                left_vf_->dissipation_coefficients() * BodyCoverVF<ftype>::elongation_matrix_;
}

template<typename ftype>
void Larynx<ftype>::DspSetup(ftype samplerate)
{
    sr_ = samplerate;
    dt_ = 1 / sr_;
    resonator_->DspSetup(samplerate);

    A0_inv_ = dt_;
}

template<typename ftype>
void Larynx<ftype>::FillMassesInterpenetrationsAndAreas()
{
    masses_interpenetrations_ = (q_(idx_next_, Eigen::seq(0, 2)).transpose() - left_vf_->rest_positions() +
                                 q_(idx_next_, Eigen::seq(3, 5)).transpose() - right_vf_->rest_positions());

    areas_below_masses_ = 0.5 * (left_vf_->lengths() + right_vf_->lengths())
                                    .cwiseProduct(softplusMatrix(-masses_interpenetrations_, kEpsilonSmooth_));
    smoothed_is_opened_ =
        (-(masses_interpenetrations_ / kEpsilonSmooth_).array().tanh().matrix() + Eigen::Vector<ftype, 3>::Ones()) / 2;
    masses_interpenetrations_derivatives_ = softplusDerivativeMatrix(masses_interpenetrations_, kEpsilonSmooth_);
    masses_interpenetrations_ = softplusMatrix(masses_interpenetrations_, kEpsilonSmooth_);
}

template<typename ftype>
void Larynx<ftype>::ComputeEffectiveAreas()
{
    area_ratio_ = areas_below_masses_(1) / (areas_below_masses_(0) + 1e-14);
    effective_surfaces_Psub_left_.setZero();
    effective_surfaces_Psup_left_.setZero();
    effective_surfaces_Psub_right_.setZero();
    effective_surfaces_Psup_right_.setZero();
    if (area_ratio_ > 1) {
        effective_surfaces_Psup_left_ =
            left_vf_->lengths().cwiseProduct(left_vf_->thicknesses()).cwiseProduct(smoothed_is_opened_);
        effective_surfaces_Psup_right_ =
            right_vf_->lengths().cwiseProduct(right_vf_->thicknesses()).cwiseProduct(smoothed_is_opened_);
    } else {
        effective_surfaces_Psub_left_(0) = left_vf_->lengths()(0) * left_vf_->thicknesses()(0) *
                                           (1 - area_ratio_ * area_ratio_) * (smoothed_is_opened_(0));
        effective_surfaces_Psup_left_(0) = left_vf_->lengths()(0) * left_vf_->thicknesses()(0) *
                                           (area_ratio_ * area_ratio_) * (smoothed_is_opened_(0));
        effective_surfaces_Psup_left_(1) =
            left_vf_->lengths()(1) * left_vf_->thicknesses()(1) * (smoothed_is_opened_(1));

        effective_surfaces_Psub_right_(0) = right_vf_->lengths()(0) * right_vf_->thicknesses()(0) *
                                            (1 - area_ratio_ * area_ratio_) * (smoothed_is_opened_(0));
        effective_surfaces_Psup_right_(0) = right_vf_->lengths()(0) * right_vf_->thicknesses()(0) *
                                            (area_ratio_ * area_ratio_) * (smoothed_is_opened_(0));
        effective_surfaces_Psup_right_(1) =
            right_vf_->lengths()(1) * right_vf_->thicknesses()(1) * (smoothed_is_opened_(1));
    }
}

template<typename ftype>
void Larynx<ftype>::ComputeNonlinearDissipationVector()
{
    area_min_ = areas_below_masses_(Eigen::seq(0, 1)).minCoeff();
    mean_flow_ = area_min_ * sqrt(2 / (kt_ * rho0_) * abs(Psub_(idx_now_) - Psup_)) * sgn(Psub_(idx_now_) - Psup_);
    Rk_ = mean_flow_ / (Psub_(idx_now_) - Psup_ + std::copysign(1e-14, Psub_(idx_now_) - Psup_));
}

template<typename ftype>
void Larynx<ftype>::ComputeSavVector()
{
    Enl_ = 0;
    Fnl_.setZero();

    // Left fold
    elongations_ = BodyCoverVF<ftype>::elongation_matrix_ * q_(idx_next_, Eigen::seq(0, 2)).transpose();

    Enl_ += 0.25 * left_vf_->eta_stiffness() *
            (left_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
             elongations_.array() * elongations_.array())
                .sum();
    Fnl_.head(3) += left_vf_->eta_stiffness() * BodyCoverVF<ftype>::elongation_matrix_.transpose() *
                    (left_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
                     elongations_.array())
                        .matrix();
    // Right fold
    elongations_ = BodyCoverVF<ftype>::elongation_matrix_ * q_(idx_next_, Eigen::seq(3, 5)).transpose();

    Enl_ += 0.25 * right_vf_->eta_stiffness() *
            (right_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
             elongations_.array() * elongations_.array())
                .sum();
    Fnl_.tail(3) += right_vf_->eta_stiffness() * BodyCoverVF<ftype>::elongation_matrix_.transpose() *
                    (right_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
                     elongations_.array())
                        .matrix();

    // Contact
    Enl_ += 0.5 * (contact_stiffness_ * (masses_interpenetrations_(0) * masses_interpenetrations_(0) +
                                         masses_interpenetrations_(1) * masses_interpenetrations_(1))) +
            contact_stiffness_ * eta_contact_stiffness_ *
                (pow(masses_interpenetrations_(0), alpha_contact_stiffness_ + 1) +
                 pow(masses_interpenetrations_(1), alpha_contact_stiffness_ + 1)) /
                (alpha_contact_stiffness_ + 1); // Contact

    Fnl_(0) += contact_stiffness_ *
               (eta_contact_stiffness_ * pow(masses_interpenetrations_(0), alpha_contact_stiffness_) +
                masses_interpenetrations_(0)) *
               masses_interpenetrations_derivatives_(0);
    Fnl_(1) += contact_stiffness_ *
               (eta_contact_stiffness_ * pow(masses_interpenetrations_(1), alpha_contact_stiffness_) +
                masses_interpenetrations_(1)) *
               masses_interpenetrations_derivatives_(1);
    Fnl_(3) += contact_stiffness_ *
               (eta_contact_stiffness_ * pow(masses_interpenetrations_(0), alpha_contact_stiffness_) +
                masses_interpenetrations_(0)) *
               masses_interpenetrations_derivatives_(0);
    Fnl_(4) += contact_stiffness_ *
               (eta_contact_stiffness_ * pow(masses_interpenetrations_(1), alpha_contact_stiffness_) +
                masses_interpenetrations_(1)) *
               masses_interpenetrations_derivatives_(1);

    g_sav_ = Fnl_ / (sqrt(2 * Enl_) + 1e-14);

    if (control_term_) {
        Enl_ = 0;

        elongations_ = BodyCoverVF<ftype>::elongation_matrix_ * 0.5 *
                       (q_(idx_next_, Eigen::seq(0, 2)) + q_(idx_now_, Eigen::seq(0, 2))).transpose();

        Enl_ += 0.25 * left_vf_->eta_stiffness() *
                (left_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
                 elongations_.array() * elongations_.array())
                    .sum();

        elongations_ = BodyCoverVF<ftype>::elongation_matrix_ * 0.5 *
                       (q_(idx_next_, Eigen::seq(3, 5)) + q_(idx_now_, Eigen::seq(3, 5))).transpose();

        Enl_ += 0.25 * right_vf_->eta_stiffness() *
                (right_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
                 elongations_.array() * elongations_.array())
                    .sum();

        masses_interpenetrations_ =
            softplusMatrix(0.5 * ((q_(idx_next_, Eigen::seq(0, 2)) + q_(idx_now_, Eigen::seq(0, 2))).transpose() -
                                  left_vf_->rest_positions() +
                                  (q_(idx_next_, Eigen::seq(3, 5)) + q_(idx_now_, Eigen::seq(3, 5))).transpose() +
                                  -right_vf_->rest_positions()),
                           kEpsilonSmooth_);

        Enl_ +=
            0.5 * (contact_stiffness_ * (masses_interpenetrations_(0) * masses_interpenetrations_(0) +
                                         masses_interpenetrations_(1) * masses_interpenetrations_(1))) +
            contact_stiffness_ * eta_contact_stiffness_ *
                (pow(masses_interpenetrations_(0), alpha_contact_stiffness_ + 1) + pow(masses_interpenetrations_(1),
                                                                                       alpha_contact_stiffness_ + 1)) /
                (alpha_contact_stiffness_ + 1); // Contact

        epsilon_sav_ = r_(idx_now_) - sqrt(2 * Enl_);
        g_sav_ = g_sav_ - (lambda_sav_ * epsilon_sav_ * dt_ *
                           (p_(idx_now_, Eigen::all).array() > 0)
                               .select(Eigen::Array<ftype, 1, 6>::Ones(), -Eigen::Array<ftype, 1, 6>::Ones()) /
                           (p_(idx_now_, Eigen::all).template lpNorm<1>() + 1e-14))
                              .matrix()
                              .transpose();
    }
}

template<typename ftype>
void Larynx<ftype>::Process(ftype Pin)
{
    auto mass_matrix_inv_left = mass_matrix_inv_.diagonal().head(3).asDiagonal();
    auto mass_matrix_inv_right = mass_matrix_inv_.diagonal().tail(3).asDiagonal();

    if (compute_powers_) {
        // Optional computations needed for power balance variables
        pmid = 0.5 * (p_(idx_now_, Eigen::all) + p_(idx_next_, Eigen::all)).transpose();
        sub_glottal_flow = -Rk_ * (Psub_(idx_next_) - Psup_) +
                           effective_surfaces_Psub_left_.transpose() * mass_matrix_inv_left * pmid.head(3) +
                           effective_surfaces_Psub_right_.transpose() * mass_matrix_inv_right * pmid.tail(3);
        dissipated_power_flow_ = Rk_ * pow(Psub_(idx_next_) - Psup_, 2);
        dissipated_power_folds_ = pmid.head(3).transpose() * mass_matrix_inv_left * dissipation_matrix_left_ *
                                  mass_matrix_inv_left * pmid.head(3);
        dissipated_power_folds_ = pmid.tail(3).transpose() * mass_matrix_inv_right * dissipation_matrix_left_ *
                                  mass_matrix_inv_right * pmid.tail(3);
        dissipated_power_ = dissipated_power_flow_ + dissipated_power_folds_;

        external_power_sub_ = sub_glottal_flow * Psub_(idx_next_);
        external_power_sup_ = sup_glottal_flow * Psup_;
        external_power_ = external_power_sub_ + external_power_sup_;
    }
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
    ComputeNonlinearDissipationVector();
    ComputeSavVector();

    // Step 3: get resonator_ feedback coefficients
    std::tie(a_resonator_, b_resonator_) = resonator_->GetIOLinearDependencyCoefficients();

    C0_feedback_ = 1 / (Rk_ + 1 / b_resonator_) *
                   (a_resonator_ / b_resonator_ + Rk_ * Psub_(idx_next_) +
                    0.5 * effective_surfaces_Psub_left_.transpose() * mass_matrix_inv_left *
                        p_(idx_now_, Eigen::seq(0, 2)).transpose() +
                    0.5 * effective_surfaces_Psub_right_.transpose() * mass_matrix_inv_right *
                        p_(idx_now_, Eigen::seq(3, 5)).transpose());
    C1_feedback_.head(3) = 1 / (Rk_ + 1 / b_resonator_) * 0.5 * mass_matrix_inv_left * effective_surfaces_Psup_left_;
    C1_feedback_.tail(3) = 1 / (Rk_ + 1 / b_resonator_) * 0.5 * mass_matrix_inv_right * effective_surfaces_Psup_right_;

    // // Step 4: solve for pnext using Woodburry
    rhs_.head(3) =
        -stiffness_matrix_left_ * q_(idx_next_, Eigen::seq(0, 2)).transpose() +
        (1 / dt_ * Eigen::Matrix<ftype, 3, 3>::Identity() - dissipation_matrix_left_ * mass_matrix_inv_left) *
            p_(idx_now_, Eigen::seq(0, 2)).transpose() -
        dt_ / 4 * g_sav_.head(3) * (g_sav_.transpose() * mass_matrix_inv_ * p_(idx_now_, Eigen::all).transpose()) -
        g_sav_.head(3) * r_[idx_now_] - effective_surfaces_Psub_left_ * Psub_(idx_next_) -
        effective_surfaces_Psup_left_ * C0_feedback_;
    rhs_.tail(3) =
        -stiffness_matrix_right_ * q_(idx_next_, Eigen::seq(3, 5)).transpose() +
        (1 / dt_ * Eigen::Matrix<ftype, 3, 3>::Identity() - dissipation_matrix_right_ * mass_matrix_inv_right) *
            p_(idx_now_, Eigen::seq(3, 5)).transpose() -
        dt_ / 4 * g_sav_.tail(3) * (g_sav_.transpose() * mass_matrix_inv_ * p_(idx_now_, Eigen::all).transpose()) -
        g_sav_.tail(3) * r_[idx_now_] - effective_surfaces_Psub_right_ * Psub_(idx_next_) -
        effective_surfaces_Psup_right_ * C0_feedback_;

    u_woodburry_.col(0) = dt_ / 4 * g_sav_;
    u_woodburry_.col(1).head(3) = effective_surfaces_Psup_left_;
    u_woodburry_.col(1).tail(3) = effective_surfaces_Psup_right_;

    v_woodburry_.row(0) = mass_matrix_inv_ * g_sav_;
    v_woodburry_.row(1) = C1_feedback_.transpose();
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

    if (compute_powers_) {
        // Optional computations needed for power balance variables
        auto qmid = 0.5 * (q_(idx_now_, Eigen::all) + q_(idx_next_, Eigen::all)).transpose();
        kinetic_energy_(idx_next_) =
            0.5 * p_(idx_now_, Eigen::seq(0, 2)) *
            (Eigen::Matrix<ftype, 3, 3>::Identity() - dt_ * dt_ / 4 * mass_matrix_inv_left * stiffness_matrix_left_ -
             dt_ / 2 * mass_matrix_inv_left * dissipation_matrix_left_) *
            mass_matrix_inv_left * p_(idx_now_, Eigen::seq(0, 2)).transpose();
        kinetic_energy_(idx_next_) +=
            0.5 * p_(idx_now_, Eigen::seq(3, 5)) *
            (Eigen::Matrix<ftype, 3, 3>::Identity() - dt_ * dt_ / 4 * mass_matrix_inv_right * stiffness_matrix_right_ -
             dt_ / 2 * mass_matrix_inv_right * dissipation_matrix_right_) *
            mass_matrix_inv_left * p_(idx_now_, Eigen::seq(3, 5)).transpose();

        potential_energy_(idx_next_) =
            0.5 * qmid.head(3).transpose() * stiffness_matrix_left_ * qmid.head(3) + 0.5 * r_(idx_now_) * r_(idx_now_);
        potential_energy_(idx_next_) =
            0.5 * qmid.tail(3).transpose() * stiffness_matrix_right_ * qmid.tail(3) + 0.5 * r_(idx_now_) * r_(idx_now_);

        stored_power_kinetic_ = (kinetic_energy_(idx_next_) - kinetic_energy_(idx_now_)) / dt_;
        stored_power_potential_ = (potential_energy_(idx_next_) - potential_energy_(idx_now_)) / dt_;
        stored_power_ = stored_power_kinetic_ + stored_power_potential_;
    }
};

// template<typename ftype>
// std::tuple<Eigen::Vector<ftype, 3>, Eigen::Vector<ftype, 3>, Eigen::Matrix<ftype, 3, 3>> Larynx<
//     ftype>::GetModalCharacteristics()
// {
//     auto mass_matrix_inv_left = mass_matrix_inv_.diagonal().head(3).asDiagonal();
//     auto mass_matrix_inv_right = mass_matrix_inv_.diagonal().tail(3).asDiagonal();
//     // Compute for the left fold
//     // We know that M, K and R are symmetric for this model.
//     // Evaluate eigenfrequencies and eigenvectors from undamped problems
//     Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ftype, 3, 3>> eig(mass_matrix_inv_left * stiffness_matrix_left_);
//     Eigen::Vector<ftype, 3> eigen_frequencies = eig.eigenvalues().cwiseMax(ftype(0)).cwiseSqrt() / (2 * M_PI);

//     Eigen::Matrix<ftype, 3, 3> Phi = eig.eigenvectors();
//     // Normalize wrt max displacement
//     for (int i = 0; i < 3; ++i) {
//         ftype norm = Phi.col(i).cwiseAbs().maxCoeff();
//         if (norm > ftype(1e-12))
//             Phi.col(i) /= norm;
//     }

//     // Modal damping ratio from projected dissipation matrix
//     // ζ_i = φᵢᵀ D φᵢ / (2 ω_i)
//     Eigen::Vector<ftype, 3> zeta = Eigen::Vector<ftype, 3>::Zero();
//     for (int i = 0; i < 3; ++i)
//         if (eigen_frequencies(i) > ftype(1e-12))
//             zeta(i) = Phi.col(i).dot(dissipation_matrix_left_ * Phi.col(i)) /
//                       (ftype(2) * ftype(2) * M_PI * eigen_frequencies(i));

//     return {eigen_frequencies, zeta, Phi};
// }
template class Larynx<float>;
template class Larynx<double>;

} // namespace tarte