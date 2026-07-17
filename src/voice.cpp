#include "voice.h"

#include <iostream>
namespace tarte {
template<VFPairModel vf_pair, typename ftype>
Voice<vf_pair, ftype>::Voice(ftype samplerate, bool yielding_walls)
{
    sr_ = samplerate;
    dt_ = 1 / sr_;

    vfs_ = std::make_shared<vf_pair>();

    resonator_ = std::make_shared<WebsterFDTD<ftype>>(sr_);
    resonator_->set_yielding_walls(yielding_walls);
    resonator_->set_time_varying_geometry(true);

    A0_inv_ = dt_;

    p_.setZero();
    q_.setZero();
    r_.setZero();
    Psub_centered_.setZero();
    Psub_.setZero();

    g_sav_.setZero();
    Fnl_.setZero();

    kinetic_energy_.setZero();
    potential_energy_.setZero();
};

template<VFPairModel vf_pair, typename ftype>
void Voice<vf_pair, ftype>::DspSetup(ftype samplerate, Articulation* art)
{
    sr_ = samplerate;
    dt_ = 1 / sr_;
    resonator_->DspSetup(samplerate, art);

    A0_inv_ = dt_;

    p_.setZero();
    q_.setZero();
    r_.setZero();
    Psub_centered_.setZero();
    Psub_.setZero();
    Psup_ = 0;

    a_resonator_ = b_resonator_ = C0_feedback_ = 0;
    sub_glottal_flow = sup_glottal_flow = 0;

    kinetic_energy_.setZero();
    potential_energy_.setZero();
}

template<VFPairModel vf_pair, typename ftype>
void Voice<vf_pair, ftype>::ComputeNonlinearDissipationCoefficient()
{
    mean_flow_ = vfs_->ComputeFlow(Psub_(idx_now_), Psup_);
    Rk_ = mean_flow_ / (Psub_(idx_now_) - Psup_ + std::copysign(1e-14, Psub_(idx_now_) - Psup_));
}

template<VFPairModel vf_pair, typename ftype>
void Voice<vf_pair, ftype>::ComputeSavVector()
{
    Enl_ = vfs_->Enl(q_(idx_next_, Eigen::placeholders::all));
    vfs_->Fnl(q_(idx_next_, Eigen::placeholders::all), Fnl_);
    g_sav_ = Fnl_ / (sqrt(2 * Enl_) + 1e-14);

    if (control_term_) {
        qmid_ = 0.5 * (q_(idx_next_, Eigen::placeholders::all) + q_(idx_now_, Eigen::placeholders::all));
        Enl_ = vfs_->Enl(qmid_, true);

        epsilon_sav_ = r_(idx_now_) - sqrt(2 * Enl_);
        g_sav_ = g_sav_ - (lambda_sav_ * epsilon_sav_ *
                           (p_(idx_now_, Eigen::placeholders::all).array() > 0)
                               .select(Eigen::Array<ftype, 1, vf_pair::get_N()>::Ones(),
                                       -Eigen::Array<ftype, 1, vf_pair::get_N()>::Ones()) /
                           ((vfs_->MinvOp(p_(idx_now_, Eigen::placeholders::all))).template lpNorm<1>() + 1e-14))
                              .matrix()
                              .transpose();
    }
}

template<VFPairModel vf_pair, typename ftype>
void Voice<vf_pair, ftype>::Process(ftype Pin)
{
    // auto mass_matrix_inv_left = mass_matrix_inv_.diagonal().head(3).asDiagonal();
    // auto mass_matrix_inv_right = mass_matrix_inv_.diagonal().tail(3).asDiagonal();

    // if (compute_powers_) {
    //     // Optional computations needed for power balance variables
    //     pmid_ = 0.5 * (p_(idx_now_, Eigen::placeholders::all) + p_(idx_next_, Eigen::placeholders::all)).transpose();
    //     sub_glottal_flow = -Rk_ * (Psub_(idx_next_) - Psup_) +
    //                        effective_surfaces_Psub_left_.transpose() * mass_matrix_inv_left * pmid_.head(3) +
    //                        effective_surfaces_Psub_right_.transpose() * mass_matrix_inv_right * pmid_.tail(3);
    //     dissipated_power_flow_ = Rk_ * pow(Psub_(idx_next_) - Psup_, 2);
    //     dissipated_power_folds_ = pmid_.head(3).transpose() * mass_matrix_inv_left * dissipation_matrix_left_ *
    //                               mass_matrix_inv_left * pmid_.head(3);
    //     dissipated_power_folds_ += pmid_.tail(3).transpose() * mass_matrix_inv_right * dissipation_matrix_left_ *
    //                                mass_matrix_inv_right * pmid_.tail(3);
    //     dissipated_power_ = dissipated_power_flow_ + dissipated_power_folds_;

    //     external_power_sub_ = sub_glottal_flow * Psub_(idx_next_);
    //     external_power_sup_ = sup_glottal_flow * Psup_;
    //     external_power_ = external_power_sub_ + external_power_sup_;
    // }
    // Step 0: Advance state
    idx_now_ = idx_next_;
    idx_next_ = (idx_now_ + 1) % 2;

    Psub_centered_(idx_next_) = Pin;                                                 // P^{n+1}
    Psub_(idx_next_) = (Psub_centered_(idx_next_) + Psub_centered_(idx_now_)) * 0.5; // P^{n+1/2}
    // Step 1: q_ update
    Minv_p_ = vfs_->MinvOp(p_(idx_now_, Eigen::placeholders::all));
    q_(idx_next_, Eigen::placeholders::all) = q_(idx_now_, Eigen::placeholders::all).transpose() + dt_ * Minv_p_;

    // Step 2: Rk_ and g_sav_ explicit computation
    vfs_->FillIntermediary(q_(idx_next_, Eigen::placeholders::all));
    vfs_->EffectiveAreas(effective_surfaces_Psub_, effective_surfaces_Psup_);
    ComputeNonlinearDissipationCoefficient();
    ComputeSavVector();

    // Step 3: get resonator_ feedback coefficients
    std::tie(a_resonator_, b_resonator_) = resonator_->GetIOLinearDependencyCoefficients();

    C0_feedback_ = 1 / (Rk_ + 1 / b_resonator_) *
                   (a_resonator_ / b_resonator_ + Rk_ * Psub_(idx_next_) +
                    0.5 * effective_surfaces_Psub_.dot(vfs_->MinvOp(p_(idx_now_, Eigen::placeholders::all))));
    C1_feedback_ = 1 / (Rk_ + 1 / b_resonator_) * 0.5 * vfs_->MinvOp(effective_surfaces_Psup_);

    // Step 4: solve for pnext using Woodburry
    vfs_->KOp(q_(idx_next_, Eigen::placeholders::all), KOp_result_);
    vfs_->ROp(Minv_p_, ROp_result_);
    rhs_ = -KOp_result_ + 1 / dt_ * p_(idx_now_, Eigen::placeholders::all).transpose() - ROp_result_ -
           dt_ / 4 * g_sav_ * (g_sav_.transpose() * Minv_p_) - g_sav_ * r_[idx_now_] -
           effective_surfaces_Psub_ * Psub_(idx_next_) - effective_surfaces_Psup_ * C0_feedback_;

    u_woodburry_.col(0) = dt_ / 4 * g_sav_;
    u_woodburry_.col(1) = effective_surfaces_Psup_;

    v_woodburry_.row(0) = vfs_->MinvOp(g_sav_);
    v_woodburry_.row(1) = C1_feedback_.transpose();
    woodburry_inv_ = (Eigen::Matrix<ftype, 2, 2>::Identity() + v_woodburry_ * A0_inv_ * u_woodburry_).inverse();

    p_(idx_next_, Eigen::placeholders::all) =
        A0_inv_ * rhs_ - A0_inv_ * A0_inv_ * u_woodburry_ * woodburry_inv_ * (v_woodburry_ * rhs_);

    // Step 5: r_ update
    r_(idx_next_) =
        r_(idx_now_) +
        0.5 * dt_ *
            g_sav_.dot(vfs_->MinvOp(p_(idx_now_, Eigen::placeholders::all) + p_(idx_next_, Eigen::placeholders::all)));

    // Step 6: Resonator update
    Psup_ = C0_feedback_ + C1_feedback_.transpose() * p_(idx_next_, Eigen::placeholders::all).transpose();
    sup_glottal_flow = (Psup_ - a_resonator_) / b_resonator_;
    // sup_glottal_flow
    //     = Rk_ * (Psub_(idx_next_) - Psup_)
    //       + 0.5 * effective_surfaces_Psup_.transpose() * mass_matrix_inv_
    //             * (p_(idx_now_, Eigen::placeholders::all) + p_(idx_next_, Eigen::placeholders::all)).transpose();
    resonator_->Process(sup_glottal_flow);

    // if (compute_powers_) {
    //     // Optional computations needed for power balance variables
    //     auto qmid =
    //         0.5 * (q_(idx_now_, Eigen::placeholders::all) + q_(idx_next_, Eigen::placeholders::all)).transpose();
    //     kinetic_energy_(idx_next_) =
    //         0.5 * p_(idx_now_, Eigen::seq(0, 2)) *
    //         (Eigen::Matrix<ftype, 3, 3>::Identity() - dt_ * dt_ / 4 * mass_matrix_inv_left * stiffness_matrix_left_ -
    //          dt_ / 2 * mass_matrix_inv_left * dissipation_matrix_left_) *
    //         mass_matrix_inv_left * p_(idx_now_, Eigen::seq(0, 2)).transpose();
    //     kinetic_energy_(idx_next_) +=
    //         0.5 * p_(idx_now_, Eigen::seq(3, 5)) *
    //         (Eigen::Matrix<ftype, 3, 3>::Identity() - dt_ * dt_ / 4 * mass_matrix_inv_right * stiffness_matrix_right_
    //         -
    //          dt_ / 2 * mass_matrix_inv_right * dissipation_matrix_right_) *
    //         mass_matrix_inv_right * p_(idx_now_, Eigen::seq(3, 5)).transpose();

    //     potential_energy_(idx_next_) =
    //         0.5 * qmid.head(3).transpose() * stiffness_matrix_left_ * qmid.head(3) + 0.5 * r_(idx_now_) *
    //         r_(idx_now_);
    //     potential_energy_(idx_next_) +=
    //         0.5 * qmid.tail(3).transpose() * stiffness_matrix_right_ * qmid.tail(3) + 0.5 * r_(idx_now_) *
    //         r_(idx_now_);

    //     stored_power_kinetic_ = (kinetic_energy_(idx_next_) - kinetic_energy_(idx_now_)) / dt_;
    //     stored_power_potential_ = (potential_energy_(idx_next_) - potential_energy_(idx_now_)) / dt_;
    //     stored_power_ = stored_power_kinetic_ + stored_power_potential_;
    // }
};

// template<VFPairModel vf_pair, typename ftype>
// std::tuple<Eigen::Vector<ftype, 3>, Eigen::Vector<ftype, 3>, Eigen::Matrix<ftype, 3, 3>> Voice<
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
template class Voice<BodyCoverPair<float>, float>;
template class Voice<BodyCoverPair<double>, double>;

} // namespace tarte