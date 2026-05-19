#include "single_reed.h"

#include <iostream>

namespace tarte {

template<typename ftype>
SingleReed<ftype>::SingleReed(float samplerate)
{
    sr_ = samplerate;
    dt_ = 1 / sr_;

    resonator_ = std::make_shared<WebsterFDTD<ftype, 300>>(sr_);

    dissipation_coefficient_ = 2 * mass_ * damping_;

    p_.setZero();
    q_.setZero();
    r_.setZero();
    Psub_centered_.setZero();
    Psub_.setZero();

    kinetic_energy_.setZero();
    potential_energy_.setZero();
};

template<typename ftype>
void SingleReed<ftype>::DspSetup(float samplerate)
{
    sr_ = samplerate;
    dt_ = 1 / sr_;

    resonator_->DspSetup(samplerate);
    resonator_->SetConstantSection(M_PI * 1e-4);

    p_.setZero();
    q_.setZero();
    r_.setZero();
    Psub_centered_.setZero();
    Psub_.setZero();

    kinetic_energy_.setZero();
    potential_energy_.setZero();
}

template<typename ftype>
void SingleReed<ftype>::fillOpeningAndInterpenetration()
{
    interpenetration_ = q_(idx_next_) - lay_position_;

    opening_ = width_ * softplus(-interpenetration_, epsilon_smooth_);
    // interpenetration_ left unsmoothed
    // interpenetration_ = softplus(interpenetration_, epsilon_smooth_);
    // interpenetration_derivative_ = softplusDerivative(interpenetration_,
    // epsilon_smooth_);
    smoothed_is_opened_ = 1;

    if (interpenetration_ >= 0) {
        interpenetration_derivative_ = 1;
    } else {
        interpenetration_ = 0;
        interpenetration_derivative_ = 0;
    }
}

template<typename ftype>
void SingleReed<ftype>::computeRk()
{
    // Here the pressure drop is evaluated with a full timestep of lag: i.e. at
    // timestep n-1/2. Some kind of evaluation of Psup_ at timestep n, r_ better
    // timestep n+1/2 nees to be used for better results.
    mean_flow_ = opening_ * sqrt(2 / (kt_ * rho_0_) * abs(Psub_(idx_now_) - Psup_)) *
                 tanh((Psub_(idx_now_) - Psup_) / epsilon_smooth_P_); // Some kind of smoothing of the sign function
    random_value_ = -1 + static_cast<float>(rand()) / static_cast<float>(RAND_MAX / (2));
    // Sign coherent noise for 0 < noise_ratio_ < 1
    noise_flow_ = mean_flow_ * (1 + noise_ratio_ * random_value_);
    R_k_ = noise_flow_ / (Psub_(idx_now_) - Psup_ + std::copysign(1e-10, Psub_(idx_now_) - Psup_));
}

template<typename ftype>
void SingleReed<ftype>::computegSAV()
{
    // Modified SAV evaluation for contact: EXP-2 from
    // Implicit and explicit schemes for energy-stable simulation of string
    // vibrations with collisions: Refinement, analysis, and comparison. 2024 van
    // Walstijn et al.

    // E_nl_ = contact_stiffness_ * pow(interpenetration_, alpha_contact_stiffness_ + 1)
    //       / (alpha_contact_stiffness_ + 1);

    // F_nl_ = contact_stiffness_ * pow(interpenetration_, alpha_contact_stiffness_)
    //       * interpenetration_derivative_;

    // g_sav_ = F_nl_ / (sqrt(2 * E_nl_) + 1e-14);

    if (interpenetration_ > 0) {
        E_nl_ =
            contact_stiffness_ * pow(interpenetration_, alpha_contact_stiffness_ + 1) / (alpha_contact_stiffness_ + 1);

        F_nl_ = contact_stiffness_ * pow(interpenetration_, alpha_contact_stiffness_) * interpenetration_derivative_;

        g_sav_ = std::copysign(1, r_(idx_now_)) * F_nl_ / (sqrt(2 * E_nl_) + 1e-14);
        if (g_sav_ * (r_(idx_now_) + g_sav_ * (q_inter_ - q_(idx_now_)) / 4) < 0) {
            rhs_ = -stiffness_ * q_(idx_next_) + (1 / dt_ - dissipation_coefficient_ / (2 * mass_)) * p_(idx_now_) +
                   surface_ * (Psub_(idx_next_) - C0_feedback_);

            p_inter_ = rhs_ / (1 / dt_ + dissipation_coefficient_ / (2 * mass_) + surface_ * C1_feedback_);
            q_inter_ = q_(idx_next_) + dt_ / mass_ * p_inter_;
            g_sav_ = -2 * r_(idx_now_) /
                     (q_inter_ - q_(idx_now_) +
                      +std::copysign(1e-12,
                                     q_inter_ - q_(idx_now_))); // q_(idx_now_) two steps before q_inter_
                                                                // g_sav_ = 0;
        }
    } else {
        rhs_ = -stiffness_ * q_(idx_next_) + (1 / dt_ - dissipation_coefficient_ / (2 * mass_)) * p_(idx_now_) +
               surface_ * (Psub_(idx_next_) - C0_feedback_);

        p_inter_ = rhs_ / (1 / dt_ + dissipation_coefficient_ / (2 * mass_) + surface_ * C1_feedback_);
        q_inter_ = q_(idx_next_) + dt_ / mass_ * p_inter_;
        g_sav_ = -2 * r_(idx_now_) /
                 (q_inter_ - q_(idx_now_) +
                  std::copysign(1e-12,
                                q_inter_ - q_(idx_now_))); // q_(idx_now_) two steps before q_inter_
    }
}

template<typename ftype>
void SingleReed<ftype>::Process(float Pin)
{
    // Optional computations needed for power balance variables
    if (compute_powers_) {
        mouth_flow_ = resonator_flow_;
        dissipated_power_flow_ = R_k_ * pow(Psub_(idx_next_) - Psup_, 2);
        ftype pmid = 0.5 * (p_(idx_now_) + p_(idx_next_));
        dissipated_power_reed_ = pow(pmid / mass_, 2) * dissipation_coefficient_;
        dissipated_power_ = dissipated_power_flow_ + dissipated_power_reed_;

        external_power_sub_ = -mouth_flow_ * Psub_(idx_next_);
        external_power_sup_ = resonator_flow_ * Psup_;
        external_power_ = external_power_sub_ + external_power_sup_;
    }

    // Step 0: Advance state
    idx_now_ = idx_next_;
    idx_next_ = (idx_now_ + 1) % 2;

    Psub_centered_(idx_next_) = Pin;                                                 // P^{n+1}
    Psub_(idx_next_) = (Psub_centered_(idx_next_) + Psub_centered_(idx_now_)) * 0.5; // P^{n+1/2}
    // Step 1: q_ update
    q_(idx_next_) = q_(idx_now_) + dt_ / mass_ * p_(idx_now_);

    // Step 2 and 3: R_k_ and g_sav_ explicit computation, resonator_ feedbac
    // coefficients
    fillOpeningAndInterpenetration();
    computeRk();
    std::tie(a_resonator_, b_resonator_) = resonator_->GetIOLinearDependencyCoefficients();
    C0_feedback_ = 1 / (R_k_ + 1 / b_resonator_) *
                   (a_resonator_ / b_resonator_ + R_k_ * Psub_(idx_next_) +
                    0.5 * surface_ * smoothed_is_opened_ / mass_ * p_(idx_now_));
    C1_feedback_ = 1 / (R_k_ + 1 / b_resonator_) * 0.5 / mass_ * surface_ * smoothed_is_opened_;

    computegSAV();

    // Step 4: solve for pnext
    rhs_ = -stiffness_ * q_(idx_next_) +
           (1 / dt_ - 0.25 * dt_ * g_sav_ * g_sav_ / mass_ - dissipation_coefficient_ / (2 * mass_)) * p_(idx_now_) +
           surface_ * smoothed_is_opened_ * (Psub_(idx_next_) - C0_feedback_) - g_sav_ * r_(idx_now_);

    p_(idx_next_) = rhs_ / (1 / dt_ + 0.25 * dt_ * g_sav_ * g_sav_ / mass_ + dissipation_coefficient_ / (2 * mass_) +
                            surface_ * smoothed_is_opened_ * C1_feedback_);

    // Step 5: r_ update
    r_(idx_next_) = r_(idx_now_) + 0.5 * dt_ * g_sav_ / mass_ * (p_(idx_now_) + p_(idx_next_));

    // Step 6: Resonator update
    Psup_ = C0_feedback_ + C1_feedback_ * p_(idx_next_);
    resonator_flow_ = (Psup_ - a_resonator_) / b_resonator_;
    // resonator_flow_ = R_k_ * (Psub_(idx_next_) - Psup_)
    //                 + 0.5 * surface_ / mass_ * (p_(idx_now_) + p_(idx_next_));
    resonator_->Process(resonator_flow_);

    // Optional computations needed for power balance variables
    if (compute_powers_) {
        ftype qmid = 0.5 * (q_(idx_now_) + q_(idx_next_));
        kinetic_energy_(idx_next_) =
            0.5 * p_(idx_now_) * (1 - dt_ * dt_ / 4 * stiffness_ / mass_) * p_(idx_now_) / mass_;
        potential_energy_(idx_next_) = 0.5 * qmid * stiffness_ * qmid + 0.5 * r_(idx_now_) * r_(idx_now_);

        stored_power_kinetic_ = (kinetic_energy_(idx_next_) - kinetic_energy_(idx_now_)) / dt_;
        stored_power_potential_ = (potential_energy_(idx_next_) - potential_energy_(idx_now_)) / dt_;
        stored_power_ = stored_power_kinetic_ + stored_power_potential_;
    }
};

template class SingleReed<float>;
template class SingleReed<double>;

} // namespace tarte