#pragma once

#include "noise_generator.h"
#include "utility/smoothing.h"
#include "vocal_folds/body_cover_vf.h"

#include <Eigen/Dense>
#include <tuple>

#include <cmath>

namespace tarte {

enum FoldIdentifier { kRight, kLeft, kBoth };

template<typename ftype>
class BodyCoverPair {
private:
    static const int N{6};
    static const int half_N{3};

    // Both folds
    std::shared_ptr<BodyCoverVF<ftype>> left_vf_, right_vf_;

    // Matrices and intermediary quantities
    Eigen::DiagonalMatrix<ftype, N> mass_matrix_inv_;                                        // Both folds
    Eigen::Matrix<ftype, half_N, half_N> stiffness_matrix_left_, dissipation_matrix_left_;   // Left fold only
    Eigen::Matrix<ftype, half_N, half_N> stiffness_matrix_right_, dissipation_matrix_right_; // Right fold only
    Eigen::Vector<ftype, 4> elongations_;

    // The following quantities need to be defined only once for both folds, as they are related to the fluid flow
    Eigen::Vector<ftype, half_N> masses_interpenetrations_, areas_below_masses_, smoothed_is_opened_,
        masses_interpenetrations_derivatives_;

    // Intermediary variables
    ftype area_ratio_, area_min_;
    ftype mean_flow_;

    /// ----------------------------------- Class specific members ----------------------------------------- ///
    // Fluid parameters
    ftype noise_ratio_{0};
    ftype rho0_{1.2}, c0_{340}, kt_{1.3};
    // Contact parameters
    ftype contact_stiffness_{15}, alpha_contact_stiffness_{1.3};

    // Smoothing at channel closure
    ftype epsilon_smooth_{1e-5};

    // Noise generator
    NoiseGenerator noise_generator_ = NoiseGenerator(NoiseColor::Pink);

    void RecomputeMatrices(bool update_from_muscles)
    {
        if (update_from_muscles) {
            left_vf_->ComputeParametersFromMuscularActivity();
            right_vf_->ComputeParametersFromMuscularActivity();
        } else {
            left_vf_->FillDissipationCoefficients();
            right_vf_->FillDissipationCoefficients();
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
                                    right_vf_->dissipation_coefficients() * BodyCoverVF<ftype>::elongation_matrix_;

        // The left vocal fold is arbitrary taken as reference to update contact stiffness
        contact_stiffness_ = 3 * left_vf_->stiffnesses().diagonal()(1); // Also arbitrary upper stiffness
    }

public:
    // Necessary interface for the solver
    using scalar_type = ftype;
    using state_type = Eigen::Vector<ftype, N>; // Either q or p
    static inline constexpr int get_N() { return N; };
    static inline constexpr int get_half_N() { return half_N; };
    BodyCoverPair()
    {
        left_vf_ = std::make_shared<BodyCoverVF<ftype>>();
        right_vf_ = std::make_shared<BodyCoverVF<ftype>>();

        RecomputeMatrices(true);
    };

    void FillIntermediary(const state_type& state_q)
    {
        left_vf_->rest_positions();
        masses_interpenetrations_ =
            (state_q.head(half_N) - left_vf_->rest_positions() + state_q.tail(half_N) - right_vf_->rest_positions());

        areas_below_masses_ = 0.5 * (left_vf_->lengths() + right_vf_->lengths())
                                        .cwiseProduct(softplusMatrix(-masses_interpenetrations_, epsilon_smooth_));
        smoothed_is_opened_ =
            (-(masses_interpenetrations_ / epsilon_smooth_).array().tanh().matrix() + Eigen::Vector<ftype, 3>::Ones()) /
            2;
        masses_interpenetrations_ =
            (masses_interpenetrations_.array() > 0).select(masses_interpenetrations_, Eigen::Vector<ftype, 3>::Zero());
    };

    void EffectiveAreas(state_type& out_area_P_sub, state_type& out_area_P_sup)
    {
        area_ratio_ = areas_below_masses_(1) / (areas_below_masses_(0) + 1e-14);
        out_area_P_sub.setZero();
        out_area_P_sup.setZero();
        if (area_ratio_ > 1) {
            out_area_P_sup.head(half_N) =
                left_vf_->lengths().cwiseProduct(left_vf_->thicknesses()).cwiseProduct(smoothed_is_opened_);
            out_area_P_sup.tail(half_N) =
                right_vf_->lengths().cwiseProduct(right_vf_->thicknesses()).cwiseProduct(smoothed_is_opened_);

            // The body mass should never be exposed to fluid pressure
            out_area_P_sup(2) = 0;
            out_area_P_sup(half_N + 2) = 0;
        } else {
            out_area_P_sub(0) = left_vf_->lengths()(0) * left_vf_->thicknesses()(0) * (1 - area_ratio_ * area_ratio_) *
                                (smoothed_is_opened_(0));
            out_area_P_sup(0) = left_vf_->lengths()(0) * left_vf_->thicknesses()(0) * (area_ratio_ * area_ratio_) *
                                (smoothed_is_opened_(0));
            out_area_P_sup(1) = left_vf_->lengths()(1) * left_vf_->thicknesses()(1) * (smoothed_is_opened_(1));

            out_area_P_sub(half_N) = right_vf_->lengths()(0) * right_vf_->thicknesses()(0) *
                                     (1 - area_ratio_ * area_ratio_) * (smoothed_is_opened_(0));

            out_area_P_sup(half_N) = right_vf_->lengths()(0) * right_vf_->thicknesses()(0) *
                                     (area_ratio_ * area_ratio_) * (smoothed_is_opened_(0));
            out_area_P_sup(half_N + 1) =
                right_vf_->lengths()(1) * right_vf_->thicknesses()(1) * (smoothed_is_opened_(1));
        }
    };

    ftype ComputeFlow(const ftype& Psub, const ftype& Psup)
    {
        area_min_ = areas_below_masses_(Eigen::seq(0, 1)).minCoeff();
        mean_flow_ = area_min_ * sqrt(2 / (kt_ * rho0_) * abs(Psub - Psup)) * sgn(Psub - Psup);
        return mean_flow_ * (1 + noise_ratio_ * noise_generator_.Process());
        ;
    };

    ftype Enl(const state_type& state_q, bool recompute_intermediary = false)
    {
        if (recompute_intermediary) {
            FillIntermediary(state_q);
        }
        // Left fold
        ftype Enl = 0;
        elongations_ = BodyCoverVF<ftype>::elongation_matrix_ * state_q.head(half_N);
        Enl += 0.25 * left_vf_->eta_stiffness() *
               (left_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
                elongations_.array() * elongations_.array())
                   .sum();
        // Right fold
        elongations_ = BodyCoverVF<ftype>::elongation_matrix_ * state_q.tail(half_N);
        Enl += 0.25 * right_vf_->eta_stiffness() *
               (right_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
                elongations_.array() * elongations_.array())
                   .sum();

        // Contact
        Enl += contact_stiffness_ *
               (pow(masses_interpenetrations_(0), alpha_contact_stiffness_ + 1) +
                pow(masses_interpenetrations_(1), alpha_contact_stiffness_ + 1)) /
               (alpha_contact_stiffness_ + 1);
        return Enl;
    };

    void Fnl(const state_type& state_q, state_type& out, bool recompute_intermediary = false)
    {
        if (recompute_intermediary) {
            FillIntermediary(state_q);
        }
        // Left fold
        elongations_ = BodyCoverVF<ftype>::elongation_matrix_ * state_q.head(half_N);
        out.head(3) = left_vf_->eta_stiffness() * BodyCoverVF<ftype>::elongation_matrix_.transpose() *
                      (left_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
                       elongations_.array())
                          .matrix();
        // Right fold
        elongations_ = BodyCoverVF<ftype>::elongation_matrix_ * state_q.tail(half_N);
        out.tail(3) = right_vf_->eta_stiffness() * BodyCoverVF<ftype>::elongation_matrix_.transpose() *
                      (right_vf_->stiffnesses().diagonal().array() * elongations_.array() * elongations_.array() *
                       elongations_.array())
                          .matrix();

        // Contact
        out(0) += contact_stiffness_ * pow(masses_interpenetrations_(0), alpha_contact_stiffness_);
        out(1) += contact_stiffness_ * pow(masses_interpenetrations_(1), alpha_contact_stiffness_);
        out(3) += contact_stiffness_ * pow(masses_interpenetrations_(0), alpha_contact_stiffness_);
        out(4) += contact_stiffness_ * pow(masses_interpenetrations_(1), alpha_contact_stiffness_);
    };

    void KOp(const state_type& state_q, state_type& out)
    {
        out.template head<half_N>().noalias() = stiffness_matrix_left_ * state_q.template head<half_N>();
        out.template tail<half_N>().noalias() = stiffness_matrix_right_ * state_q.template tail<half_N>();
    };
    void ROp(const state_type& state_p, state_type& out)
    {
        out.template head<half_N>().noalias() = dissipation_matrix_left_ * state_p.template head<half_N>();
        out.template tail<half_N>().noalias() = dissipation_matrix_right_ * state_p.template tail<half_N>();
    };
    state_type MinvOp(const state_type& state_p) { return mass_matrix_inv_ * state_p; };

    /// ----------------------------------- Control functions ----------------------------------------- ///
    void set_rest_positions(const Eigen::Vector<ftype, 3>& rest_positions, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_->set_rest_positions(rest_positions); break;
        case kRight: right_vf_->set_rest_positions(rest_positions); break;
        default:
            left_vf_->set_rest_positions(rest_positions);
            right_vf_->set_rest_positions(rest_positions);
            break;
        }
    };
    void set_rest_positions(const ftype& lower, const ftype& upper, const ftype& body, FoldIdentifier fold_id = kBoth)
    {
        Eigen::Vector<ftype, 3> v;
        v << lower, upper, body;
        set_rest_positions(v, fold_id);
    };

    void set_masses(const Eigen::Vector<ftype, 3>& masses, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_->set_masses(masses); break;
        case kRight: right_vf_->set_masses(masses); break;
        default:
            left_vf_->set_masses(masses);
            right_vf_->set_masses(masses);
            break;
        }
        RecomputeMatrices();
    };
    void set_masses(const ftype& lower, const ftype& upper, const ftype& body, FoldIdentifier fold_id = kBoth)
    {
        Eigen::Vector<ftype, 3> v;
        v << lower, upper, body;
        set_masses(v, fold_id);
    };

    void set_length(const ftype& length, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_->set_length(length); break;
        case kRight: right_vf_->set_length(length); break;
        default:
            left_vf_->set_length(length);
            right_vf_->set_length(length);
            break;
        }
    };

    void set_thicknesses(const Eigen::Vector<ftype, 3>& thicknesses, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_->set_thicknesses(thicknesses); break;
        case kRight: right_vf_->set_thicknesses(thicknesses); break;
        default:
            left_vf_->set_thicknesses(thicknesses);
            right_vf_->set_thicknesses(thicknesses);
            break;
        }
        RecomputeMatrices();
    };
    void set_thicknesses(const ftype& lower, const ftype& upper, const ftype& body, FoldIdentifier fold_id = kBoth)
    {
        Eigen::Vector<ftype, 3> v;
        v << lower, upper, body;
        set_thicknesses(v, fold_id);
    };

    void set_stiffnesses(const Eigen::Vector<ftype, 4>& stiffnesses, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_->set_stiffnesses(stiffnesses); break;
        case kRight: right_vf_->set_stiffnesses(stiffnesses); break;
        default:
            left_vf_->set_stiffnesses(stiffnesses);
            right_vf_->set_stiffnesses(stiffnesses);
            break;
        }
        RecomputeMatrices();
    };
    void set_stiffnesses(const ftype& lower,
                         const ftype& upper,
                         const ftype& body,
                         const ftype& coupling,
                         FoldIdentifier fold_id = kBoth)
    {
        Eigen::Vector<ftype, 4> v;
        v << lower, upper, body, coupling;
        set_stiffnesses(v, fold_id);
    };

    void set_eta_stiffness(const ftype& eta_stiffness, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_->set_eta_stiffness(eta_stiffness); break;
        case kRight: right_vf_->set_eta_stiffness(eta_stiffness); break;
        default:
            left_vf_->set_eta_stiffness(eta_stiffness);
            right_vf_->set_eta_stiffness(eta_stiffness);
            break;
        }
    }
    void set_contact_stiffness(const ftype& contact_stiffness)
    {
        contact_stiffness_ = std::clamp(contact_stiffness, ftype(0), ftype(1e3));
    }
    void set_alpha_contact_stiffness(const ftype& alpha_contact_stiffness)
    {
        alpha_contact_stiffness_ = std::clamp(alpha_contact_stiffness, ftype(1), ftype(5));
    }

    void set_muscles_activation(const ftype& ct_activity,
                                const ftype& ta_activity,
                                const ftype& lc_activity,
                                FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft:
            left_vf_->SetCricothyroidActivity(ct_activity);
            left_vf_->SetThyroarytenoidActivity(ta_activity);
            left_vf_->SetCricoarytenoidActivity(lc_activity);
            break;
        case kRight:
            right_vf_->SetCricothyroidActivity(ct_activity);
            right_vf_->SetThyroarytenoidActivity(ta_activity);
            right_vf_->SetCricoarytenoidActivity(lc_activity);
            break;
        default:
            left_vf_->SetCricothyroidActivity(ct_activity);
            left_vf_->SetThyroarytenoidActivity(ta_activity);
            left_vf_->SetCricoarytenoidActivity(lc_activity);

            right_vf_->SetCricothyroidActivity(ct_activity);
            right_vf_->SetThyroarytenoidActivity(ta_activity);
            right_vf_->SetCricoarytenoidActivity(lc_activity);
            break;
        }
        RecomputeMatrices(true);
    };

    // void set_xi(const ftype xi, FoldIdentifier fold_id = kBoth)
    // {
    //     switch (fold_id) {
    //     case kLeft: left_vf_.xi_ = std::clamp(xi, ftype(0), ftype(1e12)); break;
    //     case kRight: right_vf_.xi_ = std::clamp(xi, ftype(0), ftype(1e12)); break;
    //     default:
    //         left_vf_.xi_ = std::clamp(xi, ftype(0), ftype(1e12));

    //         right_vf_.xi_ = std::clamp(xi, ftype(0), ftype(1e12));
    //         break;
    //     }
    // }

    void set_c0(const ftype& c0) { c0_ = std::clamp(c0, ftype(300), ftype(380)); }
    void set_rho0(const ftype& rho0) { rho0_ = std::clamp(rho0, ftype(1), ftype(15)); }
    void set_kt(const ftype& kt) { kt_ = std::clamp(kt, ftype(1), ftype(2)); }

    void set_epsilon_smooth(const ftype& epsilon_smooth)
    {
        epsilon_smooth_ = std::clamp(epsilon_smooth, ftype(1e-6), ftype(1e-2));
    }
    void set_noise_ratio(const ftype& noise_ratio) { noise_ratio_ = std::clamp(noise_ratio, ftype(0), ftype(1)); }
    void set_gender(Gender gender)
    {
        left_vf_->set_gender(gender);
        right_vf_->set_gender(gender);
    }

    // Getters
    inline Eigen::Vector<ftype, 3> get_rest_positions(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_->rest_positions();
        } else {
            return left_vf_->rest_positions();
        }
    };

    inline Eigen::Vector<ftype, 3> get_masses(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_->masses().diagonal();
        } else {
            return left_vf_->masses().diagonal();
        }
    };

    inline Eigen::Vector<ftype, 3> get_lengths(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_->lengths();
        } else {
            return left_vf_->lengths();
        }
    };
    inline Eigen::Vector<ftype, 3> get_thicknesses(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_->thicknesses();
        } else {
            return left_vf_->thicknesses();
        }
    };
    inline Eigen::Vector<ftype, 4> get_stiffnesses(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_->stiffnesses().diagonal();
        } else {
            return left_vf_->stiffnesses().diagonal();
        }
    };
    inline ftype get_eta_stiffness(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_->eta_stiffness();
        } else {
            return left_vf_->eta_stiffness();
        }
    }

    inline ftype get_contact_stiffness() { return contact_stiffness_; }
    inline ftype get_alpha_contact_stiffness() { return alpha_contact_stiffness_; }
    // inline ftype get_xi(FoldIdentifier fold_id = kBoth)
    // {
    //     if (fold_id == kRight) {
    //         return right_vf_.xi_;
    //     } else {
    //         return left_vf_.xi_;
    //     }
    // }

    inline ftype get_c0() { return c0_; }
    inline ftype get_rho0() { return rho0_; }
    inline ftype get_kt() { return kt_; }

    inline ftype get_epsilon_smooth() { return epsilon_smooth_; }
    inline ftype get_noise_ratio() { return noise_ratio_; }

    inline Gender get_gender() { return left_vf_->get_gender(); }
};

} // namespace tarte