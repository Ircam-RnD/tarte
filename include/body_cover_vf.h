#pragma once

#include "utility/eigen_utility.h"

#include <Eigen/Dense>

#include <cmath>

namespace tarte {
enum Gender { MALE, FEMALE };

template<typename T>
class BodyCoverVF {
    // Model from [1] "Voice simulation with a body‐cover model of the vocal folds", Story and Titze, 1995
    // Muscular control based on paper [2] "Rules for controlling low-dimensional vocal fold models with muscle
    // activation", Titze and Story, 2002

private:
    // Muscle activity, normalized 0-1
    T alpha_ct_{0.5}, alpha_ta_{0.5}, alpha_lc_{0.5}; // cricothyroid, thyroarytenoid, lateral cricoarytenoid

    // Intermediary physical parameters (see section 'II. A minimal parameter set' of the paper [2])
    Gender gender_{Gender::MALE};
    T comp_striffness_, shear_stiffness_;
    T depth_, depth_cover_, depth_body_;
    T shear_mode_nodal_point_;
    T elongation_;
    T thickness_;
    Eigen::Vector<T, 3> thicknesses_{1.5e-3,
                                     1.5e-3,
                                     3e-3}; // effective thicknesses of the 3 masses for use in pressure mappings
    Eigen::Vector<T, 3> lengths_{1e-2, 1e-2, 1e-2};

    // Numerical constants for stress computations (table 1 of [2])
    static constexpr T e1_mucosa_{-0.5}, e2_mucosa_{0.35}, sigma_0_mucosa_{500}, sigma_2_mucosa_{30000}, C_mucosa_{4.4};
    static constexpr T e1_ligament_{-0.5}, e2_ligament_{0.0}, sigma_0_ligament_{400}, sigma_2_ligament_{1393},
        C_ligament_{17};
    static constexpr T e1_ta_{-0.5}, e2_ta_{-0.05}, sigma_0_ta_{1000.0}, sigma_2_ta_{1500.0}, C_ta_{6.5},
        sigma_max_ta_{105000}, e_ta_opt_{0.4}, b_ta_{1.07};

    // Numerical constants for elongation computations
    static constexpr T e_gain_{0.2}, e_torque_ratio_{3}, e_strain_factor{0.2};

    // Other numerical constants
    T density_{1040}, shear_modulus_cover_{500}, shear_modulus_muscle_{1000};

    // Geometry numerical constants
    T base_length_{1.6e-2}; // 1.6e-2 for male, 1e-2 for female
    T depth_mucosa_{0.2e-2}, depth_ligament_{0.2e-2},
        depth_muscle_{0.4e-2}; // 0.2e-2, 0.2e-2, 0.4e-2 for male, 0.15e-2, 0.15e-2, 0.3e-2 for female
    T base_thickness_{3e-3};

    // Stresses
    static T PassiveStress(const T& sigma0, const T& sigma2, const T& e1, const T& e2, const T& C, const T& e)
    {
        if (e < e1) {
            return 0;
        }
        if (e < e2) {
            return -sigma0 / e1 * (e - e1);
        }
        return -sigma0 / e1 * (e - e1) + sigma2 * (exp(C * (e - e2)) - C * (e - e2) - 1);
    }
    static T PassiveStressMucosa(const T& e)
    {
        return PassiveStress(sigma_0_mucosa_, sigma_2_mucosa_, e1_mucosa_, e2_mucosa_, C_mucosa_, e);
    };
    static T PassiveStressLigament(const T& e)
    {
        return PassiveStress(sigma_0_ligament_, sigma_2_ligament_, e1_ligament_, e2_ligament_, C_ligament_, e);
    };
    static T PassiveStressTA(const T& e) { return PassiveStress(sigma_0_ta_, sigma_2_ta_, e1_ta_, e2_ta_, C_ta_, e); };

    T sigma_body_, sigma_cover_;
    T sigma_muscle_, sigma_muscle_passive_, sigma_ligament_, sigma_mucosa_;

    // Lower level parameters: order is lower-upper-body.
    Eigen::DiagonalMatrix<T, 3> masses_{1e-5, 1e-5, 5e-5};
    Eigen::Vector<T, 3> rest_positions_{1.8e-4, 1.79e-4, 3e-3};

    Eigen::DiagonalMatrix<T, 4> stiffnesses_{5, 3.5, 20, 0.5}; // lower, upper, body, lower-upper
    T eta_stiffness_{1e6};

    T xi_{0.4};
    Eigen::DiagonalMatrix<T, 4> dissipation_coefficients_;

public:
    static inline const Eigen::Matrix<T, 4, 3> elongation_matrix_ =
        (Eigen::Matrix<T, 4, 3>() << 1, 0, -1, 0, 1, -1, 0, 0, 1, 1, -1, 0).finished();
    BodyCoverVF() { FillDissipationCoefficients(); };

    // Set muscular activity
    void SetCricothyroidActivity(const T& activity) { alpha_ct_ = std::clamp(activity, T(0), T(1)); }
    void SetThyroarytenoidActivity(const T& activity) { alpha_ta_ = std::clamp(activity, T(0), T(1)); }
    void SetCricoarytenoidActivity(const T& activity) { alpha_lc_ = std::clamp(activity, T(0), T(1)); }
    // Compute low level parameters depending on current muscular activity
    void ComputeParametersFromMuscularActivity()
    {
        // Intermediary values
        // Geometry
        shear_mode_nodal_point_ =
            (1 + alpha_ta_) / 3; // relative version, between 1/3 and 2/3 of the vocal fold thickness_
        elongation_ = e_gain_ * (e_torque_ratio_ * alpha_ct_ - alpha_ta_) - e_strain_factor * alpha_lc_;

        thickness_ = base_thickness_ / (1 + 0.8 * elongation_);
        thicknesses_(0) = thickness_ * shear_mode_nodal_point_;
        thicknesses_(1) = thickness_ * (1 - shear_mode_nodal_point_);
        thicknesses_(2) = thickness_;

        lengths_(0) = lengths_(1) = lengths_(2) = base_length_ * (1 + elongation_);

        depth_body_ = (alpha_ta_ * depth_muscle_ + 0.5 * depth_ligament_) / (1. + 0.2 * elongation_);
        depth_cover_ = (depth_mucosa_ + 0.5 * depth_ligament_) / (1 + 0.2 * elongation_);

        rest_positions_(1) = 0.25 * base_length_ * (1 - 2 * alpha_lc_);
        // rest_positions_(0) = rest_positions_(1) + thickness_ * (0.05 - 0.15 * alpha_ta_); // Convergence rule Eq. 63
        rest_positions_(0) = rest_positions_(1) + thickness_ * tan(0.0001); // "Nearly rectangular" case

        rest_positions_(2) = 3e-3; // Arbitrary, no impact on the dynamics.

        // Stresses
        sigma_ligament_ = PassiveStressLigament(elongation_);
        sigma_mucosa_ = PassiveStressMucosa(elongation_);
        sigma_muscle_passive_ = PassiveStressTA(elongation_);
        sigma_muscle_ = alpha_ta_ * sigma_max_ta_ *
                            std::max(T(0), 1 - b_ta_ * (elongation_ - e_ta_opt_) * (elongation_ - e_ta_opt_)) +
                        sigma_muscle_passive_;
        sigma_body_ = (0.5 * sigma_ligament_ * depth_ligament_ + sigma_muscle_ * depth_muscle_) / depth_body_;
        sigma_cover_ = (0.5 * sigma_ligament_ * depth_ligament_ + sigma_mucosa_ * depth_mucosa_) / depth_cover_;

        // Low level parameters
        masses_.diagonal()(0) = density_ * lengths_(0) * thicknesses_(0) * depth_cover_;
        masses_.diagonal()(1) = density_ * lengths_(1) * thicknesses_(1) * depth_cover_;
        masses_.diagonal()(2) = density_ * lengths_(0) * thickness_ * depth_body_;

        stiffnesses_.diagonal()(0) = 2 * shear_modulus_cover_ * lengths_(0) * thicknesses_(0) / depth_cover_ +
                                     M_PI * M_PI * sigma_cover_ * thicknesses_(0) * (depth_cover_ / lengths_(0));
        stiffnesses_.diagonal()(1) = 2 * shear_modulus_cover_ * lengths_(1) * thicknesses_(1) / depth_cover_ +
                                     M_PI * M_PI * sigma_cover_ * thicknesses_(1) * (depth_cover_ / lengths_(1));
        stiffnesses_.diagonal()(2) = 2 * shear_modulus_muscle_ * lengths_(0) * thickness_ / depth_body_ +
                                     M_PI * M_PI * sigma_body_ * (depth_body_ / lengths_(0)) * thickness_;
        stiffnesses_.diagonal()(3) =
            std::max(T((0.5 * shear_modulus_cover_ * (lengths_(0) * depth_cover_ / thickness_) /
                            (T(1) / T(3) - shear_mode_nodal_point_ * (1 - shear_mode_nodal_point_)) -
                        2 * shear_modulus_cover_ * (lengths_(0) * thickness_ / depth_cover_)) *
                       shear_mode_nodal_point_ * (1 - shear_mode_nodal_point_)),
                     T(0)); // coupling

        FillDissipationCoefficients();
    }

    void FillDissipationCoefficients()
    {
        dissipation_coefficients_.diagonal()(0) = 2 * xi_ * sqrt(stiffnesses_.diagonal()(0) * masses_.diagonal()(0));
        dissipation_coefficients_.diagonal()(1) = 2 * xi_ * sqrt(stiffnesses_.diagonal()(1) * masses_.diagonal()(1));
        dissipation_coefficients_.diagonal()(2) = xi_ * sqrt(stiffnesses_.diagonal()(2) * masses_.diagonal()(2));
        dissipation_coefficients_.diagonal()(3) = 0;
    }

    // Direct accessors to lower level parameters for easier read in the larynx class
    Eigen::DiagonalMatrix<T, 3> masses() { return masses_; }
    Eigen::Vector<T, 3> rest_positions() { return rest_positions_; }
    Eigen::DiagonalMatrix<T, 4> dissipation_coefficients() { return dissipation_coefficients_; }
    Eigen::Vector<T, 3> lengths() { return lengths_; }
    Eigen::Vector<T, 3> thicknesses() { return thicknesses_; }
    Eigen::DiagonalMatrix<T, 4> stiffnesses() { return stiffnesses_; }
    T eta_stiffness() { return eta_stiffness_; }

    // Lower level parameters setters
    void set_masses(const Eigen::Vector<T, 3>& masses) { masses_.diagonal() = ClipEigen(masses.array(), 1e-6, 1e-3); };
    void set_stiffnesses(const Eigen::Vector<T, 4>& stiffnesses)
    {
        stiffnesses_.diagonal() = ClipEigen(stiffnesses.array(), 1e-6, 1e6);
    };
    void set_eta_stiffness(const T& eta_stiffness) { eta_stiffness_ = std::clamp(eta_stiffness, T(0), T(1e12)); };
    void set_rest_positions(const Eigen::Vector<T, 3>& rest_positions)
    {
        rest_positions_ = ClipEigen(rest_positions.array(), T(0), T(1e-2));
    };
    void set_length(const T& length) { lengths_(0) = lengths_(1) = length; };
    void set_thicknesses(const Eigen::Vector<T, 3>& thicknesses)
    {
        thicknesses_ = ClipEigen(thicknesses.array(), 1e-5, 1e3);
    };

    // Set and get gender
    void set_gender(Gender gender)
    {
        gender_ = gender;
        switch (gender_) {
        case MALE:
            base_length_ = 1.6e-2; // 1.6e-2 for male, 1e-2 for female
            depth_mucosa_ = 0.2e-2;
            depth_ligament_ = 0.2e-2;
            depth_muscle_ = 0.4e-2; // 0.2e-2, 0.2e-2, 0.4e-2 for male, 0.15e-2, 0.15e-2, 0.3e-2 for female
            break;
        case FEMALE:
            base_length_ = 1.2e-2; // 1.6e-2 for male, 1e-2 for female
            depth_mucosa_ = 0.15e-2;
            depth_ligament_ = 0152e-2;
            depth_muscle_ = 0.3e-2; // 0.2e-2, 0.2e-2, 0.4e-2 for male, 0.15e-2, 0.15e-2, 0.3e-2 for female
            break;
        }
    }
    Gender get_gender() { return gender_; }
};

} // namespace tarte