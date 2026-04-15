#pragma once

#include <utility/eigen_utility.h>
#include <utility/maths.h>
#include <utility/smoothing.h>
#include <websterFDTD.h>

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <memory>
#include <string_view>
#include <tuple>
#include <vector>

namespace tarte {

template<typename T>
struct BodyCoverVF {
    // Vocal fold: order is lower-upper-body
    Eigen::DiagonalMatrix<T, 3> masses_{1e-5, 1e-5, 5e-5};
    Eigen::Vector<T, 3> rest_positions_{1.8e-4, 1.79e-4, 3e-3};

    Eigen::DiagonalMatrix<T, 4> stiffnesses_{5, 3.5, 20, 0.5}; // lower, upper, body, lower-upper
    T eta_stiffness_{1e6};

    T xi_{0.4};
    Eigen::DiagonalMatrix<T, 4> dissipation_coefficients_;
    static inline const Eigen::Matrix<T, 4, 3> elongation_matrix_ =
        (Eigen::Matrix<T, 4, 3>() << 1, 0, -1, 0, 1, -1, 0, 0, 1, 1, -1, 0).finished();

    BodyCoverVF() { FillDissipationCoefficients(); };
    void FillDissipationCoefficients()
    {
        dissipation_coefficients_.diagonal()(0) = 2 * xi_ * sqrt(stiffnesses_.diagonal()(0) * masses_.diagonal()(0));
        dissipation_coefficients_.diagonal()(1) = 2 * xi_ * sqrt(stiffnesses_.diagonal()(1) * masses_.diagonal()(1));
        dissipation_coefficients_.diagonal()(2) = xi_ * sqrt(stiffnesses_.diagonal()(2) * masses_.diagonal()(2));
        dissipation_coefficients_.diagonal()(3) = 0;
    }
};

enum FoldIdentifier { kRight, kLeft, kBoth };

template<typename ftype>
class Larynx {
private:
    // Folds hysical parameters are in there
    BodyCoverVF<ftype> left_vf_, right_vf_;

    // The lengths and widths of the masses are supposed to be the same for both folds
    Eigen::Vector<ftype, 3> lengths_{1.5e-3, 1.5e-3, 3e-3};
    Eigen::Vector<ftype, 3> widths_{1e-2, 1e-2, 0};

    // Contact parameters are supposed to be the same for both folds
    ftype contact_stiffness_{15}, eta_contact_stiffness_{1e6}, alpha_contact_stiffness_{3};

    ftype rho0_{1.2}, c0_{340}, kt_{1.3};

    // Matrices and intermediary quantities
    Eigen::DiagonalMatrix<ftype, 6> mass_matrix_inv_;                              // Both folds
    Eigen::Matrix<ftype, 3, 3> stiffness_matrix_left_, dissipation_matrix_left_;   // Left fold only
    Eigen::Matrix<ftype, 3, 3> stiffness_matrix_right_, dissipation_matrix_right_; // Right fold only
    Eigen::Vector<ftype, 4> elongations_;

    // The following quantities need to be defined only once for both folds, as they are realted to the fluid flow
    Eigen::Vector<ftype, 3> masses_interpenetrations_, areas_below_masses_, smoothed_is_opened_,
        masses_interpenetrations_derivatives_;
    Eigen::Vector<ftype, 3> effective_surfaces_Psub_, effective_surfaces_Psup_;

    Eigen::Vector<ftype, 6> g_sav_{0, 0, 0, 0, 0, 0}, Fnl_{0, 0, 0, 0, 0, 0};
    ftype Enl_, epsilon_sav_;
    bool const control_term_{true};
    ftype lambda_sav_{1000};

    ftype area_ratio_, area_min_;
    ftype mean_flow_, Rk_;

    ftype a_resonator_, b_resonator_;
    ftype C0_feedback_;
    Eigen::Vector<ftype, 6> C1_feedback_;

    ftype A0_inv_;
    Eigen::Vector<ftype, 6> rhs_;
    Eigen::Matrix<ftype, 2, 6> v_woodburry_;
    Eigen::Matrix<ftype, 6, 2> u_woodburry_;
    Eigen::Matrix<ftype, 2, 2> woodburry_inv_;

    ftype const kEpsilonSmooth_{1e-5};

    // State
    Eigen::Matrix<ftype, 2, 6> p_, q_;
    Eigen::Vector<ftype, 6> pmid;
    Eigen::Vector<ftype, 2> r_; // One SAV variable for both folds
    std::size_t idx_now_{0}, idx_next_{1};

    Eigen::Vector<ftype, 2> Psub_, Psub_centered_;
    ftype Psup_;

    // Power variables
    bool compute_powers_{false};
    Eigen::Vector<ftype, 2> kinetic_energy_, potential_energy_;
    ftype dissipated_power_, dissipated_power_flow_, dissipated_power_folds_;
    ftype external_power_, external_power_sub_, external_power_sup_;
    ftype stored_power_, stored_power_potential_, stored_power_kinetic_;

    // Resonator
    ftype sub_glottal_flow{0}, sup_glottal_flow{0};
    std::shared_ptr<WebsterFDTD<ftype>> resonator_;

    // Solver parameters
    ftype sr_, dt_;

    // Functions
    void RecomputeMatrices();
    void FillMassesInterpenetrationsAndAreas();
    void ComputeEffectiveAreas();
    void ComputeNonlinearDissipationVector();
    void ComputeSavVector();

public:
    Larynx(ftype samplerate, bool yielding_walls = false);
    void DspSetup(ftype sampleRate);

    void Process(ftype Pin);

    // "Listening" functions
    inline Eigen::Vector<ftype, 6> ReadFoldsDisplacement()
    {
        // Midpoint to evaluate on the same grid as inputs and momentums.
        return (q_(idx_now_, Eigen::all) + q_(idx_next_, Eigen::all)) * 0.5;
    };

    inline Eigen::Vector<ftype, 3> ReadEffectiveOpenings()
    {
        // Midpoint to evaluate on the same grid as inputs and momentums.
        return areas_below_masses_;
    };

    inline ftype ReadRadiatedPressure() { return resonator_->ReadRadiatedPressure(); }

    inline ftype ReadEpsilonSav() { return epsilon_sav_; }

    // Flow variables
    inline ftype ReadSupGlottalFlow() { return sup_glottal_flow; };
    inline ftype ReadMeanGlottalFlow() { return mean_flow_; };
    inline ftype ReadPressureDrop() { return Psub_(idx_now_) - Psup_; };

    // Power variables
    void set_compute_powers(bool compute_powers) { compute_powers_ = compute_powers; }
    std::tuple<ftype, ftype, ftype> ReadDissipatedPowers()
    {
        return {dissipated_power_, dissipated_power_flow_, dissipated_power_folds_};
    };
    std::tuple<ftype, ftype, ftype> ReadExchangedPowers()
    {
        return {external_power_, external_power_sub_, external_power_sup_};
    };
    std::tuple<ftype, ftype, ftype> ReadStoredPowers()
    {
        return {stored_power_, stored_power_kinetic_, stored_power_potential_};
    };

    // Compute and retrieve linear characteristics of the vocal fold model
    std::tuple<Eigen::Vector<ftype, 3>, Eigen::Vector<ftype, 3>, Eigen::Matrix<ftype, 3, 3>> GetModalCharacteristics();

    // Setters
    void set_rest_positions(const Eigen::Vector<ftype, 3>& rest_positions, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_.rest_positions_ = ClipEigen(rest_positions.array(), 1e-5, 1e-2); break;
        case kRight: right_vf_.rest_positions_ = ClipEigen(rest_positions.array(), 1e-5, 1e-2); break;
        default:
            left_vf_.rest_positions_ = right_vf_.rest_positions_ = ClipEigen(rest_positions.array(), 1e-5, 1e-2);
            break;
        }
    };
    void set_rest_positions(const ftype& lower, const ftype& upper, const ftype& body, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft:
            left_vf_.rest_positions_(0) = std::clamp(lower, ftype(1e-5), ftype(1e-2));
            left_vf_.rest_positions_(1) = std::clamp(upper, ftype(1e-5), ftype(1e-2));
            left_vf_.rest_positions_(2) = std::clamp(body, ftype(1e-5), ftype(1e-2));
            break;
        case kRight:
            right_vf_.rest_positions_(0) = std::clamp(lower, ftype(1e-5), ftype(1e-2));
            right_vf_.rest_positions_(1) = std::clamp(upper, ftype(1e-5), ftype(1e-2));
            right_vf_.rest_positions_(2) = std::clamp(body, ftype(1e-5), ftype(1e-2));
            break;
        default:
            left_vf_.rest_positions_(0) = std::clamp(lower, ftype(1e-5), ftype(1e-2));
            left_vf_.rest_positions_(1) = std::clamp(upper, ftype(1e-5), ftype(1e-2));
            left_vf_.rest_positions_(2) = std::clamp(body, ftype(1e-5), ftype(1e-2));
            right_vf_.rest_positions_(0) = std::clamp(lower, ftype(1e-5), ftype(1e-2));
            right_vf_.rest_positions_(1) = std::clamp(upper, ftype(1e-5), ftype(1e-2));
            right_vf_.rest_positions_(2) = std::clamp(body, ftype(1e-5), ftype(1e-2));
            break;
        }
    };

    void set_masses(const Eigen::Vector<ftype, 3>& masses, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_.masses_.diagonal() = ClipEigen(masses.array(), 1e-6, 1e-3); break;
        case kRight: right_vf_.masses_.diagonal() = ClipEigen(masses.array(), 1e-6, 1e-3); break;
        default:
            right_vf_.masses_.diagonal() = left_vf_.masses_.diagonal() = ClipEigen(masses.array(), 1e-6, 1e-3);
            break;
        }
        RecomputeMatrices();
    };
    void set_masses(const ftype& lower, const ftype& upper, const ftype& body, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft:
            left_vf_.masses_.diagonal()(0) = std::clamp(lower, ftype(1e-6), ftype(1e-3));
            left_vf_.masses_.diagonal()(1) = std::clamp(upper, ftype(1e-6), ftype(1e-3));
            left_vf_.masses_.diagonal()(2) = std::clamp(body, ftype(1e-6), ftype(1e-3));
            break;
        case kRight:
            right_vf_.masses_.diagonal()(0) = std::clamp(lower, ftype(1e-6), ftype(1e-3));
            right_vf_.masses_.diagonal()(1) = std::clamp(upper, ftype(1e-6), ftype(1e-3));
            right_vf_.masses_.diagonal()(2) = std::clamp(body, ftype(1e-6), ftype(1e-3));
            break;
        default:
            left_vf_.masses_.diagonal()(0) = std::clamp(lower, ftype(1e-6), ftype(1e-3));
            left_vf_.masses_.diagonal()(1) = std::clamp(upper, ftype(1e-6), ftype(1e-3));
            left_vf_.masses_.diagonal()(2) = std::clamp(body, ftype(1e-6), ftype(1e-3));

            right_vf_.masses_.diagonal()(0) = std::clamp(lower, ftype(1e-6), ftype(1e-3));
            right_vf_.masses_.diagonal()(1) = std::clamp(upper, ftype(1e-6), ftype(1e-3));
            right_vf_.masses_.diagonal()(2) = std::clamp(body, ftype(1e-6), ftype(1e-3));
            break;
        }
        RecomputeMatrices();
    };

    void set_lengths(const Eigen::Vector<ftype, 3>& lengths) { lengths_ = ClipEigen(lengths.array(), 1e-4, 1e-2); };
    void set_lengths(const ftype& lower, const ftype& upper, const ftype& body)
    {
        lengths_(0) = std::clamp(lower, ftype(1e-4), ftype(1e-2));
        lengths_(1) = std::clamp(upper, ftype(1e-4), ftype(1e-2));
        lengths_(2) = std::clamp(body, ftype(1e-4), ftype(1e-2));
    };

    void set_width(const ftype& width) { widths_(0) = widths_(1) = std::clamp(width, ftype(1e-3), ftype(1e-1)); };

    void set_stiffnesses(const Eigen::Vector<ftype, 4>& stiffnesses, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_.stiffnesses_.diagonal() = ClipEigen(stiffnesses.array(), 1e-6, 1e-3); break;
        case kRight: right_vf_.stiffnesses_.diagonal() = ClipEigen(stiffnesses.array(), 1e-6, 1e-3); break;
        default:
            right_vf_.stiffnesses_.diagonal() = left_vf_.stiffnesses_.diagonal() =
                ClipEigen(stiffnesses.array(), 1e-6, 1e-3);
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
        switch (fold_id) {
        case kLeft:
            left_vf_.stiffnesses_.diagonal()(0) = std::clamp(lower, ftype(0), ftype(10000));
            left_vf_.stiffnesses_.diagonal()(1) = std::clamp(upper, ftype(0), ftype(10000));
            left_vf_.stiffnesses_.diagonal()(2) = std::clamp(body, ftype(0), ftype(10000));
            left_vf_.stiffnesses_.diagonal()(3) = std::clamp(coupling, ftype(0), ftype(10000));
            break;
        case kRight:
            right_vf_.stiffnesses_.diagonal()(0) = std::clamp(lower, ftype(0), ftype(10000));
            right_vf_.stiffnesses_.diagonal()(1) = std::clamp(upper, ftype(0), ftype(10000));
            right_vf_.stiffnesses_.diagonal()(2) = std::clamp(body, ftype(0), ftype(10000));
            right_vf_.stiffnesses_.diagonal()(3) = std::clamp(coupling, ftype(0), ftype(10000));
            break;
        default:
            left_vf_.stiffnesses_.diagonal()(0) = std::clamp(lower, ftype(0), ftype(10000));
            left_vf_.stiffnesses_.diagonal()(1) = std::clamp(upper, ftype(0), ftype(10000));
            left_vf_.stiffnesses_.diagonal()(2) = std::clamp(body, ftype(0), ftype(10000));
            left_vf_.stiffnesses_.diagonal()(3) = std::clamp(coupling, ftype(0), ftype(10000));

            right_vf_.stiffnesses_.diagonal()(0) = std::clamp(lower, ftype(0), ftype(10000));
            right_vf_.stiffnesses_.diagonal()(1) = std::clamp(upper, ftype(0), ftype(10000));
            right_vf_.stiffnesses_.diagonal()(2) = std::clamp(body, ftype(0), ftype(10000));
            right_vf_.stiffnesses_.diagonal()(3) = std::clamp(coupling, ftype(0), ftype(10000));
            break;
        }
        RecomputeMatrices();
    };

    void set_eta_stiffness(const ftype& eta_stiffness, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_.eta_stiffness_ = std::clamp(eta_stiffness, ftype(0), ftype(1e12)); break;
        case kRight: right_vf_.eta_stiffness_ = std::clamp(eta_stiffness, ftype(0), ftype(1e12)); break;
        default:
            left_vf_.eta_stiffness_ = std::clamp(eta_stiffness, ftype(0), ftype(1e12));

            right_vf_.eta_stiffness_ = std::clamp(eta_stiffness, ftype(0), ftype(1e12));
            break;
        }
    }
    void set_contact_stiffness(const ftype contact_stiffness)
    {
        contact_stiffness_ = std::clamp(contact_stiffness, ftype(0), ftype(1e3));
    }
    void set_alpha_contact_stiffness(const ftype alpha_contact_stiffness)
    {
        alpha_contact_stiffness_ = std::clamp(alpha_contact_stiffness, ftype(1), ftype(5));
    }
    void set_eta_contact_stiffness(const ftype eta_contact_stiffness)
    {
        eta_contact_stiffness_ = std::clamp(eta_contact_stiffness, ftype(0), ftype(1e12));
    }
    void set_xi(const ftype xi, FoldIdentifier fold_id = kBoth)
    {
        switch (fold_id) {
        case kLeft: left_vf_.xi_ = std::clamp(xi, ftype(0), ftype(1e12)); break;
        case kRight: right_vf_.xi_ = std::clamp(xi, ftype(0), ftype(1e12)); break;
        default:
            left_vf_.xi_ = std::clamp(xi, ftype(0), ftype(1e12));

            right_vf_.xi_ = std::clamp(xi, ftype(0), ftype(1e12));
            break;
        }
    }

    void set_c0(const ftype c0) { c0_ = std::clamp(c0, ftype(300), ftype(380)); }
    void set_rho0(const ftype rho0) { rho0_ = std::clamp(rho0, ftype(1), ftype(15)); }
    void set_kt(const ftype kt) { kt_ = std::clamp(kt, ftype(1), ftype(2)); }

    // ftype set_lambda_sav() { return lambda_sav_; }

    // // Getters
    inline Eigen::Vector<ftype, 3> get_rest_positions(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_.rest_positions_;
        } else {
            return left_vf_.rest_positions_;
        }
    };

    inline Eigen::Vector<ftype, 3> get_masses(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_.masses_.diagonal();
        } else {
            return left_vf_.masses_.diagonal();
        }
    };
    inline Eigen::Vector<ftype, 3> get_lengths() { return lengths_; };
    inline Eigen::Vector<ftype, 3> get_widths() { return widths_; };
    inline Eigen::Vector<ftype, 4> get_stiffnesses(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_.stiffnesses_.diagonal();
        } else {
            return left_vf_.stiffnesses_.diagonal();
        }
    };
    inline ftype get_eta_stiffness(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_.eta_stiffness_;
        } else {
            return left_vf_.eta_stiffness_;
        }
    }

    inline ftype get_contact_stiffness() { return contact_stiffness_; }
    inline ftype get_alpha_contact_stiffness() { return alpha_contact_stiffness_; }
    inline ftype get_eta_contact_stiffness() { return eta_contact_stiffness_; }
    inline ftype get_xi(FoldIdentifier fold_id = kBoth)
    {
        if (fold_id == kRight) {
            return right_vf_.xi_;
        } else {
            return left_vf_.xi_;
        }
    }

    inline ftype get_c0() { return c0_; }
    inline ftype get_rho0() { return rho0_; }
    inline ftype get_kt() { return kt_; }

    inline ftype get_lambda_sav() { return lambda_sav_; }

    std::shared_ptr<WebsterFDTD<ftype>> get_resonator() { return resonator_; }
};

} // namespace tarte
