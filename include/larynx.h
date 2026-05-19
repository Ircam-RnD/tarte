#pragma once

#include "body_cover_vf.h"
#include "utility/eigen_utility.h"
#include "utility/maths.h"
#include "utility/smoothing.h"
#include "websterFDTD.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <memory>
#include <string_view>
#include <tuple>
#include <vector>

namespace tarte {

enum FoldIdentifier { kRight, kLeft, kBoth };

template<typename ftype>
class Larynx {
private:
    // Folds physical parameters are in there
    std::shared_ptr<BodyCoverVF<ftype>> left_vf_, right_vf_;

    // Contact parameters
    ftype contact_stiffness_{15}, eta_contact_stiffness_{1e6}, alpha_contact_stiffness_{1.3};

    ftype rho0_{1.2}, c0_{340}, kt_{1.3};

    // Matrices and intermediary quantities
    Eigen::DiagonalMatrix<ftype, 6> mass_matrix_inv_;                              // Both folds
    Eigen::Matrix<ftype, 3, 3> stiffness_matrix_left_, dissipation_matrix_left_;   // Left fold only
    Eigen::Matrix<ftype, 3, 3> stiffness_matrix_right_, dissipation_matrix_right_; // Right fold only
    Eigen::Vector<ftype, 4> elongations_;

    // The following quantities need to be defined only once for both folds, as they are related to the fluid flow
    Eigen::Vector<ftype, 3> masses_interpenetrations_, areas_below_masses_, smoothed_is_opened_,
        masses_interpenetrations_derivatives_;
    // The effective surfaces must however be fold dependant if the folds have different geometries
    Eigen::Vector<ftype, 3> effective_surfaces_Psub_left_, effective_surfaces_Psup_left_;
    Eigen::Vector<ftype, 3> effective_surfaces_Psub_right_, effective_surfaces_Psup_right_;

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

    ftype epsilon_smooth_{1e-5};

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
    void RecomputeMatrices(bool update_from_muscles = false);
    void FillMassesInterpenetrationsAndAreas();
    void ComputeEffectiveAreas();
    void ComputeNonlinearDissipationVector();
    void ComputeSavVector();

public:
    Larynx(ftype samplerate, bool yielding_walls = false);
    void DspSetup(ftype sampleRate, Articulation* art = nullptr);

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
    // std::tuple<Eigen::Vector<ftype, 3>, Eigen::Vector<ftype, 3>, Eigen::Matrix<ftype, 3, 3>>
    // GetModalCharacteristics();

    // Setters
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
    void set_eta_contact_stiffness(const ftype& eta_contact_stiffness)
    {
        eta_contact_stiffness_ = std::clamp(eta_contact_stiffness, ftype(0), ftype(1e12));
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

    void set_lambda_sav(const ftype& lambda_sav) { lambda_sav_ = std::clamp(lambda_sav, ftype(0), sr_); }
    void set_epsilon_smooth(const ftype& epsilon_smooth)
    {
        epsilon_smooth_ = std::clamp(epsilon_smooth, ftype(1e-6), ftype(1e-2));
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
    inline ftype get_eta_contact_stiffness() { return eta_contact_stiffness_; }
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

    inline ftype get_lambda_sav() { return lambda_sav_; }
    inline ftype get_epsilon_smooth() { return epsilon_smooth_; }

    std::shared_ptr<WebsterFDTD<ftype>> get_resonator() { return resonator_; }
};

} // namespace tarte
