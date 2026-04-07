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

template<typename ftype>
class Larynx {
private:
    // Physical parameters
    // Vocal folds: order is lower-upper-body
    Eigen::DiagonalMatrix<ftype, 3> masses_{1e-5, 1e-5, 5e-5};
    Eigen::Vector<ftype, 3> lengths_{1.5e-3, 1.5e-3, 3e-3};
    Eigen::Vector<ftype, 3> widths_{1e-2, 1e-2, 0};
    Eigen::Vector<ftype, 3> rest_positions_{1.8e-4, 1.79e-4, 3e-3};

    Eigen::Vector<ftype, 4> elongations_;
    Eigen::DiagonalMatrix<ftype, 4> stiffnesses_{5, 3.5, 20, 0.5}; // lower, upper, body, lower-upper
    ftype eta_stiffness_{1e6};

    ftype xi_{0.4};
    Eigen::DiagonalMatrix<ftype, 4> dissipation_coefficients_;

    ftype contact_stiffness_{15}, eta_contact_stiffness_{1e6}, alpha_contact_stiffness_{3};

    ftype rho0_{1.2}, c0_{340}, kt_{1.3};

    // Matrices and intermediary quantities
    Eigen::Matrix<ftype, 4, 3> elongation_matrix_;
    Eigen::DiagonalMatrix<ftype, 3> mass_matrix_inv_;
    Eigen::Matrix<ftype, 3, 3> stiffness_matrix_, dissipation_matrix_;
    Eigen::Vector<ftype, 3> masses_interpenetrations_, areas_below_masses_, smoothed_is_opened_,
        masses_interpenetrations_derivatives_;
    Eigen::Vector<ftype, 3> effective_surfaces_Psub_, effective_surfaces_Psup_;

    Eigen::Vector<ftype, 3> g_sav_{0, 0, 0}, Fnl_{0, 0, 0};
    ftype Enl_, epsilon_sav_;
    bool const control_term_{true};
    ftype lambda_sav_{1000};

    ftype area_ratio_, area_min_;
    ftype mean_flow_, Rk_;

    ftype a_resonator_, b_resonator_;
    ftype C0_feedback_;
    Eigen::Vector<ftype, 3> C1_feedback_;

    ftype A0_inv_;
    Eigen::Vector<ftype, 3> rhs_;
    Eigen::Matrix<ftype, 2, 3> v_woodburry_;
    Eigen::Matrix<ftype, 3, 2> u_woodburry_;
    Eigen::Matrix<ftype, 2, 2> woodburry_inv_;

    ftype const kEpsilonSmooth_{1e-5};

    // State
    Eigen::Matrix<ftype, 2, 3> p_, q_;
    Eigen::Vector<ftype, 2> r_;
    std::size_t idx_now_{0}, idx_next_{1};

    Eigen::Vector<ftype, 2> Psub_, Psub_centered_;
    ftype Psup_;

    // Power variables
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
    void FillMassesInterpenetrationsAndAreas();
    void ComputeEffectiveAreas();
    void ComputeNonlienarDissipationVector();
    void ComputeSavVector();

public:
    Larynx(ftype samplerate, bool yielding_walls = false);

    void Process(ftype Pin);

    // "Listening" functions
    inline Eigen::Vector<ftype, 3> ReadFoldsDisplacement()
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
    std::tuple<ftype, ftype, ftype> ReadDissipatedPowers()
    {
        return {dissipated_power_, dissipated_power_flow_, dissipated_power_folds_};
    };
    std::tuple<ftype, ftype, ftype> ReadExcahngedPowers()
    {
        return {external_power_, external_power_sub_, external_power_sup_};
    };
    std::tuple<ftype, ftype, ftype> ReadStoredPowers()
    {
        return {stored_power_, stored_power_kinetic_, stored_power_potential_};
    };

    // Setters
    inline void set_rest_positions(const Eigen::Vector<ftype, 3> rest_positions)
    {
        rest_positions_ = ClipEigen(rest_positions.array(), 1e-5, 1e-2);
    };

    inline void set_masses(const Eigen::Vector<ftype, 3> masses)
    {
        masses_.diagonal() = ClipEigen(masses.array(), 1e-6, 1e-3);
    };
    inline void set_lengths(const Eigen::Vector<ftype, 3> lengths)
    {
        lengths_ = ClipEigen(lengths.array(), 1e-4, 1e-2);
    };
    inline void set_width(const ftype width) { widths_(0) = widths_(1) = std::clamp(width, ftype(1e-3), ftype(1e-1)); };
    inline void set_stiffnesses(const Eigen::Vector<ftype, 4> stiffnesses)
    {
        stiffnesses_.diagonal() = ClipEigen(stiffnesses.array(), ftype(0), ftype(1e2));
    };
    inline void set_eta_stiffness(const ftype eta_stiffness)
    {
        eta_stiffness_ = std::clamp(eta_stiffness, ftype(0), ftype(1e12));
    }
    inline void set_contact_stiffness(const ftype contact_stiffness)
    {
        contact_stiffness_ = std::clamp(contact_stiffness, ftype(0), ftype(1e3));
    }
    inline void set_alpha_contact_stiffness(const ftype alpha_contact_stiffness)
    {
        alpha_contact_stiffness_ = std::clamp(alpha_contact_stiffness, ftype(1), ftype(5));
    }
    inline void set_eta_contact_stiffness(const ftype eta_contact_stiffness)
    {
        eta_contact_stiffness_ = std::clamp(eta_contact_stiffness, ftype(0), ftype(1e12));
    }
    inline void set_xi(const ftype xi) { xi_ = std::clamp(xi, ftype(0), ftype(1)); }

    inline void set_c0(const ftype c0) { c0_ = std::clamp(c0, ftype(300), ftype(380)); }
    inline void set_rho0(const ftype rho0) { rho0_ = std::clamp(rho0, ftype(1), ftype(15)); }
    inline void set_kt(const ftype kt) { kt_ = std::clamp(kt, ftype(1), ftype(2)); }

    // inline ftype set_lambda_sav() { return lambda_sav_; }

    // Getters
    inline Eigen::Vector<ftype, 3> get_rest_positions() { return rest_positions_; };

    inline Eigen::Vector<ftype, 3> get_masses() { return masses_.diagonal(); };
    inline Eigen::Vector<ftype, 3> get_lengths() { return lengths_; };
    inline Eigen::Vector<ftype, 3> get_widths() { return widths_; };
    inline Eigen::Vector<ftype, 4> get_stiffnesses() { return stiffnesses_.diagonal(); };
    inline ftype get_eta_stiffness() { return eta_stiffness_; }
    inline ftype get_contact_stiffness() { return contact_stiffness_; }
    inline ftype get_alpha_contact_stiffness() { return alpha_contact_stiffness_; }
    inline ftype get_eta_contact_stiffness() { return eta_contact_stiffness_; }
    inline ftype get_xi() { return xi_; }

    inline ftype get_c0() { return c0_; }
    inline ftype get_rho0() { return rho0_; }
    inline ftype get_kt() { return kt_; }

    inline ftype get_lambda_sav() { return lambda_sav_; }

    std::shared_ptr<WebsterFDTD<ftype>> get_resonator() { return resonator_; }
};

} // namespace tarte
