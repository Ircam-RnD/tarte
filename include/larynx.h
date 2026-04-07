#pragma once

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

    ftype contact_stiffness_{15}, eta_contact_stiffness{1e6}, alpha_contact_stiffness{3};

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
    bool control_term_{true};
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

    ftype epsilon_smooth_{1e-5};

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

    // Setters and getters

    inline Eigen::Vector<ftype, 3> get_rest_positions() { return rest_positions_; };

    std::shared_ptr<WebsterFDTD<ftype>> get_resonator() { return resonator_; }
};

} // namespace tarte
