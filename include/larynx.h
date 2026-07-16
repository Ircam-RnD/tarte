#pragma once

#include "noise_generator.h"
#include "utility/eigen_utility.h"
#include "utility/maths.h"
#include "utility/smoothing.h"
#include "vocal_folds/body_cover_vf.h"
#include "vocal_folds/vf_pair_model.h"
#include "vocal_folds/vf_pair_model_example.h"
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

template<VFPairModel vf_pair, typename ftype>
class Larynx {
private:
    // Folds physical parameters are in there
    std::shared_ptr<vf_pair> vfs_;

    // std::shared_ptr<BodyCoverVF<ftype>> left_vf_, right_vf_;

    // // Contact parameters
    // ftype contact_stiffness_{15}, alpha_contact_stiffness_{1.3};

    ftype rho0_{1.2}, c0_{340}, kt_{1.3};

    // Matrices and intermediary quantities
    // Eigen::DiagonalMatrix<ftype, 6> mass_matrix_inv_;                              // Both folds
    // Eigen::Matrix<ftype, 3, 3> stiffness_matrix_left_, dissipation_matrix_left_;   // Left fold only
    // Eigen::Matrix<ftype, 3, 3> stiffness_matrix_right_, dissipation_matrix_right_; // Right fold only
    // Eigen::Vector<ftype, 4> elongations_;

    // The following quantities need to be defined only once for both folds, as they are related to the fluid flow
    // Eigen::Vector<ftype, 3> masses_interpenetrations_, areas_below_masses_, smoothed_is_opened_,
    //     masses_interpenetrations_derivatives_;
    // The effective surfaces must however be fold dependant if the folds have different geometries
    Eigen::Vector<ftype, vf_pair::get_N()> effective_surfaces_Psub_, effective_surfaces_Psup_;

    Eigen::Vector<ftype, vf_pair::get_N()> g_sav_{0, 0, 0, 0, 0, 0}, Fnl_{0, 0, 0, 0, 0, 0};
    ftype Enl_, epsilon_sav_;
    bool const control_term_{true};
    ftype lambda_sav_{10};

    // ftype area_ratio_, area_min_;
    ftype mean_flow_, Rk_;

    ftype a_resonator_, b_resonator_;
    ftype C0_feedback_;
    Eigen::Vector<ftype, 6> C1_feedback_;

    ftype A0_inv_;
    Eigen::Vector<ftype, vf_pair::get_N()> rhs_;
    Eigen::Matrix<ftype, 2, vf_pair::get_N()> v_woodburry_;
    Eigen::Matrix<ftype, vf_pair::get_N(), 2> u_woodburry_;
    Eigen::Matrix<ftype, 2, 2> woodburry_inv_;
    Eigen::Vector<ftype, vf_pair::get_N()> Minv_p_;

    ftype epsilon_smooth_{1e-5};

    // State
    Eigen::Matrix<ftype, 2, vf_pair::get_N()> p_, q_;
    Eigen::Vector<ftype, vf_pair::get_N()> pmid_, qmid_;
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

    // Experimental noise
    NoiseGenerator noise_generator_ = NoiseGenerator(NoiseColor::Pink);
    float noise_ratio_{0.f};
    ftype noise_flow_;

    // Functions
    // void FillMassesInterpenetrationsAndAreas();
    // void ComputeEffectiveAreas();
    void ComputeNonlinearDissipationCoefficient();
    void ComputeSavVector();

public:
    Larynx(ftype samplerate, bool yielding_walls = false);
    void DspSetup(ftype sampleRate, Articulation* art = nullptr);

    void Process(ftype Pin);

    // "Listening" functions
    // inline Eigen::Vector<ftype, half_state> ReadFoldsDisplacement()
    // {
    //     // Midpoint to evaluate on the same grid as inputs and momentums.
    //     return (q_(idx_now_, Eigen::placeholders::all) + q_(idx_next_, Eigen::placeholders::all)) * 0.5;
    // };

    // inline Eigen::Vector<ftype, 3> ReadEffectiveOpenings()
    // {
    //     // Midpoint to evaluate on the same grid as inputs and momentums.
    //     return areas_below_masses_;
    // };

    // inline ftype ReadRadiatedPressure() { return resonator_->ReadRadiatedPressure(); }

    // inline ftype ReadEpsilonSav() { return epsilon_sav_; }

    // // Flow variables
    // inline ftype ReadSupGlottalFlow() { return sup_glottal_flow; };
    // inline ftype ReadMeanGlottalFlow() { return mean_flow_; };
    // inline ftype ReadPressureDrop() { return Psub_(idx_now_) - Psup_; };

    std::shared_ptr<WebsterFDTD<ftype>> get_resonator() { return resonator_; }
};

} // namespace tarte
