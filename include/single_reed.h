#pragma once
#define LARYNX_H

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
class SingleReed {
private:
    // Physical parameters
    ftype mass_{3.37e-6}, stiffness_{1821}, damping_{1500}, width_{1e-2}, surface_{1.46e-4},
        lay_position_{4e-4};                                       // Reed
    ftype contact_stiffness_{1e13}, alpha_contact_stiffness_{1.3}; // Contact
    ftype rho_0_{1.2}, c_0_{340}, kt_{1.3};                        // Fluid

    // Intermediary quantities
    ftype dissipation_coefficient_; // R_0 matrix
    ftype opening_, interpenetration_;
    ftype interpenetration_derivative_;

    ftype g_sav_{0}, E_nl_{0}, F_nl_{0};

    ftype mean_flow_, R_k_;

    ftype a_resonator_, b_resonator_;
    ftype C0_feedback_, C1_feedback_;

    ftype rhs_;

    ftype epsilon_smooth_{1e-7}, epsilon_smooth_P_{1};

    // State
    Eigen::Vector<ftype, 2> p_, q_, r_;
    ftype p_inter_, q_inter_;
    std::size_t idx_now_{0}, idx_next_{1};

    Eigen::Vector<ftype, 2> Psub_, Psub_centered_;
    ftype Psup_;

    // Power variables
    bool compute_powers_{false};
    Eigen::Vector<ftype, 2> kinetic_energy_, potential_energy_;
    ftype dissipated_power_, dissipated_power_flow_, dissipated_power_reed_;
    ftype external_power_, external_power_sub_, external_power_sup_;
    ftype stored_power_, stored_power_potential_, stored_power_kinetic_;

    // Resonator
    ftype mouth_flow_{0}, resonator_flow_{0};
    std::shared_ptr<WebsterFDTD<ftype>> resonator_;

    // Solver parameters
    ftype sr_, dt_;

    // Functions
    void fillOpeningAndInterpenetration();
    void computeRk();
    void computegSAV();

public:
    SingleReed(float samplerate);

    void Process(float Pin);

    inline ftype ReadDisplacement()
    {
        // Midpoint to evaluate on the same grid as inputs and momentums.
        return (q_(idx_now_) + q_(idx_next_)) * 0.5;
    };

    inline ftype ReadOpening() { return opening_; };

    inline ftype ReadRadiatedPressure() { return resonator_->ReadRadiatedPressure(); }
    inline ftype ReadAuxiliaryVariable() { return r_(idx_now_); };

    inline ftype get_lay_position() { return lay_position_; };

    // Power variables
    std::tuple<ftype, ftype, ftype> getCurrentDissipatedPowers()
    {
        return {dissipated_power_, dissipated_power_flow_, dissipated_power_reed_};
    };
    std::tuple<ftype, ftype, ftype> getCurrentExchangedPowers()
    {
        return {external_power_, external_power_sub_, external_power_sup_};
    };
    std::tuple<ftype, ftype, ftype> getCurrentStoredPowers()
    {
        return {stored_power_, stored_power_kinetic_, stored_power_potential_};
    };

    // Flow variables
    inline ftype getCurrentMeanFlow() { return mean_flow_; };
    inline ftype getCurrentResonatorFlow() { return resonator_flow_; };
    inline ftype getCurrentPressureDrop() { return Psub_(idx_now_) - Psup_; };
    inline ftype getCurrentState() { return Psub_(idx_now_) - Psup_; };

    std::shared_ptr<WebsterFDTD<ftype>> get_resonator() { return resonator_; }
};

} // namespace tarte
