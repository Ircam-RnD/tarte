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

    ftype smoothed_is_opened_;

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
    std::shared_ptr<WebsterFDTD<ftype, 300>> resonator_;

    // Solver parameters
    ftype sr_, dt_;
    ftype epsilon_smooth_{1e-7}, epsilon_smooth_P_{10};

    // Functions
    void fillOpeningAndInterpenetration();
    void computeRk();
    void computegSAV();

    // Experimental noise
    float noise_ratio_{0.f}, random_value_;
    ftype noise_flow_;

public:
    SingleReed(float samplerate);

    void DspSetup(float samplerate);

    void Process(float Pin);

    inline ftype ReadDisplacement()
    {
        // Midpoint to evaluate on the same grid as inputs and momentums.
        return (q_(idx_now_) + q_(idx_next_)) * 0.5;
    };

    inline ftype ReadEffectiveOpening() { return opening_; };

    inline ftype ReadRadiatedPressure() { return resonator_->ReadRadiatedPressure(); }
    inline ftype ReadAuxiliaryVariable() { return r_(idx_now_); };

    // Flow variables
    inline ftype ReadMeanMouthpieceFlow() { return mean_flow_; };
    inline ftype ReadMouthpieceFlow() { return resonator_flow_; };
    inline ftype ReadPressureDrop() { return Psub_(idx_now_) - Psup_; };
    inline ftype getCurrentState() { return Psub_(idx_now_) - Psup_; };

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

    // Setters
    void set_mass(const ftype& mass)
    {
        mass_ = std::clamp(mass, ftype(1e-5), ftype(1e-3));
        dissipation_coefficient_ = 2 * mass_ * damping_;
    }
    void set_stiffness(const ftype& stiffness) { stiffness_ = std::clamp(stiffness, ftype(1e1), ftype(1e4)); }
    void set_damping(const ftype& damping)
    {
        damping_ = std::clamp(damping, ftype(0), ftype(10000));
        dissipation_coefficient_ = 2 * mass_ * damping_;
    }
    void set_lay_position(const ftype& lay_position)
    {
        lay_position_ = std::clamp(lay_position, ftype(0), ftype(1e-2));
    }
    void set_width(const ftype& width) { width_ = std::clamp(width, ftype(1e-3), ftype(3e-2)); }
    void set_surface(const ftype& surface) { surface_ = std::clamp(surface, ftype(1e-6), ftype(1e-3)); }
    void set_contact_stiffness(const ftype& contact_stiffness)
    {
        contact_stiffness_ = std::clamp(contact_stiffness, ftype(0), ftype(1e13));
    }
    void set_epsilon_smooth(const ftype& epsilon_smooth)
    {
        epsilon_smooth_ = std::clamp(epsilon_smooth, ftype(1e-8), ftype(1e-3));
    }
    void set_noise_ratio(const ftype& noise_ratio) { noise_ratio_ = std::clamp(noise_ratio, ftype(0), ftype(1)); }

    // Getters
    ftype get_mass() { return mass_; }
    ftype get_stiffness() { return stiffness_; }
    ftype get_damping() { return damping_; }
    ftype get_width() { return width_; }
    ftype get_surface() { return surface_; }
    ftype get_lay_position() { return lay_position_; }
    ftype get_contact_stiffness() { return contact_stiffness_; }
    ftype get_epsilon_smooth() { return epsilon_smooth_; }
    ftype get_noise_ratio() { return noise_ratio_; }

    std::shared_ptr<WebsterFDTD<ftype, 300>> get_resonator() { return resonator_; }
};

} // namespace tarte
