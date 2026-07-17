#pragma once

#include <Eigen/Dense>
#include <tuple>

namespace tarte {

// This file implements an empty class that satisfies the minimum requriements of a vocal fold pair class, as checked in
// vf_pair_model.h

template<typename ftype>
class VFPairExample {
private:
    static const int N{6};
    static const int half_N{3};

public:
    using scalar_type = ftype;
    using state_type = Eigen::Vector<ftype, N>;           // Either q or p
    using half_state_type = Eigen::Vector<ftype, half_N>; // Quantities common to both folds
    static inline constexpr int get_N() { return N; };
    static inline constexpr int get_half_N() { return half_N; };
    VFPairExample() { };

    void FillIntermediary(const state_type& state_q) { };
    void EffectiveAreas(state_type& out_area_P_sub, state_type& out_area_P_sup) { };
    ftype ComputeFlow(const ftype& Psub, const ftype& Psup) { return ftype(0); };
    ftype Enl(const state_type& state_q, bool recompute_intermediary = false) { return ftype(0); };
    void Fnl(const state_type& state_q, state_type& out_area_P_sup, bool recompute_intermediary = false) { };
    void KOp(const state_type& state_q, state_type& out) { };
    void ROp(const state_type& state_p, state_type& out) { };
    state_type MinvOp(const state_type& state_p) { return state_p; };
    half_state_type ReadEffectiveOpenings() { return Eigen::Vector<ftype, half_N>::Zero(); };
};

} // namespace tarte