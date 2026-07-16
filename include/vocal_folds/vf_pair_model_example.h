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
    using state_type = Eigen::Vector<ftype, N>; // Either q or p
    static inline constexpr int get_N() { return N; };
    static inline constexpr int get_half_N() { return half_N; };
    VFPairExample() { };

    void FillIntermediary(const state_type& state_q) { };
    std::tuple<state_type, state_type> EffectiveAreas()
    {
        state_type EffectiveAreaPsub;
        state_type EffectiveAreaPsup;
        return {EffectiveAreaPsub, EffectiveAreaPsup};
    };
    ftype ComputeFlow(const ftype& deltaP) { return ftype(0); };
    ftype Enl(const state_type& state_q, bool recompute_intermediary = false) { return ftype(0); };
    state_type Fnl(const state_type& state_q, bool recompute_intermediary = false) { return state_q; };
    state_type KOp(const state_type& state_q) { return state_q; };
    state_type ROp(const state_type& state_p) { return state_p; };
    state_type MinvOp(const state_type& state_p) { return state_p; };
};

template class VFPairExample<float>;
template class VFPairExample<double>;

} // namespace tarte