#pragma once

#include <algorithm>
#include <cmath>
#include <tuple>

#include <export.h>

namespace tarte {

template<typename ftype>
class Duffing {
private:
    // Numerical epsilon value
    constexpr static ftype NUM_EPS{1e-14};

    // Linear part: system parameters (diagonal)
    ftype M, K, R0;
    // Higer level modal parameters
    ftype Amp, Omega, Decay;

    // Nonlinear parameters
    ftype eta{0};

    // Time-scheme parameters
    float sr;
    ftype dt;
    bool controlTerm{true};
    ftype lambda0{0};

    // System state
    ftype qlast, qnow, qnext;
    ftype r, rlast;

    // Intermediate vectors
    ftype RHS, LHS;

    // Nonlinear function evaluation
    ftype g, dqV;
    ftype V;

    // Drift variable
    ftype epsilon{0}, maxV{0};

    void setModalMatrices();

public:
    Duffing(float sampleRate);

    void ReinitDsp(float sampleRate);

    void computeVAndVprime(ftype q);

    void computeV(ftype q);

    std::tuple<ftype, ftype, ftype> process(ftype forceIn);

    // Physical parameters
    void setPhysicalParameters(ftype M, ftype K, ftype R);

    // Higher level modal parameters
    void setLinearParameters(ftype Amp, ftype Omega, ftype Decay);
    void setAmp(ftype Amp);
    void setFreq(ftype Freq);
    void setDecay(ftype Decay);

    void setLambda0(ftype lambda0) { this->lambda0 = std::clamp(lambda0, ftype(0), ftype(10000)); };
};

} // namespace tarte