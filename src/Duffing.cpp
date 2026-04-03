#include "Duffing.h"

namespace tarte {

template<typename ftype>
Duffing<ftype>::Duffing(float sampleRate)
{
    ReinitDsp(sampleRate);
};

template<typename ftype>
void Duffing<ftype>::ReinitDsp(float sampleRate)
{
    sr = sampleRate;
    dt = 1 / (float(sr));

    // Reinit state
    qnow = 0;
    qnext = qnow;
    qlast = qnow;

    r = 0;

    RHS = LHS = qnow;

    // Reinit nonlinear variables
    g = qnow;
    dqV = qnow;
    V = 0;

    // Reinit system matrices
    M = K = R0 = 0;

    Amp = M;
    Omega = 2 * M_PI * 100 * Amp;
    Decay = M;
};

template<typename ftype>
void Duffing<ftype>::setModalMatrices()
{
    R0 = 1 / Omega;
    M = 1 / (2 * 6.9) * Decay * R0;
    K = M * Omega * Omega;
};

template<typename ftype>
void Duffing<ftype>::setPhysicalParameters(ftype M, ftype K, ftype R)
{
    R0 = R;
    this->M = M;
    this->K = K;
};

template<typename ftype>
void Duffing<ftype>::setLinearParameters(ftype Amp, ftype Omega, ftype Decay)
{
    this->Amp = Amp;
    this->Omega = Omega;
    this->Decay = Decay;

    setModalMatrices();
};

template<typename ftype>
void Duffing<ftype>::setAmp(ftype Amp)
{
    this->Amp = Amp;
    setModalMatrices();
};

template<typename ftype>
void Duffing<ftype>::setFreq(ftype Freq)
{
    Omega = 2 * M_PI * Freq;
    setModalMatrices();
};

template<typename ftype>
void Duffing<ftype>::setDecay(ftype Decay)
{
    this->Decay = Decay;
    setModalMatrices();
};

template<typename ftype>
void Duffing<ftype>::computeVAndVprime(ftype q)
{
    V = K * eta / 4 * pow(q, 4);
    dqV = K * eta * pow(q, 3);
};

template<typename ftype>
void Duffing<ftype>::computeV(ftype q)
{
    V = K * eta / 4 * pow(q, 4);
};

template<typename ftype>
std::tuple<ftype, ftype, ftype> Duffing<ftype>::process(ftype forceIn)
{
    // Nonlinear part
    computeVAndVprime(qnow);
    g = dqV / (sqrt(2 * V) + NUM_EPS);
    if (V > maxV) {
        maxV = V;
    }

    if (controlTerm) {
        computeV((qnow + qlast) / 2);
        epsilon = r - sqrt(2 * V);
        g -= lambda0 * epsilon * dt * ((ftype(0) < (qnow - qlast)) - ((qnow - qlast) < ftype(0))) /
             (abs(qnow - qlast) + NUM_EPS);
    }

    // Linear part
    RHS = (-K + 2 * M / (dt * dt)) * (qnow) - (M / (dt * dt) - (R0) / (2 * dt)) * (qlast) + forceIn;

    // Nonlinear part
    RHS += 0.25 * g * g * qlast - g * r;

    // Solve using Shermann-Morrisson
    auto A0_inv = dt * dt / (M + dt * (R0) / 2);
    qnext = A0_inv * RHS - (A0_inv * g * A0_inv * g * RHS) / (4 + A0_inv * g * g);

    rlast = r;
    r = r + 0.5 * g * (qnext - qlast);

    qlast = qnow;
    qnow = qnext;

    return {(qnow + qlast) / 2, (qnow - qlast) / dt * M, rlast};
};

template class Duffing<double>;
template class Duffing<float>;

} // namespace tarte