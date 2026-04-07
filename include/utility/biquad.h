#pragma once

#include <algorithm>
#include <cmath>

namespace tarte {

enum BiquadMode { kLowPass, kHighPass, kBandPass, kNotch, kPeak, kLowShelf, kHighShelf };

class Biquad {
private:
    double x0_{0}, x1_{0}, x2_{0}, y0_{0}, y1_{0}, y2_{0};
    float sr_{44100};

    float freq_{100}, dB_gain_{0}, Q_{0.7};

    double A_, omega_, sin_omega_, cos_omega_, alpha_, beta_;

    double b0_{1}, b1_{0}, b2_{0}, a0_{0}, a1_{0}, a2_{0};

    BiquadMode mode{kLowPass};

    void ComputeInter()
    {
        if (mode < 4) {
            A_ = sqrt(pow(10, dB_gain_ / 20));
        } else {
            A_ = pow(10, dB_gain_ / 40);
        }
        omega_ = 2 * M_PI * freq_ / sr_;
        sin_omega_ = sin(omega_);
        cos_omega_ = cos(omega_);
        alpha_ = sin_omega_ / (2 * Q_);
        beta_ = sqrt(A_) / Q_;
    }

    void ComputeCoeffs()
    {
        switch (mode) {
        case kLowPass:
            b0_ = (1 - cos_omega_) / 2;
            b1_ = 1 - cos_omega_;
            b2_ = b0_;
            a0_ = 1 + alpha_;
            a1_ = -2 * cos_omega_;
            a2_ = 1 - alpha_;
            break;
        case kHighPass:
            b0_ = (1 + cos_omega_) / 2;
            b1_ = -1 - cos_omega_;
            b2_ = b0_;
            a0_ = 1 + alpha_;
            a1_ = -2 * cos_omega_;
            a2_ = 1 - alpha_;
            break;
        case kBandPass:
            b0_ = Q_ * alpha_;
            b1_ = 0;
            b2_ = -Q_ * alpha_;
            a0_ = 1 + alpha_;
            a1_ = -2 * cos_omega_;
            a2_ = 1 - alpha_;
            break;
        case kNotch:
            b0_ = 1;
            b1_ = -2 * cos_omega_;
            b2_ = 1;
            a0_ = 1 + alpha_;
            a1_ = -2 * cos_omega_;
            a2_ = 1 - alpha_;
            break;
        case kPeak:
            b0_ = 1 + alpha_ * A_;
            b1_ = -2 * cos_omega_;
            b2_ = A_ - alpha_ * A_;
            a0_ = 1 + alpha_ / A_;
            a1_ = -2 * cos_omega_;
            a2_ = 1 - alpha_ / A_;
            break;
        case kLowShelf:
            b0_ = A_ * ((A_ + 1) - (A_ - 1) * cos_omega_ + beta_ * sin_omega_);
            b1_ = 2 * A_ * ((A_ - 1) - (A_ + 1) * cos_omega_);
            b2_ = A_ * ((A_ + 1) - (A_ - 1) * cos_omega_ - beta_ * sin_omega_);
            a0_ = (A_ + 1) + (A_ - 1) * cos_omega_ + beta_ * sin_omega_;
            a1_ = -2 * ((A_ - 1) + (A_ + 1) * cos_omega_);
            a2_ = (A_ + 1) + (A_ - 1) * cos_omega_ - beta_ * sin_omega_;
            break;
        case kHighShelf:
            b0_ = A_ * ((A_ + 1) - (A_ - 1) * cos_omega_ + beta_ * sin_omega_);
            b1_ = -2 * A_ * ((A_ - 1) + (A_ + 1) * cos_omega_);
            b2_ = A_ * ((A_ + 1) - (A_ - 1) * cos_omega_ - beta_ * sin_omega_);
            a0_ = (A_ + 1) - (A_ - 1) * cos_omega_ + beta_ * sin_omega_;
            a1_ = 2 * ((A_ - 1) - (A_ + 1) * cos_omega_);
            a2_ = (A_ + 1) - (A_ - 1) * cos_omega_ - beta_ * sin_omega_;
            break;
        }
    }

public:
    Biquad(float sample_rate = 44100,
           BiquadMode mode = kLowPass,
           float frequency = 100,
           float dB_gain = 0,
           float Q = 0.7)
    {
        sr_ = sample_rate;
        set_mode(mode);
        set_freq(frequency);
        set_gain(dB_gain);
        set_Q(Q);
        ComputeInter();
        ComputeCoeffs();
    };

    void set_freq(float freq_)
    {
        this->freq_ = std::clamp(freq_, float(0), sr_);
        ComputeInter();
        ComputeCoeffs();
    }

    void set_gain(float dB_gain_)
    {
        this->dB_gain_ = std::clamp(dB_gain_, float(-60.0), float(20.0));
        ComputeInter();
        ComputeCoeffs();
    }

    void set_Q(float Q_)
    {
        this->Q_ = std::clamp(Q_, float(0.001), float(100.0));
        ComputeInter();
        ComputeCoeffs();
    }

    void set_mode(BiquadMode mode)
    {
        this->mode = mode;
        ComputeInter();
        ComputeCoeffs();
    }

    void InitializeState(double initial_value)
    {
        // Initialize filter states to steady-state value
        // For a step input, steady-state output equals DC gain * input
        double dc_gain = (b0_ + b1_ + b2_) / (a0_ + a1_ + a2_);
        double steady_state_outptu = dc_gain * initial_value;

        x0_ = x1_ = x2_ = initial_value;
        y0_ = y1_ = y2_ = steady_state_outptu;
    }

    double Process(double input)
    {
        x0_ = input;
        y0_ = (b0_ / a0_) * x0_ + (b1_ / a0_) * x1_ + (b2_ / a0_) * x2_ - (a1_ / a0_) * y1_ - (a2_ / a0_) * y2_;

        y2_ = y1_;
        y1_ = y0_;

        x2_ = x1_;
        x1_ = x0_;
        return y0_;
    }
};

} // namespace tarte