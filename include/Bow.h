#pragma once

#include <cmath>

namespace tarte {

enum BowMode { kMatusiak, kVigue, kTerrien };

template<class T>
class Bow {
private:
    BowMode bow_mode_ = kVigue;
    T mu_s_{0.4}, mu_c_{0.5}, mu_v_{0}, v_s_{0.1}, v_c_{0.02}, a_{100}, epsilon_{1e-4}, mu_d_{0.2}, n_{100}, vt_{0.01};

public:
    Bow() { };

    T Phi(T v_rel)
    {
        switch (bowMode) {
        case kMatusiak:
            return mu_c_ / 2 * M_PI * atan(vrel / v_c_) +
                   (mu_s_ - mu_c_) * sqrt(2 * a_) * vrel * exp(-a_ * vrel * vrel + 0.5) + mu_v_ * vrel;
            break;
        case kTerrien:
            return mu_d_ * tanh(4 * vrel / vt_) +
                   (mu_s_ - mu_d_) * vrel / vt_ / pow(0.25 * pow(vrel / vt_, 2) + 0.75, 2);
            break;
        default:
            return (mu_d_ * vrel * sqrt(vrel * vrel + epsilon_ / (n_ * n_)) +
                    2 * sqrt(mu_s_ * (mu_s_ - mu_d_)) / n_ * vrel) /
                   (vrel * vrel + 1 / (n_ * n_));
            break;
        }
    };

    // Getters
    BowMode get_bow_mode() const { return bow_mode_; }
    T get_mu_s() const { return mu_s_; }
    T get_mu_c() const { return mu_c_; }
    T get_mu_v() const { return mu_v_; }
    T get_v_s() const { return v_s_; }
    T get_v_c() const { return v_c_; }
    T get_a() const { return a_; }
    T get_epsilon() const { return epsilon_; }
    T get_mu_d() const { return mu_d_; }
    T get_n() const { return n_; }

    // Setters
    void set_bow_mode(BowMode mode) { bow_mode_ = mode; }
    void set_mu_s(T value) { mu_s_ = value; }
    void set_mu_c(T value) { mu_c_ = value; }
    void set_mu_v(T value) { mu_v_ = value; }
    void set_v_s(T value) { v_s_ = value; }
    void set_v_c(T value) { v_c_ = value; }
    void set_a(T value) { a_ = value; }
    void set_epsilon(T value) { epsilon_ = value; }
    void set_mu_d(T value) { mu_d_ = value; }
    void set_n(T value) { n_ = value; }

    // Convenience methods to set multiple Matusiak parameters at once
    void set_matusiak_parameters(T mu_s, T mu_c, T mu_v, T v_s, T v_c, T a)
    {
        mu_s_ = mu_s;
        mu_c_ = mu_c;
        mu_v_ = mu_v;
        v_s_ = v_s;
        v_c_ = v_c;
        a_ = a;
        bow_mode_ = kMatusiak;
    }

    // Convenience methods to set multiple Terrien parameters at once
    void set_terrien_parameters(T mu_s, T mu_d, T vt)
    {
        mu_s_ = mu_s;
        mu_d_ = mu_d;
        vt_ = vt;
        bow_mode_ = kTerrien;
    }

    // Convenience methods to set multiple Vigue parameters at once
    void set_vigue_parameters(T mu_s_, T mu_d, T n, T epsilon)
    {
        mu_s_ = mu_s_;
        mu_d_ = mu_d;
        n_ = n;
        epsilon_ = epsilon;
        bow_mode_ = kVigue;
    }
};
} // namespace tarte