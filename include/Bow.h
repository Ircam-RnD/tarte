#pragma once

#include <cmath>

namespace tarte {

enum BOWMODE { MATUSIAK, VIGUE, TERRIEN };

template<class T>
class Bow {
private:
    BOWMODE bowMode = VIGUE;
    T mu_s{0.4}, mu_c{0.5}, mu_v{0}, v_s{0.1}, v_c{0.02}, a{100}, epsilon{1e-4}, mu_d{0.2}, n{100}, vt{0.01};

public:
    Bow() { };

    T phi(T vrel)
    {
        switch (bowMode) {
        case MATUSIAK:
            return mu_c / 2 * M_PI * atan(vrel / v_c) +
                   (mu_s - mu_c) * sqrt(2 * a) * vrel * exp(-a * vrel * vrel + 0.5) + mu_v * vrel;
            break;
        case TERRIEN:
            return mu_d * tanh(4 * vrel / vt) + (mu_s - mu_d) * vrel / vt / pow(0.25 * pow(vrel / vt, 2) + 0.75, 2);
            break;
        default:
            return (mu_d * vrel * sqrt(vrel * vrel + epsilon / (n * n)) + 2 * sqrt(mu_s * (mu_s - mu_d)) / n * vrel) /
                   (vrel * vrel + 1 / (n * n));
            break;
        }
    };

    // Getters
    BOWMODE getBowMode() const { return bowMode; }
    T getMuS() const { return mu_s; }
    T getMuC() const { return mu_c; }
    T getMuV() const { return mu_v; }
    T getVS() const { return v_s; }
    T getVC() const { return v_c; }
    T getA() const { return a; }
    T getEpsilon() const { return epsilon; }
    T getMuD() const { return mu_d; }
    T getN() const { return n; }

    // Setters
    void setBowMode(BOWMODE mode) { bowMode = mode; }
    void setMuS(T value) { mu_s = value; }
    void setMuC(T value) { mu_c = value; }
    void setMuV(T value) { mu_v = value; }
    void setVS(T value) { v_s = value; }
    void setVC(T value) { v_c = value; }
    void setA(T value) { a = value; }
    void setEpsilon(T value) { epsilon = value; }
    void setMuD(T value) { mu_d = value; }
    void setN(T value) { n = value; }

    // Convenience methods to set multiple Matusiak parameters at once
    void setMatusiakParams(T mu_s_val, T mu_c_val, T mu_v_val, T v_s_val, T v_c_val, T a_val)
    {
        mu_s = mu_s_val;
        mu_c = mu_c_val;
        mu_v = mu_v_val;
        v_s = v_s_val;
        v_c = v_c_val;
        a = a_val;
        bowMode = MATUSIAK;
    }

    // Convenience methods to set multiple Terrien parameters at once
    void setTerrienParams(T mu_s_val, T mu_d_val, T vt_val)
    {
        mu_s = mu_s_val;
        mu_d = mu_d_val;
        vt = vt_val;
        bowMode = TERRIEN;
    }

    // Convenience methods to set multiple Vigue parameters at once
    void setVigueParams(T mu_s_val, T mu_d_val, T n_val, T epsilon_val)
    {
        mu_s = mu_s_val;
        mu_d = mu_d_val;
        n = n_val;
        epsilon = epsilon_val;
        bowMode = VIGUE;
    }
};
} // namespace tarte