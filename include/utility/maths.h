#pragma once
#include <vector>

namespace tarte {
template<typename T>
void normalize_vector(std::vector<T>& vec)
{
    T max = 0;
    for (int i = 0; i < vec.size(); i++) {
        T abs_val = abs(vec.at(i));
        if (abs_val > max)
            max = abs_val;
    }
    for (int i = 0; i < vec.size(); i++) {
        vec.at(i) /= max;
    }
};

template<typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}
} // namespace tarte