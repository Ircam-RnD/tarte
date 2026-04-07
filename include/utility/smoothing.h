#pragma once

#include <Eigen/Dense>

namespace tarte {

template<typename derived, typename ftype>
auto smoothRampMatrix(const Eigen::MatrixBase<derived>& x, ftype epsilon = 1)
{
    const auto t = ((x.array()) / (epsilon)).cwiseMax(0).cwiseMin(1);
    return (x.array() > epsilon).select(x.array(), (-t * t * t + 2 * t * t) * epsilon).matrix();
}

template<typename ftype>
ftype softplus(const ftype& x, ftype epsilon = 1, ftype threshold = 40)
{
    if (x / epsilon > threshold) {
        return x;
    } else {
        return epsilon * log1p(exp(x / epsilon));
    }
}

template<typename ftype>
ftype softplusDerivative(const ftype& x, ftype epsilon = 1, ftype threshold = 40)
{
    if (x / epsilon > threshold) {
        return 1;
    } else {
        return exp(x / epsilon) / (1 + exp(x / epsilon));
    }
}

template<typename derived, typename ftype>
auto softplusMatrix(const Eigen::MatrixBase<derived>& x, ftype epsilon = 1, ftype threshold = 40)
{
    return ((x / epsilon).array() > threshold).select(x, (x / epsilon).array().exp().log1p().matrix() * epsilon);
}

template<typename derived, typename ftype>
auto softplusDerivativeMatrix(const Eigen::MatrixBase<derived>& x, ftype epsilon = 1, ftype threshold = 40)
{
    return ((x / epsilon).array() > threshold)
        .select(Eigen::ArrayX<ftype>::Ones(x.rows(), x.cols()),
                (x / epsilon).array().exp() / (1 + (x / epsilon).array().exp()));
}
} // namespace tarte
