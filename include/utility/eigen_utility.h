#pragma once

#include <Eigen/Dense>

namespace tarte {

template<typename Derived, typename ftype>
auto ClipEigen(Eigen::Ref<const Eigen::ArrayX<Derived>> array, const ftype& min, const ftype& max)
{
    return array.cwiseMin(max).cwiseMax(min);
}

template<typename Derived, typename ftype>
auto ClipEigen(const Eigen::ArrayBase<Derived>& array, const ftype& min, const ftype& max)
{
    return array.cwiseMin(max).cwiseMax(min);
}

template<typename Derived1, typename Derived2>
auto SafeSetEigen(Eigen::ArrayBase<Derived1>& array, Eigen::Ref<const Eigen::ArrayX<Derived2>> array2)
{
    int minDim = std::min(array.size(), array2.size());
    array.head(minDim) = array2.head(minDim);
}

template<typename Derived1, typename Derived2>
auto SafeSetEigen(Eigen::ArrayBase<Derived1>& array, const Eigen::ArrayBase<Derived2>& array2)
{
    int minDim = std::min(array.size(), array2.size());
    array.head(minDim) = array2.head(minDim);
}

} // namespace tarte