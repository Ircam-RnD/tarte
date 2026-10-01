

# File eigen\_utility.h

[**File List**](files.md) **>** [**include**](dir_d44c64559bbebec7f509842c48db8b23.md) **>** [**utility**](dir_3a4d35156996fb70540e167b84a39bd1.md) **>** [**eigen\_utility.h**](eigen__utility_8h.md)

[Go to the documentation of this file](eigen__utility_8h.md)


```C++
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
```


