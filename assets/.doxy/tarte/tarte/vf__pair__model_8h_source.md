

# File vf\_pair\_model.h

[**File List**](files.md) **>** [**include**](dir_d44c64559bbebec7f509842c48db8b23.md) **>** [**vocal\_folds**](dir_35d2a1600e0c11c335c654e0ec7e535e.md) **>** [**vf\_pair\_model.h**](vf__pair__model_8h.md)

[Go to the documentation of this file](vf__pair__model_8h.md)


```C++
template<typename T>
concept VFPairModel = requires(T model,
                               const typename T::state_type& state,
                               typename T::state_type& out_state,
                               const typename T::scalar_type& P,
                               bool recompute)
{
    typename T::scalar_type;
    typename T::state_type;

    {T::get_N()}->std::convertible_to<int>;
    {T::get_half_N()}->std::convertible_to<int>;

    {model.FillIntermediary(state)}->std::same_as<void>;
    {model.EffectiveAreas(out_state, out_state)}->std::same_as<void>;
    {model.ComputeFlow(P, P)}->std::same_as<typename T::scalar_type>;

    {model.Enl(state)}->std::same_as<typename T::scalar_type>;
    {model.Enl(state, recompute)}->std::same_as<typename T::scalar_type>;

    {model.Fnl(state, out_state)}->std::same_as<void>;
    {model.Fnl(state, out_state, recompute)}->std::same_as<void>;

    {model.KOp(state, out_state)}->std::same_as<void>;
    {model.ROp(state, out_state)}->std::same_as<void>;
    {model.MinvOp(state)}->std::same_as<typename T::state_type>;

    {model.ReadEffectiveOpenings()}->std::same_as<typename T::half_state_type>;
};
```


