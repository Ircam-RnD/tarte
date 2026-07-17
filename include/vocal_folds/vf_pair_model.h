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