#include <concepts>
#include <tuple>

namespace tarte {

// This file implements a concept for the vocal fold pair model to be provided to the voice simulation. Basically, it
// serves as a compile time check of the valitidy of a vocal fold model, hopefully yielding clearer errors at compile
// time.

template<typename T>
concept VFPairModel =
    requires(T model, const typename T::state_type& state, const typename T::scalar_type& deltaP, bool recompute)
{
    // Required nested types
    typename T::scalar_type;
    typename T::state_type;

    // Static introspection
    {T::get_N()}->std::convertible_to<int>;
    {T::get_half_N()}->std::convertible_to<int>;

    // Core physics interface
    {model.FillIntermediary(state)}->std::same_as<void>;

    {model.EffectiveAreas()}->std::same_as<std::tuple<typename T::state_type, typename T::state_type>>;

    {model.ComputeFlow(deltaP)}->std::same_as<typename T::scalar_type>;

    {model.Enl(state)}->std::same_as<typename T::scalar_type>;
    {model.Enl(state, recompute)}->std::same_as<typename T::scalar_type>;

    {model.Fnl(state)}->std::same_as<typename T::state_type>;
    {model.Fnl(state, recompute)}->std::same_as<typename T::state_type>;

    {model.KOp(state)}->std::same_as<typename T::state_type>;
    {model.ROp(state)}->std::same_as<typename T::state_type>;
    {model.MinvOp(state)}->std::same_as<typename T::state_type>;
};

} // namespace tarte