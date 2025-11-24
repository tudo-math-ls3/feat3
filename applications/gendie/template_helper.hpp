#pragma once
#include <type_traits>

namespace Gendie
{
  namespace Intern
  {
    // check if fbm fields are supported
    template<typename, typename = void>
    constexpr bool supports_fbm = false;

    template<typename T>
    constexpr bool supports_fbm<T,
    std::void_t<decltype(T::fbm_support)>> = T::fbm_support;

    // check if coloring field exists
    template<typename, typename = void>
    constexpr bool supports_element_coloring = false;

    template<typename T>
    constexpr bool supports_element_coloring<T,
    std::void_t<decltype(T::LevelType::element_coloring)>> = true;

    // check if pressure scaling flag is supported
    template<typename, typename = void>
    constexpr bool pressure_scaling_on = false;

    template<typename T>
    constexpr bool pressure_scaling_on<T,
    std::void_t<decltype(T::scale_pressure)>> = T::scale_pressure;
  }

}