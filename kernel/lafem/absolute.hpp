#pragma once
#ifndef KERNEL_LAFEM_ABSOLUTE_HPP
#define KERNEL_LAFEM_ABSOLUTE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/container.hpp>

#include <cmath>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Absolute value calculations.
     *
     * \tparam DT_ The datatype to be used.
     *
     * This class is used to generate absolute values.
     *
     * \author Dirk Ribbrock
     */
    template<typename DT_>
    struct Absolute
    {
        /**
         * \brief Calculate absolute value.
         *
         * \param[in] val The value to be processed
         *
         * \returns The absolute value of the input.
         */
      static DT_ value(DT_ val)
      {
        return Math::abs(val);
      }
    };

    /// \cond
    template<>
    struct Absolute<Index>
    {
      static Index value(Index val)
      {
        return val;
      }

    };
    /// \endcond

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ABSOLUTE_HPP
