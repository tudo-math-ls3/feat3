#pragma once
#ifndef KERNEL_SPACE_RANNACHER_TUREK_VARIANT_HPP
#define KERNEL_SPACE_RANNACHER_TUREK_VARIANT_HPP 1

// includes, FEAST
#include <kernel/space/base.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace RannacherTurek
    {
      /**
       * \brief Variant namespace for Rannacher-Turek elements
       */
      namespace Variant
      {
        /**
         * \brief Standard Non-Parametric variant
         */
        struct StdNonPar
        {
          /*enum
          {
            non_par = 1,
            bubble = 0
          };*/

          static String name()
          {
            return "StdNonPar";
          }
        };
      } // namespace Variant
    } // namespace RannacherTurek
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_RANNACHER_TUREK_VARIANT_HPP
