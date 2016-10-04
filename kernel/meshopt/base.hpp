#pragma once
#ifndef KERNEL_MESHOPT_BASE_HPP
#define KERNEL_MESHOPT_BASE_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>

/// \cond internal
namespace FEAT
{
  namespace Meshopt
  {
    namespace Intern
    {
      /**
       * \brief Information about the FE space the transformation belongs to
       *
       * \tparam Trafo_
       * The transformation
       */
      template<typename Trafo_>
      struct TrafoFE
      {
#ifdef DOXYGEN
        /// Finite Element space the trafo belongs to
        typedef ... Space;
#endif
      };

      template<typename Mesh_>
      struct TrafoFE<Trafo::Standard::Mapping<Mesh_>>
      {
        typedef Space::Lagrange1::Element<Trafo::Standard::Mapping<Mesh_>> Space;
      };
    } // namespace Intern
  } // namespace Meshopt
} // namespace FEAT
/// \endcond
#endif // KERNEL_MESHOPT_BASE_HPP
