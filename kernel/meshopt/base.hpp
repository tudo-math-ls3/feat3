// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>

namespace FEAT
{
  /**
   * \brief Namespace for everything mesh optimizer related
   *
   * Mesh optimizers in general need parts of Geometry (i.e. meshes), Trafo, Space (because FE knowledge is
   * required), Assembly to assemble systems of equations, and LAFEM to solve these equations.
   *
   * If possible, access them through their respective control classes.
   *
   */
  namespace Meshopt
  {
    /// \cond internal
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
      /// \endcond
    } // namespace Intern
  } // namespace Meshopt
} // namespace FEAT
