// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Argyris
    {
      /**
       * \brief Argyris Dof-Traits class template.
       *
       * \author Peter Zajac
       */
      template<typename Shape_, int dim_>
      struct DofTraits
      {
        static constexpr int count = 0;
      };

      template<>
      struct DofTraits<Shape::Simplex<2>, 0>
      {
        /// 6 dofs per vertex
        static constexpr int count = 6;
      };

      template<>
      struct DofTraits<Shape::Simplex<2>, 1>
      {
        /// 1 dof per edge
        static constexpr int count = 1;
      };
    } // namespace Argyris
  } // namespace Space
} // namespace FEAT
