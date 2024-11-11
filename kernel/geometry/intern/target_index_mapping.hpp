// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/geometry/intern/congruency_mapping.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Target-Index Mapping class template
       *
       * \todo detailed description
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        int face_dim_>
      class TargetIndexMapping
      {
        static_assert(face_dim_ < Shape_::dimension, "invalid cell dimension");
        static_assert(face_dim_ >= 0, "invalid face dimension");

      protected:
        template<
          typename Inner_,
          typename Outer_>
        class CompIndexMap
        {
        protected:
          const Inner_& _inner;
          const Outer_& _outer;

        public:
          explicit CompIndexMap(
            const Inner_& inner,
            const Outer_& outer)
              :
            _inner(inner),
            _outer(outer)
          {
          }

          Index operator[](int i) const
          {
            return _outer[_inner[i]];
          }
        }; // class CompIndexMap

        /// cell orientation code
        int _cell_orient;

      public:
        template<
          typename TrgVerts_,
          typename SrcVerts_,
          typename VertIdx_>
        explicit TargetIndexMapping(
          const TrgVerts_& target_verts,
          const SrcVerts_& source_verts,
          const VertIdx_& vert_idx)
        {
          typedef CongruencySampler<Shape_> SamplerType;

          CompIndexMap<SrcVerts_, VertIdx_> comp_index_map(source_verts, vert_idx);

          _cell_orient = SamplerType::compare(comp_index_map, target_verts);
        }

        Index map(int face) const
        {
          return Index(CongruencyMapping<Shape_, face_dim_>::map(_cell_orient, face));
        }
      }; // class TargetIndexMapping
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
