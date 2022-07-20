// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_FACET_FLIPPER_HPP
#define KERNEL_GEOMETRY_FACET_FLIPPER_HPP 1

// includes, FEAT
#include <kernel/util/assertion.hpp>
#include <kernel/geometry/intern/congruency_sampler.hpp>
#include <kernel/geometry/intern/congruency_mapping.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/geometry/index_set.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      class FacetCongruencySampler
      {
      public:
        static constexpr int facet_dim = Shape_::dimension-1;
        typedef typename Shape::FaceTraits<Shape_, facet_dim>::ShapeType FacetType;

        template<typename Outer_>
        class CompIndexMap
        {
        protected:
          const Outer_& _outer;
          int _idx;

        public:
          explicit CompIndexMap(const Outer_& outer, int idx) : _outer(outer), _idx(idx) {}

          Index operator[](int i) const
          {
            return _outer[Geometry::Intern::FaceIndexMapping<Shape_, facet_dim, 0>::map(_idx,i)];
          }
        }; // class CompIndexMap

        template<typename VertsAtFace_, typename VertsAtShape_, typename FacesAtShape_>
        static int compare(const VertsAtFace_& verts_at_face, const VertsAtShape_& verts_at_shape,
          const FacesAtShape_& faces_at_shape, const Index shape_idx, const int local_face)
        {
          CompIndexMap<typename VertsAtShape_::IndexTupleType> cim(verts_at_shape[shape_idx], local_face);
          // get code and orientation
          int code = CongruencySampler<FacetType>::compare(verts_at_face[faces_at_shape(shape_idx, local_face)], cim);
          int loc_ori = CongruencySampler<FacetType>::orientation(code);
          int ref_ori = Shape::ReferenceCell<Shape_>::facet_orientation(local_face);
          return loc_ori * ref_ori;
        }
      }; // class FacetCongruencySampler<...>

      template<typename Facet_, int codim_ = Facet_::dimension>
      class FacetFlipRecurser
      {
      public:
        static constexpr int subdim = Facet_::dimension - codim_;

        static void flip(IndexSetWrapper<Facet_, Facet_::dimension-1>& isw, Index ifacet)
        {
          FacetFlipRecurser<Facet_, codim_-1>::flip(isw, ifacet);
          auto& idx_set = isw.template get_index_set<subdim>();
          Intern::CongruencyMapping<Facet_, subdim>::flip(idx_set[ifacet]);
        }
      };

      template<typename Facet_>
      class FacetFlipRecurser<Facet_, 0>
      {
      public:
        template<typename ISW_>
        static void flip(ISW_&, Index)
        {
          // nothing to do here
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Helper class for reorienting boundary facets to be positive
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class FacetFlipper
    {
    public:
      static_assert(Shape_::dimension > 1, "invalid shape dimension");

      /**
       * \brief Reorients the boundary facets in the index set holder to be positive
       */
      static void reorient(IndexSetHolder<Shape_>& ish)
      {
        static constexpr int shape_dim = Shape_::dimension;
        typedef typename Shape::FaceTraits<Shape_, shape_dim-1>::ShapeType FacetType;

        // get the three index sets:
        auto& facet_isw = ish.template get_index_set_wrapper<shape_dim-1>();
        const auto& verts_at_face = ish.template get_index_set<shape_dim-1, 0>();
        const auto& faces_at_elem = ish.template get_index_set<shape_dim, shape_dim-1>();
        const auto& verts_at_elem = ish.template get_index_set<shape_dim, 0>();

        // get the counts
        const Index num_faces = verts_at_face.get_num_entities();
        const Index num_elems = verts_at_elem.get_num_entities();

        // loop over all elements and determine all boundary facets
        std::vector<int> nfaces(num_faces, 0);
        for(Index ielem(0); ielem < num_elems; ++ielem)
        {
          for(int loc_face(0); loc_face < faces_at_elem.num_indices; ++loc_face)
          {
            ++nfaces.at(faces_at_elem(ielem, loc_face));
          }
        }

        // loop over all elements again
        for(Index ielem(0); ielem < num_elems; ++ielem)
        {
          // loop over all local facets of this element
          for(int loc_face(0); loc_face < faces_at_elem.num_indices; ++loc_face)
          {
            // only consider boundary faces = faces with 1 adjacent element
            const Index iface = faces_at_elem(ielem, loc_face);
            if(nfaces.at(iface) > 1)
              continue;

            // get the local facet's orientation
            int ori = Intern::FacetCongruencySampler<Shape_>::compare(
              verts_at_face, verts_at_elem, faces_at_elem, ielem, loc_face);

            // if the facet is oriented negatively, we have to flip it
            if(ori == 1)
              continue;

            // orientation must be +1 or -1
            XASSERTM(ori == -1, "invalid facet orientation");

            // okay, we have to flip all indices of this facet
            Intern::FacetFlipRecurser<FacetType>::flip(facet_isw, iface);
          }
        }
      }
    }; // class FacetFlipper
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_FACET_FLIPPER_HPP
