#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_PATCH_INDEX_MAPPING_HPP
#define KERNEL_GEOMETRY_INTERN_PATCH_INDEX_MAPPING_HPP 1

// includes, FEAT
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/target_set.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      struct PatchIndexMappingHelper
      {
        template<int num_indices_>
        static void apply(
          IndexSet<num_indices_>& iso,
          const IndexSet<num_indices_>& isi,
          const TargetSet& tsc,
          Index* tsf)
        {
          const Index num_cells(tsc.get_num_entities());
          for(Index i(0); i < num_cells; ++i)
          {
            const Index isx(tsc[i]);
            for(Index j(0); j < Index(num_indices_); ++j)
            {
              iso(i,j) = tsf[isi(isx, j)];
            }
          }
        }
      };

      template<
        typename Shape_,
        int face_dim_ = Shape_::dimension - 1>
      struct PatchIndexMappingWrapper
      {
        static void apply(
          IndexSetWrapper<Shape_, face_dim_>& iso,
          const IndexSetWrapper<Shape_, face_dim_>& isi,
          const TargetSetHolder<Shape_>& tsc,
          Index** tsf)
        {
          // recurse down
          PatchIndexMappingWrapper<Shape_, face_dim_ - 1>::apply(iso, isi, tsc, tsf);

          // call helper
          PatchIndexMappingHelper::apply(
            iso.template get_index_set<face_dim_>(),
            isi.template get_index_set<face_dim_>(),
            tsc.template get_target_set<Shape_::dimension>(),
            tsf[face_dim_]);
        }
      };

      template<typename Shape_>
      struct PatchIndexMappingWrapper<Shape_, 0>
      {
        static void apply(
          IndexSetWrapper<Shape_, 0>& iso,
          const IndexSetWrapper<Shape_, 0>& isi,
          const TargetSetHolder<Shape_>& tsc,
          Index** tsf)
        {
          // call helper
          PatchIndexMappingHelper::apply(
            iso.template get_index_set<0>(),
            isi.template get_index_set<0>(),
            tsc.template get_target_set<Shape_::dimension>(),
            tsf[0]);
        }
      };

      template<typename Shape_>
      struct PatchIndexMapping
      {
        static constexpr int shape_dim = Shape_::dimension;
        typedef typename Shape::FaceTraits<Shape_, shape_dim - 1>::ShapeType FaceType;

        static void apply(
          IndexSetHolder<Shape_>& iso,
          const IndexSetHolder<Shape_>& isi,
          const TargetSetHolder<Shape_>& tsc,
          const Index num_entities[])
        {
          Index *tsf[shape_dim + 1];
          _build_tsf(tsf, tsc, num_entities);
          _apply(iso, isi, tsc, &tsf[0]);
          for(int i(0); i <= shape_dim; ++i)
          {
            delete [] (tsf[i]);
          }
        }

      //protected:
        static void _build_tsf(
          Index** tsf,
          const TargetSetHolder<Shape_>& tsh,
          const Index num_entities[])
        {
          // recurse down
          PatchIndexMapping<FaceType>::_build_tsf(tsf, tsh, num_entities);

          tsf[shape_dim] = new Index[num_entities[shape_dim]];
          Index* tsi(tsf[shape_dim]);
          const TargetSet& ts(tsh.template get_target_set<shape_dim>());
          const Index ni(tsh.get_num_entities(shape_dim));
          for(Index i(0); i < ni; ++i)
          {
            tsi[ts[i]] = i;
          }
        }

        static void _apply(
          IndexSetHolder<Shape_>& iso,
          const IndexSetHolder<Shape_>& isi,
          const TargetSetHolder<Shape_>& tsc,
          Index** tsf)
        {
          // recurse down
          PatchIndexMapping<FaceType>::_apply(iso, isi, tsc, tsf);

          // call wrapper
          PatchIndexMappingWrapper<Shape_>::apply(
            iso.template get_index_set_wrapper<shape_dim>(),
            isi.template get_index_set_wrapper<shape_dim>(), tsc, tsf);
        }
      };

      template<>
      struct PatchIndexMapping<Shape::Vertex>
      {
      //protected:
        static void _build_tsf(
          Index** tsf,
          const TargetSetHolder<Shape::Vertex>& tsh,
          const Index num_entities[])
        {
          tsf[0] = new Index[num_entities[0]];
          Index* tsi(tsf[0]);
          const TargetSet& ts(tsh.get_target_set<0>());
          const Index ni(tsh.get_num_entities(0));
          for(Index i(0); i < ni; ++i)
          {
            tsi[ts[i]] = i;
          }
        }

        template<typename ISO_, typename ISI_, typename TSC_>
        static void _apply(ISO_&, const ISI_&, const TSC_&, Index**)
        {
          // do nothing
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_PATCH_INDEX_MAPPING_HPP
