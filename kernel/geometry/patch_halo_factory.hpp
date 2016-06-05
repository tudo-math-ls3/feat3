#pragma once
#ifndef KERNEL_GEOMETRY_PATCH_HALO_FACTORY_HPP
#define KERNEL_GEOMETRY_PATCH_HALO_FACTORY_HPP 1

// includes, FEAT
#include <kernel/geometry/mesh_part.hpp>

// includes, system
#include <map>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct PatchHaloWrapper;
    }

    /// \endcond

    /// \todo Documentation
    template<typename Mesh_>
    class PatchHaloFactory DOXY({});

    /// \cond internal
    template<typename Shape_, int num_coords_, int stride_, typename Coord_>
    class PatchHaloFactory< MeshPart<ConformalMesh<Shape_, num_coords_, stride_, Coord_> > > :
      public Factory< MeshPart<ConformalMesh<Shape_, num_coords_, stride_, Coord_> > >
    {
    public:
      /// mesh type
      typedef MeshPart<ConformalMesh<Shape_, num_coords_, stride_, Coord_> > MeshType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    private:
      const MeshType& _patch_set;
      const MeshType& _halo_set;

      typedef std::map<Index, Index> IndexMap;

      IndexMap _ics[Shape_::dimension + 1];

    public:
      explicit PatchHaloFactory(const MeshType& patch_set, const MeshType& halo_set) :
        _patch_set(patch_set),
        _halo_set(halo_set)
      {
        Intern::PatchHaloWrapper<Shape_>::build(_ics, _patch_set.get_target_set_holder());
      }

      virtual Index get_num_entities(int dim) override
      {
        return _halo_set.get_num_entities(dim);
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        Intern::PatchHaloWrapper<Shape_>::map(_ics, _halo_set.get_target_set_holder(), target_set_holder);
      }

      virtual void fill_attribute_sets(typename MeshType::MeshAttributeContainer& DOXY(attribute_container)) override
      {
        // do nothing as the object has no attribute sets
      }

      virtual void fill_index_sets(typename MeshType::IndexSetHolderType*& DOXY(index_set_holder)) override
      {
        // do nothing as the object has no index sets
      }

    }; // class PatchHaloFactory<MeshPart<...>>
    /// \endcond

    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct PatchHaloHelper
      {
        static void build(std::map<Index, Index>& idx_map, const TargetSet& trg_set)
        {
          const Index num_entities(trg_set.get_num_entities());
          for(Index i(0); i < num_entities; ++i)
          {
            bool okay = idx_map.insert(std::make_pair(trg_set[i], i)).second;
#if defined DEBUG
            ASSERT(okay);
#else
            (void)okay;
#endif
          }
        }

        static void map(const std::map<Index, Index>& ics, const TargetSet& trg_in, TargetSet& trg_out)
        {
          const Index n(trg_in.get_num_entities());
          for(Index i(0); i < n; ++i)
          {
            std::map<Index, Index>::const_iterator it(ics.find(trg_in[i]));
            ASSERT(it != ics.end());
            trg_out[i] = it->second;
          }
        }
      };

      template<typename Shape_>
      struct PatchHaloWrapper
      {
        static void build(std::map<Index, Index>* idx_map, const TargetSetHolder<Shape_>& trg_holder)
        {
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType FacetType;
          PatchHaloWrapper<FacetType>::build(idx_map, trg_holder);
          PatchHaloHelper<Shape_>::build(idx_map[Shape_::dimension],
            trg_holder.template get_target_set<Shape_::dimension>());
        }

        static void map(std::map<Index, Index>* idx_map, const TargetSetHolder<Shape_>& trg_in, TargetSetHolder<Shape_>& trg_out)
        {
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType FacetType;
          PatchHaloWrapper<FacetType>::map(idx_map, trg_in, trg_out);
          PatchHaloHelper<Shape_>::map(idx_map[Shape_::dimension],
            trg_in.template get_target_set<Shape_::dimension>(),
            trg_out.template get_target_set<Shape_::dimension>());
        }
      };


      template<>
      struct PatchHaloWrapper<Shape::Vertex>
      {
        static void build(std::map<Index, Index>* idx_map, const TargetSetHolder<Shape::Vertex>& trg_holder)
        {
          PatchHaloHelper<Shape::Vertex>::build(idx_map[0], trg_holder.get_target_set<0>());
        }

        static void map(std::map<Index, Index>* idx_map, const TargetSetHolder<Shape::Vertex>& trg_in, TargetSetHolder<Shape::Vertex>& trg_out)
        {
          PatchHaloHelper<Shape::Vertex>::map(idx_map[0], trg_in.get_target_set<0>(), trg_out.get_target_set<0>());
        }
      };
    }
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PATCH_HALO_FACTORY_HPP
