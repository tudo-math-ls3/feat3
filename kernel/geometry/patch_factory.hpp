#pragma once
#ifndef KERNEL_GEOMETRY_PATCH_FACTORY_HPP
#define KERNEL_GEOMETRY_PATCH_FACTORY_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/patch_index_mapping.hpp>

namespace FEAST
{
  namespace Geometry
  {
    template<typename Mesh_>
    class PatchFactory;

    template<
      typename Shape_,
      int num_coords_,
      int stride_,
      typename Coord_>
    class PatchFactory<ConformalMesh<Shape_, num_coords_, stride_, Coord_> > :
      public Factory<ConformalMesh<Shape_, num_coords_, stride_, Coord_> >
    {
    public:
      /// mesh typedef
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
      typedef MeshPart<MeshType> CellSetType;

      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    protected:
      const MeshType& _base_mesh;
      const CellSetType& _patch_set;

    public:
      explicit PatchFactory(const MeshType& base_mesh, const CellSetType& patch_set) :
        _base_mesh(base_mesh),
        _patch_set(patch_set)
      {
      }

      virtual ~PatchFactory()
      {
      }

      virtual Index get_num_entities(int dim)
      {
        return _patch_set.get_num_entities(dim);
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set)
      {
        // fetch vertex-target-indices of the patch mesh
        const typename CellSetType::template TargetSet<0>::Type& idx(_patch_set.template get_target_set<0>());
        typedef typename VertexSetType::VertexReference VertexRef;
        typedef typename VertexSetType::ConstVertexReference VertexConstRef;

        // fetch base-mesh vertex set
        const VertexSetType& vertex_set_in(_base_mesh.get_vertex_set());

        // loop over all vertices
        const Index num_verts(idx.get_num_entities());
        for(Index i(0); i < num_verts; ++i)
        {
          VertexRef vo(vertex_set[i]);
          VertexConstRef vi(vertex_set_in[idx[i]]);
          for(int j(0); j < num_coords_; ++j)
          {
            vo[j] = vi[j];
          }
        }
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
      {
        // set up base-mesh num_entities
        Index num_entities[Shape_::dimension + 1];
        for(int i(0); i <= Shape_::dimension; ++i)
        {
          num_entities[i] = _base_mesh.get_num_entities(i);
        }

        // build index sets
        Intern::PatchIndexMapping<Shape_>::apply(
          index_set_holder,
          _base_mesh.get_index_set_holder(),
          _patch_set.get_target_set_holder(),
          num_entities);
      }
    };
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_PATCH_FACTORY_HPP
