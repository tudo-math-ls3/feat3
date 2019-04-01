// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_PATCH_MESH_FACTORY_HPP
#define KERNEL_GEOMETRY_PATCH_MESH_FACTORY_HPP 1

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/patch_index_mapping.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \todo Documentation
    template<typename Mesh_>
    class PatchMeshFactory;

    /// \cond internal
    template<
      typename Shape_,
      int num_coords_,
      typename Coord_>
    class PatchMeshFactory<ConformalMesh<Shape_, num_coords_, Coord_> > :
      public Factory<ConformalMesh<Shape_, num_coords_, Coord_> >
    {
    public:
      /// mesh typedef
      typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
      typedef MeshPart<MeshType> MeshPartType;

      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    protected:
      const MeshType& _base_mesh;
      const MeshPartType& _patch_part;

    public:
      explicit PatchMeshFactory(const MeshType& base_mesh, const MeshPartType& patch_part) :
        _base_mesh(base_mesh),
        _patch_part(patch_part)
      {
      }

      virtual ~PatchMeshFactory()
      {
      }

      virtual Index get_num_entities(int dim) override
      {
        return _patch_part.get_num_entities(dim);
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        // fetch vertex-target-indices of the patch mesh
        const typename MeshPartType::template TargetSet<0>::Type& idx(_patch_part.template get_target_set<0>());
        typedef typename VertexSetType::VertexType VertexType;

        // fetch base-mesh vertex set
        const VertexSetType& vertex_set_in(_base_mesh.get_vertex_set());

        // loop over all vertices
        const Index num_verts(idx.get_num_entities());
        for(Index i(0); i < num_verts; ++i)
        {
          VertexType& vo(vertex_set[i]);
          const VertexType& vi(vertex_set_in[idx[i]]);
          for(int j(0); j < num_coords_; ++j)
          {
            vo[j] = vi[j];
          }
        }
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
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
          _patch_part.get_target_set_holder(),
          num_entities);
      }
    }; // class PatchMeshFactory<ConformalMesh<Shape_, num_coords_, Coord_> >
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PATCH_MESH_FACTORY_HPP
