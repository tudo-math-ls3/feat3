// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/macro_index_mapping.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \todo Documentation
    template<typename Mesh_>
    class MacroFactory;

    /// \cond internal
    template<
      typename Shape_,
      int num_coords_,
      typename Coord_>
    class MacroFactory<ConformalMesh<Shape_, num_coords_, Coord_> > :
      public Factory<ConformalMesh<Shape_, num_coords_, Coord_> >
    {
    public:
      /// mesh typedef
      typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    protected:
      const MeshType& _base_mesh;
      const Index _cell_idx;

    public:
      explicit MacroFactory(const MeshType& base_mesh, Index cell_idx) :
        _base_mesh(base_mesh),
        _cell_idx(cell_idx)
      {
        XASSERTM(cell_idx < base_mesh.get_num_entities(Shape_::dimension), "cell index out-of-bounds");
      }

      virtual Index get_num_entities(int dim) override
      {
        return Index(Intern::DynamicNumFaces<Shape_>::value(dim));
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        // fetch the cell's vertex indices
        typedef typename MeshType::template IndexSet<Shape_::dimension, 0>::Type IndexSetType;
        const IndexSetType& idx(_base_mesh.template get_index_set<Shape_::dimension, 0>());

        // fetch base-mesh vertex set
        const VertexSetType& vertex_set_in(_base_mesh.get_vertex_set());

        // loop over all vertices
        for(Index i(0); i < Index(IndexSetType::num_indices); ++i)
        {
          vertex_set[i] = vertex_set_in[idx(_cell_idx, i)];
        }
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
      {
        Intern::MacroIndexWrapper<Shape_>::build(index_set_holder);
      }
    };

    template<typename BaseMesh_>
    class MacroFactory<MeshPart<BaseMesh_> > :
      public Factory<MeshPart<BaseMesh_> >
    {
    public:
      /// mesh typedef
      typedef MeshPart<BaseMesh_> MeshType;
      typedef typename MeshType::ShapeType ShapeType;

      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
      typedef typename MeshType::AttributeSetContainer AttributeSetContainer;

    protected:
      const BaseMesh_& _base_mesh;
      const Index _cell_idx;

    public:
      explicit MacroFactory(const BaseMesh_& base_mesh, Index cell_idx) :
        _base_mesh(base_mesh),
        _cell_idx(cell_idx)
      {
        XASSERTM(cell_idx < base_mesh.get_num_entities(ShapeType::dimension), "cell index out-of-bounds");
      }

      virtual Index get_num_entities(int dim) override
      {
        return Index(Intern::DynamicNumFaces<ShapeType>::value(dim));
      }

      virtual void fill_attribute_sets(AttributeSetContainer&) override
      {
        // nothing to do here
      }

      virtual void fill_index_sets(std::unique_ptr<IndexSetHolderType>&) override
      {
        // nothing to do here
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        // set cell index
        target_set_holder.template get_target_set<ShapeType::dimension>()[0] = _cell_idx;
        // fill remaining indices
        Intern::MacroTargetWrapper<ShapeType>::build(target_set_holder, _base_mesh.get_index_set_holder(), _cell_idx);
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
