#pragma once
#ifndef KERNEL_GEOMETRY_SHAPE_CONVERT_FACTORY_HPP
#define KERNEL_GEOMETRY_SHAPE_CONVERT_FACTORY_HPP 1

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/index_calculator.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/shape_convert_index.hpp>
#include <kernel/geometry/intern/shape_convert_target.hpp>
#include <kernel/geometry/intern/shape_convert_traits.hpp>
#include <kernel/geometry/intern/shape_convert_vertex.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Shape-Conversion Mesh Factory class template
     *
     * This class template acts as a conversion factory for conformal Simplex and Hypercube meshes,
     * i.e. this factory will convert a Simplex<n> mesh into a Hypercube<n> mesh or vice versa by
     * conforming element subdivision.
     *
     * \tparam Mesh_
     * The type that the input mesh is to be converted to.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
#ifndef DOXYGEN
    class ShapeConvertFactory;
#else
    class ShapeConvertFactory :
      public Factory<Mesh_>
    {
    public:
      /// the output mesh type
      typedef Mesh_ MeshType;
      /// the input mesh type
      typedef ... OtherMeshType;

      /**
       * \brief Constructor.
       *
       * \param[in] other_mesh
       * The input mesh that is to be converted.
       */
      explicit ShapeConvertFactory(const OtherMeshType& other_mesh);
    };
#endif // DOXYGEN

    /**
     * \brief ShapeConvertFactory implementation for ConformalMesh
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      int num_coords_,
      typename Coord_>
    class ShapeConvertFactory<ConformalMesh<Shape_, num_coords_, Coord_> > :
      public Factory<ConformalMesh<Shape_, num_coords_, Coord_> >
    {
    public:
      typedef Factory<ConformalMesh<Shape_, num_coords_, Coord_> > BaseClass;
      typedef Shape_ ShapeType;
      typedef typename Intern::OtherShape<Shape_>::Type OtherShapeType;
      typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
      typedef ConformalMesh<OtherShapeType, num_coords_, Coord_> OtherMeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      /// shape dimension
      static constexpr int shape_dim = ShapeType::dimension;

    protected:
      const OtherMeshType& _other_mesh;
      /// number of entities for coarse mesh
      Index _num_entities_in[shape_dim + 1];
      /// number of entities for fine mesh
      Index _num_entities_out[shape_dim + 1];

    public:
      explicit ShapeConvertFactory(const OtherMeshType& other_mesh) :
        _other_mesh(other_mesh)
      {
        // get number of entities in input mesh
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities_out[i] = _num_entities_in[i] = _other_mesh.get_num_entities(i);
        }

        // calculate number of entities in output mesh
        Intern::EntityCountWrapper<Intern::ShapeConvertTraits, OtherShapeType>::query(_num_entities_out);
      }

      virtual Index get_num_entities(int dim) override
      {
        return _num_entities_out[dim];
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        // call wrapper
        Intern::ShapeConvertVertexWrapper<OtherShapeType, VertexSetType>
          ::refine(vertex_set, _other_mesh.get_vertex_set(), _other_mesh.get_index_set_holder());
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
      {
        // call wrapper to build essential index sets
        Intern::ShapeConvertIndexWrapper<ShapeType>
          ::refine(index_set_holder, _num_entities_in, _other_mesh.get_index_set_holder());

        // build redundant index sets
        RedundantIndexSetBuilder<ShapeType>::compute(index_set_holder);
      }
    }; // class ShapeConvertFactory<ConformalMesh<...>>

    /**
     * \brief ShapeConvertFactory implementation for ConformalMesh
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename Coord_>
    class ShapeConvertFactory<MeshPart<ConformalMesh<Shape_, Shape_::dimension, Coord_> > > :
      public Factory<MeshPart<ConformalMesh<Shape_, Shape_::dimension, Coord_> > >
    {
    public:
      typedef Factory<MeshPart<ConformalMesh<Shape_, Shape_::dimension, Coord_> > > BaseClass;
      typedef Shape_ ShapeType;
      typedef typename Intern::OtherShape<Shape_>::Type OtherShapeType;
      typedef MeshPart<ConformalMesh<Shape_, Shape_::dimension, Coord_> > MeshType;
      typedef MeshPart<ConformalMesh<OtherShapeType, Shape_::dimension, Coord_> > OtherMeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index set holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

      /// shape dimension
      static constexpr int shape_dim = ShapeType::dimension;

    protected:
      const OtherMeshType& _other_mesh;
      /// number of entities for coarse mesh
      Index _num_entities_in[shape_dim + 1];
      /// number of entities for fine mesh
      Index _num_entities_out[shape_dim + 1];

    public:
      explicit ShapeConvertFactory(const OtherMeshType& other_mesh) :
        _other_mesh(other_mesh)
      {
        // get number of entities in input mesh
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities_out[i] = _num_entities_in[i] = _other_mesh.get_num_entities(i);
        }

        // calculate number of entities in output mesh
        Intern::EntityCountWrapper<Intern::ShapeConvertTraits, OtherShapeType>::query(_num_entities_out);

      }

      virtual Index get_num_entities(int dim)
      {
        return _num_entities_out[dim];
      }

      virtual int get_num_coords()
      {
        return _other_mesh.get_vertex_set().get_num_coords();
      }

      virtual int get_vertex_stride()
      {
        return _other_mesh.get_vertex_set().get_stride();
      }

      virtual void fill_vertex_set(VertexSetType& vertex_set)
      {
        // call wrapper
        Intern::ShapeConvertVertexWrapper<OtherShapeType, VertexSetType>
          ::refine(vertex_set, _other_mesh.get_vertex_set(), _other_mesh.get_index_set_holder());
      }

      virtual void fill_index_sets(IndexSetHolderType& index_set_holder)
      {
        // call wrapper to build essential index sets
        Intern::ShapeConvertIndexWrapper<ShapeType>
          ::refine(index_set_holder, _num_entities_in, _other_mesh.get_index_set_holder());

        // build redundant index sets
        RedundantIndexSetBuilder<ShapeType>::compute(index_set_holder);
      }

      virtual void fill_target_sets(TargetSetHolderType& /*target_set_holder*/)
      {
        // TODO
      }
    }; // class ShapeConvertFactory<MeshPart<ConformalMesh<...>>>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_SHAPE_CONVERT_FACTORY_HPP
