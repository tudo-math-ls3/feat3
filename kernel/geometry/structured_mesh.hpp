// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_STRUCTURED_MESH_HPP
#define KERNEL_GEOMETRY_STRUCTURED_MESH_HPP 1

// includes, FEAT
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/struct_index_set.hpp>
#include <kernel/geometry/intern/structured_vertex_refiner.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Structured mesh class template
     *
     * \todo detailed documentation
     * \todo define index set type
     * \todo Implement neighbors lookup
     *
     * \tparam shape_dim_
     * The dimension of the shape (Hypercube) to be used for the mesh.
     *
     * \author Peter Zajac
     */
    template<
      int shape_dim_,
      int num_coords_ = shape_dim_,
      typename Coord_ = Real>
    class StructuredMesh
    {
      static_assert(shape_dim_ > 0, "invalid shape dimension");
      static_assert(num_coords_ >= shape_dim_, "invalid number of coordinates");

    public:
      /// Shape type
      typedef Shape::Hypercube<shape_dim_> ShapeType;

      /// Coordinate type
      typedef Coord_ CoordType;

      /// vertex set type
      typedef VertexSet<num_coords_, Coord_> VertexSetType;

      /// index set holder type
      typedef StructIndexSetHolder<shape_dim_> IndexSetHolderType;

      /// shape dimension
      static constexpr int shape_dim = shape_dim_;
      /// world dimension
      static constexpr int world_dim = num_coords_;

      /// the mesh is structured
      static constexpr bool is_structured = true;

      /**
       * \brief Index set type class template
       *
       * This nested class template is used to define the return type of the StructuredMesh::get_index_set()
       * function template.
       *
       * \tparam cell_dim_, face_dim_
       * The cell and face dimension parameters as passed to the StructuredMesh::get_index_set() function template.
       */
      template<
        int cell_dim_,
        int face_dim_>
      struct IndexSet
      {
        static_assert(cell_dim_ <= shape_dim, "invalid cell dimension");
        static_assert(face_dim_ < cell_dim_, "invalid face/cell dimension");
        static_assert(face_dim_ >= 0, "invalid face dimension");

        /// index set type
        typedef StructIndexSet<shape_dim_, cell_dim_, face_dim_> Type;
      }; // struct IndexSet<...>

    protected:
      /// number of slices for each direction
      Index _num_slices[shape_dim];
      /// number of entities for each dimension
      Index _num_entities[shape_dim + 1];

      /// the vertex set of the mesh.
      VertexSetType _vertex_set;

      /// index set holder
      IndexSetHolderType _index_set_holder;

    private:
      StructuredMesh& operator=(const StructuredMesh&);

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_slices
       * An array of length at least #shape_dim holding the number of slices for each direction.
       * Must not be \c nullptr.
       */
      explicit StructuredMesh(const Index num_slices[]) :
        _vertex_set(Intern::StructNumEntities<shape_dim_, 0>::compute(num_slices)),
        _index_set_holder(num_slices)
      {
        XASSERT(num_slices != nullptr);

        // store slice counts
        for(int i(0); i < shape_dim; ++i)
        {
          XASSERT(num_slices[i] > 0);
          _num_slices[i] = num_slices[i];
        }

        // calculate number of entities
        Intern::StructNumEntitiesWrapper<shape_dim_>::compute(_num_entities, _num_slices);
      }

      /**
       * \brief Factory constructor
       *
       * \param[in] factory
       * The factory that is to be used to create the mesh.
       */
      explicit StructuredMesh(Factory<StructuredMesh>& factory) :
        _vertex_set(
          Intern::StructNumEntities<shape_dim_, 0>::compute(
            Intern::NumSlicesWrapper<shape_dim>(factory).num_slices)),
        _index_set_holder(Intern::NumSlicesWrapper<shape_dim>(factory).num_slices)
      {
        // store slice count
        Intern::NumSlicesWrapper<shape_dim>::apply(factory, _num_slices);

        // calculate number of entities
        Intern::StructNumEntitiesWrapper<shape_dim_>::compute(_num_entities, _num_slices);

        // fill vertex set
        factory.fill_vertex_set(_vertex_set);
      }

      /**
       * \brief Copy Constructor
       *
       * \param[in] other
       * The structured mesh that is to be copied.
       */
      StructuredMesh(const StructuredMesh& other) :
        _vertex_set(other.get_vertex_set()),
        _index_set_holder(other.get_num_slices())
      {
        for(int i(0); i < shape_dim_; ++i)
        {
          _num_slices[i] = other.get_num_slices(i);
        }

        // calculate number of entities
        Intern::StructNumEntitiesWrapper<shape_dim_>::compute(_num_entities, _num_slices);
      }

      /// virtual destructor
      virtual ~StructuredMesh()
      {
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        return _vertex_set.bytes() + _index_set_holder.bytes();
      }

      /**
       * \brief Returns the number of slices.
       *
       * \param[in] dir
       * The direction of the slice whose count is to be returned. Must be 0 <= \p dir < #shape_dim.
       *
       * \returns
       * The number of slices in direction \p dir.
       */
      Index get_num_slices(int dir) const
      {
        XASSERT(dir >= 0);
        XASSERT(dir < shape_dim);
        return _num_slices[dir];
      }

      /// \brief Returns the num_slices array.
      const Index* get_num_slices() const
      {
        return _num_slices;
      }

      /**
       * \brief Returns the number of entities.
       *
       * \param[in] dim
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= #shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      Index get_num_entities(int dim) const
      {
        XASSERT(dim >= 0);
        XASSERT(dim <= shape_dim);
        return _num_entities[dim];
      }

      /**
       * \brief Returns the number of vertices.
       * \returns The number of vertices.
       */
      Index get_num_vertices() const
      {
        return _num_entities[0];
      }

      /**
       * \brief Returns the number of elements.
       * \returns The number of elements.
       */
      Index get_num_elements() const
      {
        return _num_entities[shape_dim];
      }

      /// Returns a reference to the vertex set of the mesh.
      VertexSetType& get_vertex_set()
      {
        return _vertex_set;
      }

      /** \copydoc get_vertex_set() */
      const VertexSetType& get_vertex_set() const
      {
        return _vertex_set;
      }

      /**
       * \brief Returns a reference to an index set.
       *
       * \tparam cell_dim_
       * The dimension of the entity whose index set is to be returned.
       *
       * \tparam face_dim_
       * The dimension of the face that the index set refers to.
       *
       * \returns
       * A reference to the index set.
       */
      template<
        int cell_dim_,
        int face_dim_>
      const StructIndexSet<shape_dim_, cell_dim_, face_dim_>& get_index_set() const
      {
        return _index_set_holder.template get_index_set<cell_dim_, face_dim_>();
      }

      /// \cond internal
      IndexSetHolderType& get_index_set_holder()
      {
        return _index_set_holder;
      }

      const IndexSetHolderType& get_index_set_holder() const
      {
        return _index_set_holder;
      }

      IndexSetHolderType* get_topology()
      {
        return &_index_set_holder;
      }

      const IndexSetHolderType* get_topology() const
      {
        return &_index_set_holder;
      }
      /// \endcond

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "StructuredMesh<...>";
      }
    }; // class StructuredMesh<...>

    /* ************************************************************************************************************* */

    /**
     * \brief Factory specialization for StructuredMesh class template
     *
     * \author Peter Zajac
     */
    template<
      int shape_dim_,
      int num_coords_,
      typename Coord_>
    class Factory< StructuredMesh<shape_dim_, num_coords_, Coord_> >
    {
    public:
      /// mesh type
      typedef StructuredMesh<shape_dim_, num_coords_, Coord_> MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;

    public:
      /// virtual destructor
      virtual ~Factory()
      {
      }

      /**
       * \brief Returns the number of slices
       *
       * \param[in] dir
       * The direction of the slice whose count is to be returned.
       *
       * \returns
       * The number of slices in direction \p dir.
       */
      virtual Index get_num_slices(int dir) = 0;

      /**
       * \brief Fills the vertex set.
       *
       * \param[in,out] vertex_set
       * The vertex set whose coordinates are to be filled.
       */
      virtual void fill_vertex_set(VertexSetType& vertex_set) = 0;
    }; // class Factory<StructuredMesh<...>>

    /* ************************************************************************************************************* */

    /**
     * \brief StandardRefinery implementation for StructuredMesh
     *
     * \author Peter Zajac
     */
    template<
      int shape_dim_,
      int num_coords_,
      typename Coord_>
    class StandardRefinery<StructuredMesh<shape_dim_, num_coords_, Coord_> > :
      public Factory< StructuredMesh<shape_dim_, num_coords_, Coord_> >
    {
    public:
      /// mesh type
      typedef StructuredMesh<shape_dim_, num_coords_, Coord_> MeshType;
      /// shape type
      typedef typename MeshType::ShapeType ShapeType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;

      /// shape dimension
      static constexpr int shape_dim = ShapeType::dimension;

    protected:
      /// coarse mesh reference
      const MeshType& _coarse_mesh;
      /// number of slices in coarse mesh
      Index _num_slices_coarse[shape_dim];

    public:
      /**
       * \brief Constructor
       *
       * \param[in] coarse_mesh
       * A reference to the coarse mesh that is to be refined.
       */
      explicit StandardRefinery(const MeshType& coarse_mesh) :
        _coarse_mesh(coarse_mesh)
      {
        for(int i(0); i < shape_dim; ++i)
        {
          _num_slices_coarse[i] = coarse_mesh.get_num_slices(i);
        }
      }

      /**
       * \brief Returns the number of slices
       *
       * \param[in] dir
       * The direction of the slice whose count is to be returned.
       *
       * \returns
       * The number of slices in direction \p dir.
       */
      virtual Index get_num_slices(int dir) override
      {
        return 2 * _coarse_mesh.get_num_slices(dir);
      }

      /**
       * \brief Fills the vertex set.
       *
       * \param[in,out] vertex_set
       * The vertex set whose coordinates are to be filled.
       */
      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        // refine vertices
        Intern::StructuredVertexRefiner<ShapeType, VertexSetType>
          ::refine(vertex_set, _coarse_mesh.get_vertex_set(), _num_slices_coarse);
      }
    }; // class StandardRefinery<StructuredMesh<...>>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_STRUCTURED_MESH_HPP
