#pragma once
#ifndef KERNEL_GEOMETRY_STRUCTURED_MESH_HPP
#define KERNEL_GEOMETRY_STRUCTURED_MESH_HPP 1

// includes, FEAST
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/intern/structured_vertex_refiner.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      // helper class to calculate the number of entities from the number of slices
      template<int dim_>
      struct StructCalcNumEntities;

      template<>
      struct StructCalcNumEntities<1>
      {
        static void apply(Index num_entities[], const Index num_slices[])
        {
          num_entities[0] = num_slices[0] + 1;
          num_entities[1] = num_slices[0];
        }

        static Index num_verts(const Index num_slices[])
        {
          return num_slices[0] + 1;
        }
      };

      template<>
      struct StructCalcNumEntities<2>
      {
        static void apply(Index num_entities[], const Index num_slices[])
        {
          num_entities[0] = (num_slices[0] + 1) * (num_slices[1] + 1);
          num_entities[1] = (num_slices[0] + 1) * num_slices[1] + num_slices[0] * (num_slices[1] + 1);
          num_entities[2] = num_slices[0] * num_slices[1];
        }

        static Index num_verts(const Index num_slices[])
        {
          return (num_slices[0] + 1) * (num_slices[1] + 1);
        }
      };

      template<>
      struct StructCalcNumEntities<3>
      {
        static void apply(Index num_entities[], const Index num_slices[])
        {
          num_entities[0] = (num_slices[0] + 1) * (num_slices[1] + 1) * (num_slices[2] + 1);
          num_entities[1] =
            num_slices[0] * (num_slices[1] + 1) * (num_slices[2] + 1) +
            num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) +
            num_slices[2] * (num_slices[0] + 1) * (num_slices[1] + 1);
          num_entities[2] =
            (num_slices[0] + 1) * num_slices[1] * num_slices[2] +
            (num_slices[1] + 1) * num_slices[0] * num_slices[2] +
            (num_slices[2] + 1) * num_slices[0] * num_slices[1];
          num_entities[3] = num_slices[0] * num_slices[1] * num_slices[2];
        }

        static Index num_verts(const Index num_slices[])
        {
          return (num_slices[0] + 1) * (num_slices[1] + 1) * (num_slices[2] + 1);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Structured mesh class template
     *
     * \todo detailed documentation
     * \todo define index set type
     *
     * \tparam shape_dim_
     * The dimension of the shape (Hypercube) to be used for the mesh.
     *
     * \author Peter Zajac
     */
    template<
      int shape_dim_,
      int num_coords_ = shape_dim_,
      int stride_ = num_coords_,
      typename Coord_ = Real>
    class StructuredMesh
    {
      static_assert(shape_dim_ > 0, "invalid shape dimension");
      static_assert(num_coords_ >= shape_dim_, "invalid number of coordinates");
      static_assert(stride_ >= num_coords_, "invalid stride");

    public:
      /// Shape type
      typedef Shape::Hypercube<shape_dim_> ShapeType;

      /// vertex set type
      typedef VertexSetFixed<num_coords_, stride_, Coord_> VertexSetType;

      /// dummy enum
      enum
      {
        /// shape dimension
        shape_dim = shape_dim_
      };

    protected:
      /// number of slices for each direction
      Index _num_slices[shape_dim];
      /// number of entities for each dimension
      Index _num_entities[shape_dim + 1];

      /// the vertex set of the mesh.
      VertexSetType _vertex_set;

    private:
      StructuredMesh(const StructuredMesh&);
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
        _vertex_set(Intern::StructCalcNumEntities<shape_dim>::num_verts(num_slices))
      {
        CONTEXT(name() + "::StructuredMesh()");
        ASSERT_(num_slices != nullptr);

        // store slice counts
        for(int i(0); i < shape_dim; ++i)
        {
          ASSERT_(num_slices[i] > 0);
          _num_slices[i] = num_slices[i];
        }

        // calculate number of enitites
        Intern::StructCalcNumEntities<shape_dim>::apply(_num_entities, _num_slices);
      }

      /**
       * \brief Factory constructor
       *
       * \param[in] factory
       * The factory that is to be used to create the mesh.
       */
      explicit StructuredMesh(const Factory<StructuredMesh>& factory) :
        _vertex_set(
          Intern::StructCalcNumEntities<shape_dim>::num_verts(
            Intern::NumSlicesWrapper<shape_dim>(factory).num_slices))
      {
        CONTEXT(name() + "::StructuredMesh() [factory]");

        // store slice count
        Intern::NumSlicesWrapper<shape_dim>::apply(factory, _num_slices);

        // calculate number of enitites
        Intern::StructCalcNumEntities<shape_dim>::apply(_num_entities, _num_slices);

        // fill vertex set
        factory.fill_vertex_set(_vertex_set);
      }

      /// virtual destructor
      virtual ~StructuredMesh()
      {
        CONTEXT(name() + "::~StructuredMesh()");
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
        CONTEXT(name() + "::get_num_slices()");
        ASSERT_(dir >= 0);
        ASSERT_(dir < shape_dim);
        return _num_slices[dir];
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
        CONTEXT(name() + "::get_num_entities()");
        ASSERT_(dim >= 0);
        ASSERT_(dim <= shape_dim);
        return _num_entities[dim];
      }

      /// Returns a reference to the vertex set of the mesh.
      VertexSetType& get_vertex_set()
      {
        CONTEXT(name() + "::get_vertex_set()");
        return _vertex_set;
      }

      /** \copydoc get_vertex_set() */
      const VertexSetType& get_vertex_set() const
      {
        CONTEXT(name() + "::get_vertex_set() [const]");
        return _vertex_set;
      }

      /**
       * \brief Refines the mesh.
       *
       * This function applies the standard refinement algorithm onto the mesh and returns the refined mesh.
       *
       * \returns
       * A pointer to the refined mesh.
       */
      StructuredMesh* refine() const
      {
        CONTEXT(name() + "::refine()");

        return new StructuredMesh(StandardRefinery<StructuredMesh>(*this));
      }

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
     * \brief Factory specialisation for StructuredMesh class template
     *
     * \author Peter Zajac
     */
    template<
      int shape_dim_,
      int num_coords_,
      int stride_,
      typename Coord_>
    class Factory< StructuredMesh<shape_dim_, num_coords_, stride_, Coord_> >
    {
    public:
      /// mesh type
      typedef StructuredMesh<shape_dim_, num_coords_, stride_, Coord_> MeshType;
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
      virtual Index get_num_slices(int dir) const = 0;

      /**
       * \brief Fills the vertex set.
       *
       * \param[in,out] vertex_set
       * The vertex set whose coordinates are to be filled.
       */
      virtual void fill_vertex_set(VertexSetType& vertex_set) const = 0;
    }; // class Factory<StructuredMesh<...>>

    /* ************************************************************************************************************* */

    /**
     * \brief StandardRefinery implementation for StructuredMesh
     */
    template<
      int shape_dim_,
      int num_coords_,
      int stride_,
      typename Coord_>
    class StandardRefinery<StructuredMesh<shape_dim_, num_coords_, stride_, Coord_>, Nil> :
      public Factory< StructuredMesh<shape_dim_, num_coords_, stride_, Coord_> >
    {
    public:
      /// mesh type
      typedef StructuredMesh<shape_dim_, num_coords_, stride_, Coord_> MeshType;
      /// shape type
      typedef typename MeshType::ShapeType ShapeType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;

      /// dummy enum
      enum
      {
        /// shape dimension
        shape_dim = ShapeType::dimension
      };

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
      virtual Index get_num_slices(int dir) const
      {
        return 2 * _coarse_mesh.get_num_slices(dir);
      }

      /**
       * \brief Fills the vertex set.
       *
       * \param[in,out] vertex_set
       * The vertex set whose coordinates are to be filled.
       */
      virtual void fill_vertex_set(VertexSetType& vertex_set) const
      {
        // refine vertices
        Intern::StructuredVertexRefiner<ShapeType, VertexSetType>
          ::refine(vertex_set, _coarse_mesh.get_vertex_set(), _num_slices_coarse);
      }
    }; // class StandardRefinery<StructuredMesh<...>>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_STRUCTURED_MESH_HPP