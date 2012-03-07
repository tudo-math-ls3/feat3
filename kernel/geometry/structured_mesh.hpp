#pragma once
#ifndef KERNEL_GEOMETRY_STRUCTURED_MESH_HPP
#define KERNEL_GEOMETRY_STRUCTURED_MESH_HPP 1

// includes, FEAST
#include <kernel/geometry/vertex_set.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Standard structured mesh policy.
     *
     * This class defines a default policy for the StructuredMesh class template.
     *
     * \tparam shape_dim_
     * The dimension of the shape (Hypercube) to be used for the mesh.
     *
     * \tparam VertexSet_
     * The vertex set class to be used by the mesh. By default, VertexSetFixed is used.
     *
     * \author Peter Zajac
     */
    template<
      int shape_dim_,
      typename VertexSet_ = VertexSetFixed<shape_dim_> >
    struct StructuredMeshPolicy
    {
      /// shape type; must always be a Hypercube shape
      typedef Shape::Hypercube<shape_dim_> ShapeType;

      /// Vertex set traits type
      typedef VertexSet_ VertexSetType;
    }; // struct StructuredMeshPolicy<...>

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

        static Index nverts(const Index num_slices[])
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
          num_entities[0] =
            num_slices[0] * (num_slices[1] + 1) * (num_slices[2] + 1) +
            num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) +
            num_slices[2] * (num_slices[0] + 1) * (num_slices[1] + 1);
          num_entities[0] =
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
     * \author Peter Zajac
     */
    template<typename Policy_>
    class StructuredMesh
    {
    public:
      /// policy type
      typedef Policy_ PolicyType;

      /// Shape type
      typedef typename PolicyType::ShapeType ShapeType;

      /// vertex set type
      typedef typename PolicyType::VertexSetType VertexSetType;

      /// dummy enum
      enum
      {
        /// shape dimension
        shape_dim = ShapeType::dimension,
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
       *
       * \param[in] num_coords
       * The number of coordinates per vertex. This parameter is passed to the constructor of the vertex set.
       *
       * \param[in] vertex_stride
       * The vertex stride. This parameter is passed to the constructor of the vertex set.
       */
      explicit StructuredMesh(
        const Index num_slices[],
        int num_coords = shape_dim,
        int vertex_stride = 0)
          :
        _vertex_set(Intern::StructCalcNumEntities<shape_dim>::num_verts(num_slices), num_coords, vertex_stride)
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
       * The nubmer of entities of dimension \p dim.
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

      /// Returns the name of the class.
      static String name()
      {
        return "StructuredMesh<...>";
      }
    }; // class StructuredMesh<...>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_STRUCTURED_MESH_HPP
