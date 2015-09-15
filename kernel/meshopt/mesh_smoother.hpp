#pragma once
#ifndef KERNEL_MESHOPT_MESH_SMOOTHER_HPP
#define KERNEL_MESHOPT_MESH_SMOOTHER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>

namespace FEAST
{
  /**
   * \brief Namespace for everything mesh optimiser / mesh smoother related
   *
   * Because mesh smoothers in general need part of Geometry (i.e. meshes), Trafo, Space (because FE knowledge is
   * required), Assembly to assemble systems of equations, and LAFEM to solve these equations.
   *
   */
  namespace Meshopt
  {

    /**
     * \brief Baseclass for mesh optimisation algorithms
     *
     * This abstract class is the baseclass for all mesh optimisation algorithms, which can be direct, variational
     * or something entirely different.
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \author Jordi Paul
     *
     */
    template<typename MeshType_>
    class MeshSmoother
    {
      public:
        /// Type of the mesh to optimise
        typedef MeshType_ MeshType;
        /// Our datatype
        typedef typename MeshType::CoordType CoordType;
        /// The shape type
        typedef typename MeshType::ShapeType ShapeType;
        /// Type for the vectors to hold coordinates etc.
        typedef LAFEM::DenseVectorBlocked<Mem::Main, CoordType, Index, MeshType::world_dim> VertexVectorType;

      public:
        /// The mesh for the underlying transformation
        MeshType& _mesh;
        /// Coordinates, used for setting new boundary values etc.
        VertexVectorType _coords;

      public:
        /// \brief Constructor
        explicit MeshSmoother(MeshType& mesh_) :
          _mesh(mesh_),
          _coords(mesh_.get_num_entities(0), CoordType(0))
          {
          }

        /// \brief Virtual destructor
        virtual ~MeshSmoother()
        {
        }

        /// \brief Initialises parts of the MeshSmoother not set in in the constructor
        virtual void init()
        {
          get_coords();
        }

        /// \brief Gets the coordinates from the underlying mesh and saves them in _coords.
        virtual void get_coords()
        {
          const typename MeshType::VertexSetType& vertex_set = _mesh.get_vertex_set();

          for(Index i(0); i < _mesh.get_num_entities(0); ++i)
            _coords(i, vertex_set[i]);
        }

        /// \brief Sets the coordinates in the underlying mesh to _coords.
        virtual void set_coords()
        {
          typename MeshType::VertexSetType& vertex_set = _mesh.get_vertex_set();

          for(Index i(0); i < _mesh.get_num_entities(0); ++i)
            vertex_set[i] = _coords(i);
        }

        /// \brief Optimises the mesh according to the criteria implemented in the mesh smoother.
        virtual void optimise() = 0;

        /// \brief Prepares the mesh optimiser for application
        virtual void prepare() = 0;

    }; // class MeshSmoother

    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Information about the FE space the transformation belongs to
       *
       * \tparam Trafo_
       * The transformation
       */
      template<typename Trafo_>
      struct TrafoFE
      {
#ifdef DOXYGEN
        /// Finite Element space the trafo belongs to
        typedef ... Space;
#endif
      };

      template<typename Mesh_>
      struct TrafoFE<Trafo::Standard::Mapping<Mesh_>>
      {
        typedef Space::Lagrange1::Element<Trafo::Standard::Mapping<Mesh_>> Space;
      };
    }
    /// \endcond

  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_MESH_SMOOTHER_HPP
