#pragma once
#ifndef KERNEL_MESHOPT_MESH_QUALITY_FUNCTIONAL_HPP
#define KERNEL_MESHOPT_MESH_QUALITY_FUNCTIONAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/util/mpi_cout.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/meshopt/base.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>

namespace FEAT
{
  /**
   * \brief Namespace for everything mesh optimiser related
   *
   * Mesh optimisers in general need parts of Geometry (i.e. meshes), Trafo, Space (because FE knowledge is
   * required), Assembly to assemble systems of equations, and LAFEM to solve these equations.
   *
   * If possible, access them through their respective control classes.
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
    class MeshQualityFunctional
    {
      public:
        /// Type of the mesh to optimise
        typedef MeshType_ MeshType;
        /// Our datatype
        typedef typename MeshType::CoordType CoordType;
        /// The shape type
        typedef typename MeshType::ShapeType ShapeType;
        /// Type for the vectors to hold coordinates etc.
        typedef LAFEM::DenseVectorBlocked<Mem::Main, CoordType, Index, MeshType::world_dim> CoordsBufferType;

      public:
        /// The mesh for the underlying transformation
        Geometry::RootMeshNode<MeshType>* _mesh_node;
        /// Coordinates, used for setting new boundary values etc.
        CoordsBufferType _coords_buffer;

      protected:
        /// Counter for number of function evaluations
        Index _num_func_evals;
        /// Counter for number of gradient evaluations
        Index _num_grad_evals;
        /// Counter for number of hessian evaluations
        Index _num_hess_evals;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] mesh_node_
         * The RootMeshNode this mesh optimiser will refer to.
         *
         */
        explicit MeshQualityFunctional(Geometry::RootMeshNode<MeshType>* mesh_node_) :
          _mesh_node(mesh_node_),
          _coords_buffer(mesh_node_->get_mesh()->get_num_entities(0), CoordType(0)),
          _num_func_evals(0),
          _num_grad_evals(0),
          _num_hess_evals(0)
          {
            XASSERTM(mesh_node_ != nullptr, "MeshNode must not be nullptr.");

            mesh_to_buffer();
          }

        explicit MeshQualityFunctional():
          _mesh_node(nullptr),
          _coords_buffer(),
          _num_func_evals(0),
          _num_grad_evals(0),
          _num_hess_evals(0)
          {
          }

        /// \brief Virtual destructor
        virtual ~MeshQualityFunctional()
        {
          // Just set it to nullptr since we did not allocate the memory for this
          _mesh_node = nullptr;
          _coords_buffer.clear();
        }

        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        virtual String name() const
        {
          return "MeshQualityFunctional<"+MeshType_::name()+">";
        }

        /// \returns The root mesh
        MeshType* get_mesh()
        {
          XASSERT(_mesh_node != nullptr);
          return _mesh_node->get_mesh();
        }

        /// \returns The root mesh as const pointer
        const MeshType* get_mesh() const
        {
          XASSERT(_mesh_node != nullptr);
          return _mesh_node->get_mesh();
        }

        /// \brief Gets the coordinates from the underlying mesh and saves them in _coords_buffer.
        virtual void mesh_to_buffer()
        {
          const typename MeshType::VertexSetType& vertex_set = get_mesh()->get_vertex_set();

          for(Index i(0); i < get_mesh()->get_num_entities(0); ++i)
            _coords_buffer(i, vertex_set[i]);
        }

        /// \brief Sets the coordinates in the underlying mesh to _coords_buffer.
        virtual void buffer_to_mesh()
        {
          typename MeshType::VertexSetType& vertex_set = get_mesh()->get_vertex_set();

          for(Index i(0); i < get_mesh()->get_num_entities(0); ++i)
            vertex_set[i] = _coords_buffer(i);
        }

        /**
         * \brief Gets the coords buffer
         *
         * \returns A reference to the coords buffer for manipulation.
         */
        CoordsBufferType& get_coords()
        {
          return _coords_buffer;
        }

        /**
         * \brief Gets the coords buffer
         *
         * \returns A const reference to the coords buffer for manipulation.
         */
        const CoordsBufferType& get_coords() const
        {
          return _coords_buffer;
        }

        /**
         * \returns The number of times the functional value was computed.
         */
        Index get_num_func_evals() const
        {
          return _num_func_evals;
        }

        /**
         * \returns The number of times the gradient was computed.
         */
        Index get_num_grad_evals() const
        {
          return _num_grad_evals;
        }

        /**
         * \returns The number of times the hessian was computed.
         */
        Index get_num_hess_evals() const
        {
          return _num_hess_evals;
        }

        /**
         * \brief Resets all evaluation counts
         */
        void reset_num_evals()
        {
          _num_func_evals = Index(0);
          _num_grad_evals = Index(0);
          _num_hess_evals = Index(0);
        }

        /**
         * \brief Performs one-time initialisations
         *
         * Because of the whole class inheritance hierarchy, some one time initialisations cannot be performed in the
         * constructors.
         */
        virtual void init() = 0;

        ///// \brief Prepares the mesh optimiser for application
        //virtual void prepare(VectorTypeR&) = 0;

        /**
         * \brief Computes a quality indicator concerning the cell sizes
         *
         * \returns The relative cell size quality indicator.
         */
        virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max,
        CoordType& vol_min, CoordType& vol_max) const = 0;

    }; // class MeshQualityFunctional

  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_MESH_QUALITY_FUNCTIONAL_HPP
