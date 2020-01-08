// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_MESHOPT_MESH_QUALITY_FUNCTIONAL_HPP
#define KERNEL_MESHOPT_MESH_QUALITY_FUNCTIONAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/meshopt/base.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>

namespace FEAT
{
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

        /**
         * \brief Empty standard constructor
         *
         * This is needed because derived classes like DuDvFunctional are used in Global::Matrix, which needs an
         * empty standard constructor.
         *
         */
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
        static String name()
        {
          return "MeshQualityFunctional<"+MeshType_::name()+">";
        }

        /**
         * \brief Adds relevant quantities of this object to a VTK exporter
         *
         * \param[in, out] exporter
         * The exporter to add our data to.
         *
         */
        virtual void add_to_vtk_exporter(Geometry::ExportVTK<MeshType>& DOXY(exporter)) const
        {
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

        /// \returns The root mesh node
        Geometry::RootMeshNode<MeshType>* get_mesh_node() const
        {
          return _mesh_node;
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
         * \param[out] lambda_min
         * Minimum of the optimal cell size lambda over all cells
         *
         * \param[out] lambda_max
         * Maximum of the optimal cell size lambda over all cells
         *
         * \param[out] vol_min
         * Minimum cell volume
         *
         * \param[out] vol_max
         * Maximum cell volume
         *
         * \param[out] vol
         * Total volume of the domain given by the mesh
         *
         * In a truly optimal mesh (consisting ONLY of reference cells of the right size), every cell's volume is
         * exactly lambda(cell). This is especially the goal for r-adaptivity.
         * So in an optimal mesh,
         * \f[
         *   \forall K \in \mathcal{T}_h: |K|/|\Omega| = \lambda(K)
         * \f]
         * so we compute the 1-norm of the vector
         * \f$(v)_i = \left| \frac{|K_i|}{\sum_j |K_j|} - \lambda(K_i) \right| \f$.
         *
         * \returns The relative cell size quality indicator.
         *
         * \note lambda_min, lambda_max, vol_min, and vol_max are all volume fractions.
         *
         */
        virtual CoordType compute_cell_size_defect(CoordType& DOXY(lambda_min), CoordType& DOXY(lambda_max),
        CoordType& DOXY(vol_min), CoordType& DOXY(vol_max), CoordType& DOXY(vol)) const = 0;

    }; // class MeshQualityFunctional

  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_MESH_QUALITY_FUNCTIONAL_HPP
