#pragma once
#ifndef KERNEL_GEOMETRY_DUDV_SMOOTHER_HPP
#define KERNEL_GEOMETRY_DUDV_SMOOTHER_HPP 1

#include <kernel/archs.hpp>

#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/common_operators.hpp>            // for DuDvOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/proto_solver.hpp>
#include <kernel/lafem/sparse_matrix_csr_blocked.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/geometry/mesh_smoother/mesh_smoother.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Mesh optimiser based on minimisation of harmonic energy
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \author Jordi Paul
     *
     */
    template<typename DataType_, typename MemType_, typename TrafoType_>
    class DuDvSmoother :
      public MeshSmoother<DataType_, MemType_, TrafoType_>
    {
      public:
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;

        /// Our base class
        typedef MeshSmoother<DataType, MemType, TrafoType> BaseClass;

        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;

        /// FE space the transformation lives in
        typedef Space::Lagrange1::Element<TrafoType> SpaceType;

        /// Block data types
        typedef Tiny::Vector<DataType, MeshType::shape_dim> TinyVectorType;
        /// Vector type for the coupled system for the coordinates
        typedef LAFEM::DenseVectorBlocked<MemType, DataType, Index, MeshType::shape_dim> VectorType;

        /// Type for the system matrix
        typedef LAFEM::SparseMatrixCSRBlocked<MemType, DataType, Index, MeshType::shape_dim, MeshType::shape_dim>
          MatrixType;
        /// Filter for out system
        typedef LAFEM::UnitFilterBlocked<MemType, DataType, Index, MeshType::shape_dim> FilterType;

      protected:
        /// Transformation
        SpaceType _trafo_space;
        /// System matrix
        MatrixType _sys_matrix;
        /// Right hand side for the mesh problem
        VectorType _vec_rhs;
        /// Cubature factory, for P1/Q1 transformations in 2d degree 5 is enough
        Cubature::DynamicFactory _cubature_factory;
        /// Filter for the mesh problem
        FilterType _filter;
        /// Assembler for the boundary values for the mesh problem
        Assembly::UnitFilterAssembler<MeshType> _dirichlet_asm;

      public:
        /// Constructor
        explicit DuDvSmoother(TrafoType& trafo_) :
          BaseClass(trafo_),
          _trafo_space(trafo_),
          _sys_matrix(),
          _vec_rhs(trafo_.get_mesh().get_num_entities(0),DataType(0)),
          _cubature_factory("auto-degree:5"),
          _filter(trafo_.get_mesh().get_num_entities(0)),
          _dirichlet_asm()
          {
            /// Type for the boundary mesh
            typedef typename Geometry::MeshPart<MeshType> BoundaryType;
            /// Factory for the boundary mesh
            typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

            // Get the boundary set
            BoundaryFactoryType boundary_factory(this->_mesh);
            BoundaryType boundary(boundary_factory);

            _dirichlet_asm.add_mesh_part(boundary);

          }

        virtual ~DuDvSmoother()
        {
        }

        /// \brief Initialises parts of the MeshSmoother not set in in the constructor
        virtual void init() override
        {
          BaseClass::init();

          Assembly::SymbolicMatrixAssembler<>::assemble1(_sys_matrix, _trafo_space);
        }


        /// \brief Optimises the mesh according to the criteria implemented in the mesh smoother.
        virtual void optimise() override
        {
          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));

          prepare();

          VectorType coords_blocked(nvertices, DataType(0));
          for (Index i(0); i < nvertices; ++i)
          {
            Tiny::Vector<DataType, MeshType::world_dim> tmp;

            for (int d(0); d < MeshType::world_dim; ++d)
               tmp(Index(d)) = this->_coords[d](i);

            coords_blocked(i, tmp);
          }

          _vec_rhs.format();

          _dirichlet_asm.assemble(_filter, _trafo_space, coords_blocked);

          _filter.filter_rhs(_vec_rhs);
          _filter.filter_sol(coords_blocked);

          // Some preconditioners/solvers do not care about filtered matrices (SSOR, PCG), but other do (UMFPACK) or
          // might (ILU)
          _filter.filter_mat(_sys_matrix);
          // Create a SSOR preconditioner
          // This is not implemented for SparseMatrixCSRBlocked yet, so instead we use NO preconditioner at all
          //LAFEM::PreconWrapper<MatrixType, LAFEM::DiagonalPreconditioner> precond();//_sys_matrix);
          // Create a PCG solver
          LAFEM::PCGSolver<MatrixType, FilterType> solver(_sys_matrix, _filter, nullptr);//&precond);
          // Enable convergence plot
          solver.set_plot(false);
          solver.set_max_iter(5000);
          // Initialise the solver
          solver.init();

          // Solve the system and correct coords_blocked
          solver.correct(coords_blocked, _vec_rhs);

          // Release the solver
          solver.done();

          // Copy back the new coordinates from the blocked vector
          for (Index i(0); i < nvertices; ++i)
          {
            for (int d(0); d < MeshType::world_dim; ++d)
              this->_coords[d](i, coords_blocked(i)(Index(d)));
          }

          // Copy back the coordinates to the underlying mesh
          this->set_coords();

          return;
        }

        /// \copydoc MeshSmoother::prepare()
        virtual void prepare() override
        {
          _sys_matrix.format();

          // For the CSRBlocked version to come
          Assembly::Common::DuDvOperatorBlocked<MeshType::shape_dim> my_operator;
          Assembly::BilinearOperatorAssembler::assemble_block_matrix1(
            _sys_matrix,           // the matrix that receives the assembled operator
            my_operator, // the operator that is to be assembled
            _trafo_space,            // the finite element space in use
            _cubature_factory  // the cubature factory to be used for integration
            );

          return;
        }

    }; // class DuDvSmoother

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_DUDV_SMOOTHER_HPP
