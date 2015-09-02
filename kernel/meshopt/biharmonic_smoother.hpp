#pragma once
#ifndef KERNEL_MESHOPT_BIHARMONIC_SMOOTHER_HPP
#define KERNEL_MESHOPT_BIHARMONIC_SMOOTHER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/meshopt/mesh_smoother.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /// \cond Internal
    namespace Intern
    {
      template<typename Mem_>
      struct BHMSSolverParameters
      {
      };
    }
    /// \endcond

    /**
     * \brief Mesh optimiser based on minimisation of biharmonic energy
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
    template<typename Mem_, typename DT_, typename IT_, typename TrafoType_>
    class BiharmonicSmoother :
      public MeshSmoother<typename TrafoType_::MeshType>
    {
      public:
        /// Memory architecture
        typedef Mem_ MemType;
        /// Our datatype
        typedef DT_ DataType;
        /// Our index type
        typedef IT_ IndexType;
        /// Type for the transformation
        typedef TrafoType_ TrafoType;

        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Our base class
        typedef MeshSmoother<MeshType> BaseClass;

        /// Type of the system matrix blocks for the assembly
        typedef LAFEM::SparseMatrixCSR<Mem::Main, CoordType, Index> SubMatrixType;
        /// TYpe of the system matrix for the assembly
        typedef LAFEM::SaddlePointMatrix<SubMatrixType> MatrixType;

        /// Type of the system matrix for the solver
        typedef typename Intern::BHMSSolverParameters<Mem_>::template SolverMatrixType<DT_, IT_> SolverMatrixType;
        /// Type for the vectors to hold coordinates etc.
        typedef LAFEM::DenseVector<Mem_, DT_, IT_> SolverSubVectorType;
        // Chose your poison: UnitFilter means bogus UnitFilter values for the Laplacian, if they are zero then the
        // whole method is equivalent to just using the Laplace smoother.
        // Or use a filter for Neuman BVs for the Laplacian. Unfortunately, only homogeneous Neumann BVs are
        // implemented now, and this does not make any sense for the mesh coordinate distribution point of view.
        /// Filter type for the first component of the block system, meaning the normal derivative of the coordinates
        typedef LAFEM::NoneFilter<Mem_, DT_, IT_> SubFilterType0;
        /// Filter type for the second component of the block system, meaning the coordinates
        typedef LAFEM::UnitFilter<Mem_, DT_, IT_> SubFilterType1;
        /// Filter type for the block system
        typedef LAFEM::TupleFilter<SubFilterType0, SubFilterType1> FilterType;

        /// Finite Element space for the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space TrafoSpace;
        /// Maximum polynomial degree
        // 2 * (degree of trafo) for both trial and test spaces, +1 for safety reasons
        // This could be decreased by the degree of the operator, i.e. 2 for Du:Dv
        static constexpr int _local_degree = 4*TrafoSpace::local_degree + 1;

      protected:
        /// Transformation
        TrafoSpace _trafo_space;
        /// System matrix
        MatrixType _sys_matrix;
        /// Right hand side for the mesh problem
        typename SolverMatrixType::VectorTypeL _vec_rhs;
        /// Cubature factory, for P1/Q1 transformations in 2d degree 5 is enough
        Cubature::DynamicFactory _cubature_factory;
        /// Filter for the mesh problem
        FilterType _filter;
        /// Assembler for the boundary values for the mesh problem
        Assembly::UnitFilterAssembler<MeshType> _dirichlet_asm;

      public:
        /// Constructor
        explicit BiharmonicSmoother(TrafoType& trafo_) :
          BaseClass(trafo_.get_mesh()),
          _trafo_space(trafo_),
          _sys_matrix(),
          _vec_rhs(),
          _cubature_factory("auto-degree:"+stringify(int(_local_degree))),
          _filter(),
          _dirichlet_asm()
          {
            /// Type for the boundary mesh
            typedef typename Geometry::MeshPart<MeshType> BoundaryType;
            /// Factory for the boundary mesh
            typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

            // Total number of vertices in the mesh
            Index nvertices(this->_mesh.get_num_entities(0));

            // Get the boundary set
            BoundaryFactoryType boundary_factory(this->_mesh);
            BoundaryType boundary(boundary_factory);

            // This does not have to be set if we use a NoneFilter
            //_filter.template at<0>() = std::move(SubFilterType0(nvertices));
            _filter.template at<1>() = std::move(SubFilterType1(nvertices));

            _vec_rhs.template at<0>() = std::move(SolverSubVectorType(nvertices));
            _vec_rhs.template at<1>() = std::move(SolverSubVectorType(nvertices));

            _dirichlet_asm.add_mesh_part(boundary);

          }

        virtual ~BiharmonicSmoother()
        {
        }

        /// \brief Initialises parts of the MeshSmoother not set in in the constructor
        virtual void init() override
        {
          BaseClass::init();
          // Assemble the A block of the SaddlePointMatrix, in this case it will be the mass matrix
          Assembly::SymbolicMatrixAssembler<>::assemble1(_sys_matrix.block_a(), _trafo_space);
          // The other blocks have the same sparsity pattern, just different values
          _sys_matrix.block_b() = _sys_matrix.block_a().clone(LAFEM::CloneMode::Weak);
          _sys_matrix.block_d() = _sys_matrix.block_a().clone(LAFEM::CloneMode::Weak);
        }

        /// \brief Optimises the mesh according to the criteria implemented in the mesh smoother.
        virtual void optimise() override
        {
          // Calling prepare() assembles the system matrix numerically
          prepare();

          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));

          // Convert matrix for the solver. If the types are the same, this just shallow copying.
          SolverMatrixType solver_matrix;
          solver_matrix.convert(_sys_matrix);

          // Create solution vector for the solver
          typename SolverMatrixType::VectorTypeR sol;
          sol.template at<0>() = std::move(SolverSubVectorType(nvertices));
          sol.template at<1>() = std::move(SolverSubVectorType(nvertices));

          // Some preconditioners/solvers do not care about filtered matrices (SSOR, PCG), but other do (UMFPACK) or
          // might (ILU)
          //_filter.filter_mat(solver_matrix);

          // Create the solver
          auto solver = Intern::BHMSSolverParameters<Mem_>::solver(solver_matrix, _filter);
          // Initialise the solver
          solver->init();
          solver->set_plot(true);

          for(int d(0); d < MeshType::world_dim; ++d)
          {
            sol.template at<0>().format();
            // Copy the d-th component of _coords to sol
            for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
              sol.template at<1>()(i, this->_coords(i)(d));

            // Assemble Dirichlet boundary conditions from sol, as this->_coords contains the new boundary
            // coordinates
            _dirichlet_asm.assemble(_filter.template at<1>(), _trafo_space, sol.template at<1>());

            _vec_rhs.format();
            _filter.filter_rhs(_vec_rhs);

            // Correct our initial solution vector
            solver->correct(sol, _vec_rhs);

            // Convert the solution back (so we have one big download from Mem::CUDA)
            typename SubMatrixType::VectorTypeR sol_backcopy;
            sol_backcopy.convert(sol.template at<1>());

            // Copy solution to _coords
            for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
            {
              Tiny::Vector<CoordType, MeshType::world_dim> tmp(this->_coords(i));
              tmp(d) = sol_backcopy(i);
              this->_coords(i, tmp);
            }
          }

          // Release the solver
          solver->done();

          // Copy back the coordinates to the underlying mesh
          this->set_coords();

          return;
        }

        /// \copydoc MeshSmoother::prepare()
        virtual void prepare() override
        {
          //_sys_matrix.format();
          //Assembly::Common::LaplaceOperator laplace_operator;
          //Assembly::BilinearOperatorAssembler::assemble_matrix1(
          //  _sys_matrix,           // the matrix that receives the assembled operator
          //  laplace_operator, // the operator that is to be assembled
          //  _trafo_space,            // the finite element space in use
          //  _cubature_factory  // the cubature factory to be used for integration
          //  );

          _sys_matrix.block_a().format();
          Assembly::Common::IdentityOperator identity_operator;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(
            _sys_matrix.block_a(),           // the matrix that receives the assembled operator
            identity_operator, // the operator that is to be assembled
            _trafo_space,            // the finite element space in use
            _cubature_factory,  // the cubature factory to be used for integration
            DataType(-1)
            );

          _sys_matrix.block_b().format();
          Assembly::Common::LaplaceOperator laplace_operator;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(
            _sys_matrix.block_b(),           // the matrix that receives the assembled operator
            laplace_operator, // the operator that is to be assembled
            _trafo_space,            // the finite element space in use
            _cubature_factory  // the cubature factory to be used for integration
            );

          _sys_matrix.block_d() = _sys_matrix.block_b().shared();

          return;
        }
    }; // class BiharmonicSmoother

    /// \cond Internal
    namespace Intern
    {
      template<>
      struct BHMSSolverParameters<Mem::Main>
      {
        template<typename DT_, typename IT_>
        using SolverMatrixType = LAFEM::SaddlePointMatrix<LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>>;

        template<typename Matrix_, typename Filter_>
        static std::shared_ptr<Solver::PreconditionedIterativeSolver<Matrix_, Filter_>> solver(Matrix_& matrix, Filter_& filter)
        {
          // auto precond = Solver::new_ssor_precond(matrix, filter);
          auto my_solver = Solver::new_pcg<Matrix_, Filter_>(matrix, filter, nullptr);
          return my_solver;
        }
      };
#ifdef FEAST_BACKENDS_CUDA
      template<>
      struct BHMSSolverParameters<Mem::CUDA>
      {
        template<typename DT_, typename IT_>
        using SolverMatrixType = LAFEM::SparseMatrixELL<Mem::CUDA, DT_, IT_>;

        template<typename Matrix_, typename Filter_>
        static std::shared_ptr<Solver::IterativeSolver<Matrix_, Filter_>> solver(Matrix_& matrix, Filter_& filter)
        {
          auto precond = Solver::new_jacobi_precond(matrix, filter);
          auto my_solver = Solver::new_pcg<Matrix_, Filter_>(matrix, filter, precond);
          return my_solver;
        }
      };
#endif
    }
    /// \endcond


  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_BIHARMONIC_SMOOTHER_HPP
