#pragma once
#ifndef KERNEL_MESHOPT_DUDV_SMOOTHER_HPP
#define KERNEL_MESHOPT_DUDV_SMOOTHER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/common_operators.hpp>            // for DuDvOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/sparse_matrix_csr_blocked.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/meshopt/mesh_smoother.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/space/lagrange1/element.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /// \cond Internal
    namespace Intern
    {
      template<typename Mem_>
      struct DuDvMSSolverParameters
      {
      };
    }
    /// \endcond
    /**
     * \brief Mesh optimiser based on minimisation of harmonic energy
     *
     * \tparam Mem_
     * Memory architecture for the solver (not the mesh)
     *
     * \tparam DT_
     * Data type for the solver (not the mesh)
     *
     * \tparam IT_
     * Index type for the solver (not the mesh)
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \author Jordi Paul
     *
     */
    template<typename Mem_, typename DT_, typename IT_, typename TrafoType_>
    class DuDvSmoother :
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

        /// Type of the system matrix for the assembly
        typedef LAFEM::SparseMatrixCSRBlocked<Mem::Main, CoordType, Index, MeshType::world_dim, MeshType::world_dim>
          MatrixType;
        /// Type of the system matrix for the solver
        typedef typename Intern::DuDvMSSolverParameters<Mem_>::template
          SolverMatrixType<DT_, IT_, MeshType::world_dim> SolverMatrixType;
        /// Filter for Dirichlet boundary conditions
        typedef LAFEM::UnitFilterBlocked<Mem_, DT_, IT_, MeshType::world_dim> DirichletFilterType;
        /// Filter for slip boundary conditions
        typedef LAFEM::SlipFilter<Mem_, DT_, IT_, MeshType::world_dim> SlipFilterType;
        /// Combined filter
        typedef LAFEM::FilterChain<SlipFilterType, DirichletFilterType> FilterType;

        /// Finite Element space for the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space TrafoSpace;
        /// Maximum polynomial degree
        // 2 * (degree of trafo) for both trial and test spaces, +1 for safety reasons
        // This could be decreased by the degree of the operator, i.e. 2 for Du:Dv
        static constexpr int _local_degree = 4*TrafoSpace::local_degree + 1;

      protected:
        /// Transformation
        TrafoType _trafo;
        /// FE space the transformation is from
        TrafoSpace _trafo_space;
        /// System matrix
        MatrixType _sys_matrix;
        /// Right hand side for the mesh problem
        typename SolverMatrixType::VectorTypeL _vec_rhs;
        /// Cubature factory, for P1/Q1 transformations in 2d degree 5 is enough
        Cubature::DynamicFactory _cubature_factory;
        /// Filter for the mesh problem
        FilterType _filter;
        /// Assembler for the Dirchlet boundary values for the mesh problem
        Assembly::UnitFilterAssembler<MeshType> _dirichlet_asm;
        /// Assembler for the slip boundary values for the mesh problem
        Assembly::SlipFilterAssembler<MeshType> _slip_asm;
        /// List of boundary identifiers for enforcing Dirichlet boundary conditions
        std::deque<String> _dirichlet_list;
        /// List of boundary identifiers for enforcing slip boundary conditions
        std::deque<String> _slip_list;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] rmn_
         * The RootMeshNode representing the tree of root mesh, all of its MeshParts and Charts
         *
         * \param[in] dirichlet_list_
         * List of boundary identifiers for enforcing Dirichlet boundary conditions, can be empty
         *
         * \param[in] slip_list_
         * List of boundary identifiers for enforcing slip boundary conditions, can be empty
         *
         */
        explicit DuDvSmoother(Geometry::RootMeshNode<MeshType>* rmn_,
        std::deque<String>& dirichlet_list_, std::deque<String>& slip_list_) :
          BaseClass(rmn_),
          _trafo(*(rmn_->get_mesh())),
          _trafo_space(_trafo),
          _sys_matrix(),
          _vec_rhs(rmn_->get_mesh()->get_num_entities(0),DataType(0)),
          _cubature_factory("auto-degree:"+stringify(int(_local_degree))),
          _filter(),
          _dirichlet_asm(),
          _slip_asm(*(rmn_->get_mesh())),
          _dirichlet_list(dirichlet_list_),
          _slip_list(slip_list_)
          {
            // Add all specified mesh parts to the Dirichlet filter assember
            for(auto& it : this->_dirichlet_list)
            {
              auto* mpp = this->_mesh_node->find_mesh_part(it);
              if(mpp != nullptr)
                _dirichlet_asm.add_mesh_part(*mpp);
            }

            // Add all specified mesh parts to the slip filter assember
            for(auto& it : this->_slip_list)
            {
              auto* mpp = this->_mesh_node->find_mesh_part(it);
              if(mpp != nullptr)
                _slip_asm.add_mesh_part(*mpp);
            }

          }

        /// \brief Virtual destructor
        virtual ~DuDvSmoother()
        {
        }

        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        static String name()
        {
          return "DuDvSmoother<"+MeshType::name()+">";
        }

        /**
         * \brief Performs one-time initialisations
         *
         * This is not done in the constructor for the case that the system matrix gets overwritten by a derived
         * class, so the unused system matrix of THIS class is not assembled symbolically
         */
        virtual void init() override
        {
          // Assemble the homogeneous slip filter
          _slip_asm.assemble(_filter.template at<0>(), _trafo_space);
          // Symbolically assemble the system matrix
          Assembly::SymbolicMatrixAssembler<>::assemble1(_sys_matrix, _trafo_space);
        }

        /// \brief Optimises the mesh according to the criteria implemented in the mesh smoother.
        virtual void optimise() override
        {
          prepare();

          // Declare solution vector and matrix for the solver
          typename SolverMatrixType::VectorTypeL sol;
          SolverMatrixType solver_matrix;

          // Convert matrix and solution for the solver. If the types are the same, this just shallow copying.
          sol.convert(this->_coords);
          solver_matrix.convert(_sys_matrix);

          // Assemble Dirichlet boundary conditions from sol, as this->_coords contains the new boundary coordinates
          _dirichlet_asm.assemble(_filter.template at<1>(), _trafo_space, sol);

          _vec_rhs.format();

          _filter.filter_rhs(_vec_rhs);

          // Create the solver
          auto solver = Intern::DuDvMSSolverParameters<Mem_>::solver(solver_matrix, _filter);
          // Initialise the solver
          solver->init();
          // Uncomment this to see the convergence behaviour
          //solver->set_plot(true);

          // Set a tighter relative convergence criterion than the default
          solver->set_tol_rel(Math::pow<DataType>(Math::eps<DataType>(), DataType(0.75)));

          // Solve the system and correct coords_blocked
          solver->correct(sol, _vec_rhs);

          // Release the solver
          solver->done();

          // Convert the solution
          this->_coords.convert(sol);

          // Copy back the coordinates to the underlying mesh
          this->set_coords();

          return;
        }

        /// \copydoc MeshSmoother::prepare()
        virtual void prepare() override
        {
          // Adapt all slip boundaries
          for(auto& it : _slip_list)
            this->_mesh_node->adapt_by_name(it);

          // Assemble homogeneous slip boundary conditions, as the outer normal could have changed
          _slip_asm.assemble(_filter.template at<0>(), _trafo_space);

          _sys_matrix.format();

          Assembly::Common::DuDvOperatorBlocked<MeshType::world_dim> my_operator;

          Assembly::BilinearOperatorAssembler::assemble_block_matrix1(
            _sys_matrix,           // the matrix that receives the assembled operator
            my_operator, // the operator that is to be assembled
            _trafo_space,            // the finite element space in use
            _cubature_factory  // the cubature factory to be used for integration
            );

          return;
        }

    }; // class DuDvSmoother

    /// \cond Internal
    namespace Intern
    {
      template<>
      struct DuDvMSSolverParameters<Mem::Main>
      {
        template<typename DT_, typename IT_, int BlockSize_>
        using SolverMatrixType = LAFEM::SparseMatrixCSRBlocked<Mem::Main, DT_, IT_, BlockSize_, BlockSize_>;

        template<typename Matrix_, typename Filter_>
        static std::shared_ptr<Solver::IterativeSolver<typename Matrix_::VectorTypeR>> solver(Matrix_& matrix, Filter_& filter)
        {
          auto my_solver = Solver::new_pcg<Matrix_, Filter_>(matrix, filter);
          return my_solver;
        }
      };
#ifdef FEAST_BACKENDS_CUDA
      /**
       * Dummy specialisation with nonsensical stuff, do not use. Here
       *
       * Since not all filter operations are implemented for UnitFilterBlocked on CUDA, this will not work (yet).
       * There is not blocked matrix type for CUDA yet, and SparseMatrixCSR is a nonsensical format for CUDA.
       * There is no preconditioner yet.
       *
       */
      template<>
      struct DuDvMSSolverParameters<Mem::CUDA>
      {
        template<typename DT_, typename IT_, int BlockSize_>
        using SolverMatrixType = LAFEM::SparseMatrixCSRBlocked<Mem::CUDA, DT_, IT_, BlockSize_, BlockSize_>;

        template<typename Matrix_, typename Filter_>
        static std::shared_ptr<Solver::IterativeSolvertypename Matrix_::VectorTypeR>> solver(Matrix_& matrix, Filter_& filter)
        {
          auto precond = Solver::new_jacobi_precond(matrix, filter);
          auto my_solver = Solver::new_pcg<Matrix_, Filter_>(matrix, filter, precond);
          return my_solver;
        }
      };
#endif
    } // namespace Intern
    /// \endcond

  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_DUDV_SMOOTHER_HPP
