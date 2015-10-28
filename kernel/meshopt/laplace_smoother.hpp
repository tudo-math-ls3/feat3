#pragma once
#ifndef KERNEL_MESHOPT_LAPLACE_SMOOTHER_HPP
#define KERNEL_MESHOPT_LAPLACE_SMOOTHER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/meshopt/mesh_smoother.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/ssor_precond.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /// \cond Internal
    namespace Intern
    {
      template<typename Mem_>
      struct LMSSolverParameters
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
    class LaplaceSmoother :
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
        typedef LAFEM::SparseMatrixCSR<Mem::Main, CoordType, Index> MatrixType;
        /// Type of the system matrix for the solver
        typedef typename Intern::LMSSolverParameters<Mem_>::template SolverMatrixType<DT_, IT_> SolverMatrixType;
        /// Filter for our system
        typedef LAFEM::UnitFilter<Mem_, DT_, IT_> FilterType;

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
        /// Assembler for the boundary values for the mesh problem
        Assembly::UnitFilterAssembler<MeshType> _dirichlet_asm;
        /// List of boundary identifiers for enforcing Dirichlet boundary conditions
        std::deque<String>& _dirichlet_list;

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
         * \note As the components are decoupled, slip BCs are not supported and slip_list_ is here just for
         * interface compatibility
         *
         */
        explicit LaplaceSmoother(Geometry::RootMeshNode<MeshType>* rmn_,
        std::deque<String>& dirichlet_list_, std::deque<String>& DOXY(slip_list_)) :
          BaseClass(rmn_),
          _trafo(*(rmn_->get_mesh())),
          _trafo_space(_trafo),
          _sys_matrix(),
          _vec_rhs(rmn_->get_mesh()->get_num_entities(0),DataType(0)),
          _cubature_factory("auto-degree:"+stringify(int(_local_degree))),
          _filter(),
          _dirichlet_asm(),
          _dirichlet_list(dirichlet_list_)
          {
            // Add all boundaries specified to the dirichlet assembler
            for(auto& it : this->_dirichlet_list)
            {
              auto* mpp = this->_mesh_node->find_mesh_part(it);
              if(mpp != nullptr)
                _dirichlet_asm.add_mesh_part(*mpp);
            }

          }

        /// \brief Virtual destructor
        virtual ~LaplaceSmoother()
        {
        }

        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        static String name()
        {
          return "LaplaceSmoother<"+MeshType::name()+">";
        }

        /**
         * \brief Performs one-time initialisations
         *
         * This is not done in the constructor for the case that the system matrix gets overwritten by a derived
         * class, so the unused system matrix of THIS class is not assembled symbolically
         */
        virtual void init() override
        {
          // Symbolically assemble the system matrix
          Assembly::SymbolicMatrixAssembler<>::assemble1(_sys_matrix, _trafo_space);
        }

        /// \brief Optimises the mesh according to the criteria implemented in the mesh smoother.
        virtual void optimise() override
        {
          prepare();

          // Convert matrix for the solver. If the types are the same, this just shallow copying.
          SolverMatrixType solver_matrix;
          solver_matrix.convert(_sys_matrix);

          // Create solution vector for the solver
          typename SolverMatrixType::VectorTypeR sol(this->get_mesh()->get_num_entities(0));
          // Some preconditioners/solvers do not care about filtered matrices (SSOR, PCG), but other do (UMFPACK) or
          // might (ILU)
          // Assemble bogus homogeneous filter so we can filter the matrix
          // _dirichlet_asm.assemble(_filter, _trafo_space);
          //_filter.filter_mat(solver_matrix);

          // Create the solver
          auto solver = Intern::LMSSolverParameters<Mem_>::solver(solver_matrix, _filter);
          // Initialise the solver
          solver->init();
          // Set a tighter relative convergence criterion than the default
          solver->set_tol_rel(Math::pow<DataType>(Math::eps<DataType>(), DataType(0.75)));

          for(int d(0); d < MeshType::world_dim; ++d)
          {
            // Copy the d-th component of _coords to sol
            for(Index i(0); i < this->get_mesh()->get_num_entities(0); ++i)
              sol(i, this->_coords(i)(d));

            // Assemble Dirichlet boundary conditions from sol, as this->_coords contains the new boundary
            // coordinates
            _dirichlet_asm.assemble(_filter, _trafo_space, sol);

            _vec_rhs.format();
            _filter.filter_rhs(_vec_rhs);

            // Correct our initial solution vector
            auto st = solver->correct(sol, _vec_rhs);

            // Print solver summary
            std::cout << "Component " << d << ": " << solver->get_plot_name() << ": " << st << ", "
            << solver->get_num_iter() << " its, defect initial/final: " << stringify_fp_sci(solver->get_def_initial())
            << " / " << stringify_fp_sci(solver->get_def_final()) << std::endl;

            // Convert the solution back (so we have one big download from Mem::CUDA)
            typename MatrixType::VectorTypeR sol_backcopy;
            sol_backcopy.convert(sol);

            // Copy solution to _coords
            for(Index i(0); i < this->get_mesh()->get_num_entities(0); ++i)
            {
              Tiny::Vector<CoordType, MeshType::world_dim> tmp(this->_coords(i));
              tmp(d) = sol_backcopy(i);
              this->_coords(i, tmp);
            }
          }

          // Print solver summary
          std::cout << "Component " << d << ": " << solver->get_plot_name() << ": " << st << ", "
          << solver->get_num_iter() << " its, defect initial/final: " << stringify_fp_sci(solver->get_def_initial())
          << " / " << stringify_fp_sci(solver->get_def_final()) << std::endl;

          // Release the solver
          solver->done();

          // Copy back the coordinates to the underlying mesh
          this->set_coords();

          return;
        }

        /// \copydoc MeshSmoother::prepare()
        virtual void prepare() override
        {
          _sys_matrix.format();
          Assembly::Common::LaplaceOperator laplace_operator;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(
            _sys_matrix,           // the matrix that receives the assembled operator
            laplace_operator, // the operator that is to be assembled
            _trafo_space,            // the finite element space in use
            _cubature_factory  // the cubature factory to be used for integration
            );

          return;
        }

    }; // class LaplaceSmoother

    /// \cond Internal
    namespace Intern
    {
      template<>
      struct LMSSolverParameters<Mem::Main>
      {
        template<typename DT_, typename IT_>
        using SolverMatrixType = LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>;

        template<typename Matrix_, typename Filter_>
        static std::shared_ptr<Solver::IterativeSolver<typename Matrix_::VectorTypeR>> solver(Matrix_& matrix, Filter_& filter)
        {
          auto precond = Solver::new_ssor_precond(matrix, filter);
          auto my_solver = Solver::new_pcg<Matrix_, Filter_>(matrix, filter, precond);
          return my_solver;
        }
      };
#ifdef FEAST_BACKENDS_CUDA
      template<>
      struct LMSSolverParameters<Mem::CUDA>
      {
        template<typename DT_, typename IT_>
        using SolverMatrixType = LAFEM::SparseMatrixELL<Mem::CUDA, DT_, IT_>;

        template<typename Matrix_, typename Filter_>
        static std::shared_ptr<Solver::IterativeSolver<typename Matrix_::VectorTypeR>> solver(Matrix_& matrix, Filter_& filter)
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
#endif // KERNEL_MESHOPT_LAPLACE_SMOOTHER_HPP
