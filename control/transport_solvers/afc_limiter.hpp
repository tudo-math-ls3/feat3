// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>

#include <memory>

namespace FEAT
{
  namespace Control
  {
    namespace Intern
    {
      template<typename MatrixType_>
      struct MatrixTypeHelper
      {
        static constexpr bool is_global = false;
        typedef MatrixType_ LocalMatrixType;
        typedef LAFEM::UnitFilter<typename MatrixType_::DataType> FilterType;
        typedef LAFEM::UnitFilter<typename MatrixType_::DataType> LocalFilterType;
      };

      template<typename LocalMatrixType_, typename RowMirror_, typename ColMirror_>
      struct MatrixTypeHelper<Global::Matrix<LocalMatrixType_, RowMirror_, ColMirror_>>
      {
        static constexpr bool is_global = true;
        typedef LocalMatrixType_ LocalMatrixType;
        typedef Global::Filter<LAFEM::UnitFilter<typename LocalMatrixType_::DataType>, RowMirror_> FilterType;
        typedef LAFEM::UnitFilter<typename LocalMatrixType_::DataType> LocalFilterType;
        typedef RowMirror_ RowMirror;
        typedef ColMirror_ ColMirror;
      };
    }

    /**
     * \brief An algebraic flux correction solver for the instationary scalar advection, diffusion, reaction equation
     *
     * This interface solves the equation given by
     * \f$ M u + K u + L u + N u = f \\ u = g on \delta_+\Omega  \f$
     * where
     * \f$ M \f$ is the (lumped) mass matrix
     * \f$ K \f$ is a convetion matrix
     * \f$ L \f$ is a diffusion matrix
     * \f$ N \f$ is a scaled mass matrix
     *
     * \tparam MatrixType_ The matrix type to be used
     *
     * \todo Test on distorted mesh if parallel gets same result if serial
     * \todo add implicit MCL
     */
    template<typename MatrixType_>
    class AFCTransportSolver
    {
    public:
      typedef MatrixType_ MatrixType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef typename MatrixType::IndexType IndexType;
      static constexpr bool is_global = Intern::MatrixTypeHelper<MatrixType_>::is_global;
      typedef typename Intern::MatrixTypeHelper<MatrixType_>::LocalMatrixType LocalMatrixType;
      typedef typename LocalMatrixType::VectorTypeR LocalVectorType;

      typedef typename Intern::MatrixTypeHelper<MatrixType_>::FilterType FilterType;
      typedef typename Intern::MatrixTypeHelper<MatrixType_>::LocalFilterType LocalFilterType;

    protected:
      // our internal (to be filtered) matrices
      MatrixType  _matrix_mass, _matrix_conv, _matrix_stiff, _matrix_react;
      // Helper matrices
      MatrixType _matrix_art_diff;
      MatrixType _matrix_hlp;
      MatrixType _matrix_hlp_t; // MCL
      // L = M_L + theta * dt * (K - D + L + N + B_(L)) (Matrices evaluated at n + 1)
      MatrixType _matrix_left; // FCT impl, MCL impl
      // R = M_L - (1 - theta) * dt * (K - D + L + N + B_(L)) (Matrices evaluated at n)
      MatrixType _matrix_right;
      LocalMatrixType _matrix_schwarz;

      //our solver for the implicit diffusion part
      std::shared_ptr<Solver::IterativeSolver<VectorType>> _solver;

      // helper vectors
      VectorType _vec_mass_lump;
      VectorType _vec_mass_lump_inverse;
      // weak boundary conditions
      VectorType _vec_bc;

      // strong boundary conditions
      FilterType _filter;
      LocalFilterType _loc_filter;

      // known size? -> std::array<VectorType, k>, anyway reserve enough space beforehand
      std::vector<VectorType> _vec_hlp;

      // which type of theta scheme do we want
      DataType _theta = 0;
      // stepsize for the time stepping
      DataType _stepsize = DataType(-1);
      // which limiter do we want? 1: FCT, 2: MCL
      int _limiter = 1;
      // how do we want the boundary conditions to be applied? 1: weak, 2: strong, 3: both
      int _bc_config = 0;
      // config for solver:
      int _lin_solver = 1;
      int _precond = 1;

      // we need those for synchronising the limiter
      const Global::Gate<LocalVectorType, typename Intern::MatrixTypeHelper<MatrixType_>::RowMirror>* _row_gate;
      const Global::Gate<LocalVectorType, typename Intern::MatrixTypeHelper<MatrixType_>::ColMirror>* _col_gate;
      std::shared_ptr<Global::SynchMatrix<LocalMatrixType, typename Intern::MatrixTypeHelper<MatrixType_>::RowMirror>> _synch_matrix_ticket;
      LocalMatrixType _matrix_freq;

      bool _matrix_mass_set = false;
      bool _matrix_conv_set = false;
      bool _matrix_stiff_set = false;
      bool _matrix_react_set = false;
      bool _init_symbolic_called = false;
      bool _init_stepsize_called = false;
      bool _init_time_stepping_called = false;
      // should the diffusion (i.e. stiffnes matrix) also be limited?
      bool _limit_matrix_stiff = false;

      // stopwatches for time measuring
      StopWatch watch_solve, watch_init;

      // we need those to iterate over our matrices all the time
      IndexType _num_rows;

      void _set_loc_filter()
      {
        if constexpr(is_global)
        {
          _loc_filter = _filter.local().clone(LAFEM::CloneMode::Shallow);
        }
        else
        {
          _loc_filter = _filter.clone(LAFEM::CloneMode::Shallow);
        }
      }

      void _set_glob_filter()
      {
        if constexpr(is_global)
        {
          _filter.local() = _loc_filter.clone(LAFEM::CloneMode::Shallow);
        }
        else
        {
          _filter = _loc_filter.clone(LAFEM::CloneMode::Shallow);
        }
      }

    public:
      // Rule of 5
      AFCTransportSolver() = default;
      AFCTransportSolver(const AFCTransportSolver&) = delete;
      AFCTransportSolver& operator=(const AFCTransportSolver&) = delete;
      AFCTransportSolver(AFCTransportSolver&&) = default;
      AFCTransportSolver& operator=(AFCTransportSolver&&) = default;
      ~AFCTransportSolver() = default;

      LocalFilterType& get_filter()
      {
        return _loc_filter;
      }

      const LocalFilterType& get_filter() const
      {
        return _loc_filter;
      }

      void set_filter(const LocalFilterType& filter)
      {
        _loc_filter = filter.clone(LAFEM::CloneMode::Shallow);
        _set_glob_filter();
        _bc_config = 2;
      }

      /**
       * \brief Sets vector for weak boundary conditions
       *
       * That vector needs to be assembled according to
       * \int_{\partial\Omega} -min(v*n, 0) * f(x) * \phi dx
       *
       */
      void set_weak_bc(const VectorType& vec)
      {
        XASSERTM(_init_symbolic_called, "Need to call init_symbolic() first");
        _vec_bc.copy(vec);
        _bc_config = 1;
      }

      /**
       * \brief Returns the mass matrix \int_\Omega u \psi dx
       */
      MatrixType& get_mass_matrix()
      {
        _matrix_mass_set = true;
        return _matrix_mass;
      }

      /**
       * \brief Returns the convection matrix \int_\Omega (v\cdot \nabla) u \psi dx
       */
      MatrixType& get_conv_matrix()
      {
        _matrix_conv_set = true;
        return _matrix_conv;
      }

      /**
       * \brief Returns the sitffness matrix \int_\Omega \nabla u \cdot \nabla \psi dx
       */
      MatrixType& get_stiff_matrix()
      {
        _matrix_stiff_set = true;
        return _matrix_stiff;
      }

      /**
       * \brief Returns the reaction matrix \int_\Omega u\psi dx
       */
      MatrixType& get_react_matrix()
      {
        _matrix_react_set = true;
        return _matrix_react;
      }

      /**
       * \brief Wether or not to also limit the stiffness matrix
       *
       * Depending on mesh quality and diffusion tensor it may happen that the stiffnes matrix
       * has positive off diagonal elements, which can lead to the low order scheme not being
       * bound preserving. To counteract this, you can choose to incorporate the sitffness matrix
       * into the inital calculation of the artifical diffusion. This can lead to an order
       * reduction, but the scheme shouldn't violate the initial boundaries
       *
       * Should be set before calculating the artificial diffusion, i.e. before calling
       * init_stepsize()
       *
       * TODO: what happens to mcl when using this?
       */
      void limit_stiff_matrix(bool limit)
      {
        _limit_matrix_stiff = limit;
      }

      /**
       * \brief changes physical diffusion matrix to make it an M-Matrix
       *
       * remove positive off diagonal elements by substracting antidiuffsive part L^+, recieving
       * \tilde{L} = L - L^+, which should be an M-Matrix. Sort of a last resort, if limiting
       * the physical diffusion doesn't work properly. Will lead to more diffusive, possibly
       * inconsistent results
       *
       * TODO: 0 to 1 to 0 synchronization?
       */
      void stiff_mat_to_m()
      {
        XASSERTM(_init_symbolic_called, "Need to call init_symbolic() first");
        XASSERTM(_matrix_stiff_set, "stiffnes matrix not set");
        auto mat_l_plus = _matrix_stiff.local().clone(LAFEM::CloneMode::Weak);
        mat_l_plus.format();
        auto mat_stiff_t = _matrix_stiff.local().transpose();

        const Memory::TypedView<IndexType> row_ptr_view = _matrix_mass.local().row_ptr_view_r();
        const Memory::TypedView<IndexType> col_idx_view = _matrix_mass.local().col_idx_view_r();
        const Memory::TypedView<DataType> val_stiff = _matrix_stiff.local().val_view_r();
        const Memory::TypedView<DataType> val_stiff_t = mat_stiff_t.val_view_r();
        Memory::TypedView<DataType> val_plus = mat_l_plus.val_view_rw();
        for(IndexType row = 0; row < _num_rows; ++row)
        {
          const IndexType start = row_ptr_view[row];
          const IndexType end = row_ptr_view[row + 1];

          IndexType diag = end;
          DataType lump = DataType(0);

          for(IndexType k = start; k < end; ++k)
          {
            if(col_idx_view[k] == row)
              diag = k;
            else
            {
              auto val = Math::max(DataType(0), Math::max(val_stiff[k], val_stiff_t[k]));
              val_plus[k] += val;
              lump += val;
            }
          }
          val_plus[diag] += -lump;
        }
        // 0 to 1 to 0???!!!
        _matrix_stiff.local().axpy(mat_l_plus, -1);
      }

      /**
       * \brief Sets the time stepping parameter
       *
       * theta = 0: explicit Euler (or Heun's method when using that apply call)
       * theta = 1/2: Crank-Nicolsin
       * theta = 1: implicit Euler
       */
      void set_theta(DataType theta)
      {
        XASSERTM(!_init_symbolic_called, "theta can only be set before init_symbolic or after done_symbolic has been called");
        XASSERTM((theta >= 0) && (theta <= 1), "theta must be in the intervall [0, 1]");
        _theta = theta;
      }

      /**
       * \brief Choose the limiter
       *
       * 1: FCT (explicit and implicit)
       * 2: MCL (only explicit)
       */
      void set_limiter(int limiter)
      {
        XASSERTM(!_init_symbolic_called, "Limiter can only be set before init_symbolic or after done_symbolic has been called");
        XASSERTM((limiter == 1) || (limiter == 2), "Limiter must be 1 or 2");
        _limiter = limiter;
      }

      /**
       * \brief Choose a linear solver for implict time stepping
       *
       * 1: BiCGStab
       * 2: FGMRES
       */
      void set_lin_solver(IndexType solver)
      {
        XASSERTM(!_init_symbolic_called, "Linear solver can only be set before init_symbolic or after done_symbolic has been called");
        XASSERTM((solver == 1) || (solver == 2), "solver must be 1 or2");
        _lin_solver = solver;
      }

      /**
       * \brief Choose a linear solver for implict time stepping
       *
       * 1: Jacobi
       * 2: SSOR
       * 3: ILU
       */
      void set_precond(IndexType precond)
      {
        XASSERTM(!_init_symbolic_called, "Preconditioner can only be set before init_symbolic or after done_symbolic has been called");
        XASSERTM((precond == 1) || (precond == 2) || (precond == 3), "Preconditioner must be 1, 2 or 3");
        _precond = precond;
      }

      /**
       * \brief Initializes solver
       *
       * Calls init_symbolic() and then init_numeric().
       *
       * \note The main matrices _matrix_stiff, _matrix_mass, _matrix_diff and _matrix_react have to be
       *       allocated to their correct sizes beforehand
       */
      void init()
      {
        this->init_symbolic();
        this->init_numeric();
      }

      /**
       * \brief Initializes symbolic structure of solver
       *
       * This method performs all necessary allocations for the internal solving steps.
       *
       * \note The main matrices _matrix_stiff, _matrix_mass, _matrix_diff, _matrix_react have to be
       *       allocated to their correct sizes beforehand
       */
      void init_symbolic()
      {
        XASSERTM(_matrix_mass_set && _matrix_conv_set, "Massmatrix and convection matrix must be set!");
        watch_init.start();

        // symbolic init, structure of matrices must not be changend after this

        // matrix structures for internal matrices
        _matrix_hlp = _matrix_mass.clone(LAFEM::CloneMode::Layout);
        _matrix_right = _matrix_mass.clone(LAFEM::CloneMode::Layout);
        _matrix_art_diff = _matrix_mass.clone(LAFEM::CloneMode::Layout);

        // only need matrix_left for implicit time stepping
        if(_theta > 0)
          _matrix_left = _matrix_mass.clone(LAFEM::CloneMode::Layout);

        // need those all the time to iterate over matrixes
        _num_rows = _matrix_mass.local().num_rows();

        // we need one to clone our help vectors from
        _vec_mass_lump = _matrix_mass.create_vector_l();
        //now we also need way more help vectors, 5 for explicit MCL, 8 for FCT and 7 for implicit MCL
        for(int i(0); i < 8; ++i) // everyone need at least 4 vecs
          _vec_hlp.push_back(_vec_mass_lump.clone(LAFEM::CloneMode::Layout));
        // also initialize vec_bc and set it to zero in case its not set by user
        _vec_bc = _vec_mass_lump.clone(LAFEM::CloneMode::Layout);
        _vec_bc.format();

        // the stuff for synchronising
        _row_gate = _matrix_mass.get_row_gate();
        _col_gate = _matrix_mass.get_col_gate();

        _synch_matrix_ticket = std::make_shared<Global::SynchMatrix<LocalMatrixType, typename Intern::MatrixTypeHelper<MatrixType_>::RowMirror>>(
          *_row_gate->_comm, _row_gate->_ranks, _row_gate->_mirrors, _col_gate->_mirrors);

        // matrix for 1 to 0 conversion of matrices
        _matrix_freq = _matrix_mass.local().clone(LAFEM::CloneMode::Layout);
        _matrix_freq.format(DataType(1.0));

        // sync it to know which matrix entry lives on how many patches
        _synch_matrix_ticket->init(_matrix_freq);
        _synch_matrix_ticket->exec(_matrix_freq);

        // solver for implicit time stepping
        if(_theta > 0)
        {
          _matrix_schwarz = _matrix_left.convert_to_1();

          _setup_lin_solver();

          _solver->set_plot_mode(Solver::PlotMode::none);
          _solver->set_tol_abs(1e-6);
          _solver->init_symbolic();
        }

        watch_init.stop();

        _init_symbolic_called = true;
      }

      /**
       * \brief Calculated the cfl condition
       *
       * Returns the cfl condition for the timestep d_t <= 1/(1 - theta) * m_i/(l_ii + k_ii - d_ii);
       * L = stiffnes matrix, K = Convection
       *
       * \note Will return 0 if theta = 1
       * \note Artificial diffusion needs to be assembled first
       */
      DataType cfl_cond()
      {
        if(_theta == 1)
          return 0;
        VectorType& vec_cfl = _vec_hlp.at(1);
        _get_cfl_vec(vec_cfl);
        // first we have to sync it (since its made from the type-0 matrices)
        // vec_cfl.sync_0(); // or we don't, since we call all operations on global matrices

        // vec_mass_lump needs to be synced here!!!
        vec_cfl.component_product(vec_cfl, _vec_mass_lump);

        return vec_cfl.min_abs_element();
      }

      /**
       * \brief Initializes the time stepizse
       *
       * Calculates the artificial diffusion and returns the maximum timestepsize according to the cfl condition
       *
       * \returns stepsize
       * maximum timestepsize
       *
       * \note
       * The convection, mass (and stiffness and reaction, if wanted) matrix must be set and numerically assembled
       */
      DataType init_stepsize()
      {
        XASSERTM(_init_symbolic_called, "Need to call init_symbolic() first");
        watch_init.start();

        // all this stuff is unnecessary for those matrices that didn't change after the last init_stepsize call...
        _matrix_sync_0_1_0(_matrix_mass);
        _matrix_sync_0_1_0(_matrix_conv);
        if(_matrix_stiff_set)
          _matrix_sync_0_1_0(_matrix_stiff);
        if(_matrix_react_set)
          _matrix_sync_0_1_0(_matrix_react);

        _vec_mass_lump = _matrix_mass.lump_rows();
        _calc_artificial_diffusion();
        _matrix_sync_0_1_0(_matrix_art_diff);

        DataType cfl = cfl_cond();

        _init_stepsize_called = true;
        watch_init.stop();

        return cfl;
      }

      /**
       * \brief Initializes the time stepping scheme
       *
       * Builds matrices, vectors and solver needed for a limiting step
       *
       * \param[in] stepsize
       * the chosen timestep size, must be greater 0 and smaller inf
       *
       * \note
       * init_stepsize must be called beforehand
       */
      void init_time_stepping(DataType stepsize)
      {
        XASSERTM(_init_stepsize_called, "Need to call init_stepsize() first");
        XASSERTM(stepsize > 0, "stepsize has to be greater than zero! If using implicit euler as time stepping scheme (theta = 1), you need to specify a stepsize yourself!");

        watch_init.start();

        _stepsize = stepsize;

        // assemble matrix_right fpr calculation of right hand side for low order solution
        // we start with R = K - D + L + N
        _matrix_right.local().copy(_matrix_conv.local());
        _matrix_right.local().axpy(_matrix_art_diff.local(), DataType(-1.0));
        if(_matrix_stiff_set)
        {
          // if(_bc_config == 2)
          // {
          //   _filter.local().filter_mat(_matrix_stiff.local());

          //   auto& vec_dof_fltr = _filter.local().get_filter_vector();
          //   auto* idx_fltr = vec_dof_fltr.indices();

          //   DataType* val_stiff = _matrix_stiff.local().val();
          //   for(Index i = 0; i < vec_dof_fltr.used_elements(); ++i)
          //   {
          //     const IndexType row = idx_fltr[i];
          //     for(IndexType j(_row_ptr[row]); j < _row_ptr[row + 1]; ++j)
          //       val_stiff[j] = DataType(0);
          //   }
          // }

          _matrix_right.local().axpy(_matrix_stiff.local());
        }
        if(_matrix_react_set)
        {
          // _filter.local().filter_mat(_matrix_react.local());
          _matrix_right.local().axpy(_matrix_react.local());
        }

        // then we scale according to our time scheme and limiter
        DataType scale_right = (_theta - DataType(1.0)) * (_limiter == 1 ? _stepsize : 1);
        if(_theta > 0)
          _matrix_left.local().scale(_matrix_right.local(), _stepsize * _theta);
        _matrix_right.local().scale(_matrix_right.local(), scale_right);

        if(_limiter == 1)
        {
          // we need to add the lmped mass matrix onto the other matrices,
          // therefore we need the vector to be type 0
          auto vec_mass_lump_type_0 = _matrix_mass.lump_rows(false);
          const Memory::TypedView<DataType> lump_vec_view(vec_mass_lump_type_0.local().elements_view_r());
          const Memory::TypedView<IndexType> row_ptr_view(_matrix_right.local().row_ptr_view_r());
          const Memory::TypedView<IndexType> col_idx_view(_matrix_right.local().col_idx_view_r());
          Memory::TypedView<DataType> val_right = _matrix_right.local().val_view_rw();
          Memory::TypedView<DataType> val_left = _matrix_left.local().val_view_rw();
          for(IndexType i = 0; i < _num_rows; ++i)
          {
            for(IndexType j = row_ptr_view[i]; j < row_ptr_view[i+1]; ++j)
            {
              if(col_idx_view[j] == i)
              {
                val_right[j] += lump_vec_view(i);
                if(_theta > 0)
                  val_left[j] += lump_vec_view(i);
              }
            }
          }
        }

        // _filter.local().filter_mat(_matrix_right.local());

        if(_theta > 0)
        {
          _matrix_schwarz = _matrix_left.convert_to_1();
          _solver->init_numeric();
        }
        else
        {
          // if theta = 0, we solve with the lumped mass matrix (which is diagonal)
          _vec_mass_lump_inverse = _vec_mass_lump.clone(LAFEM::CloneMode::Layout);
          _vec_mass_lump_inverse.component_invert(_vec_mass_lump);
        }

        _init_time_stepping_called = true;
        watch_init.stop();
      }

      /**
       * \brief Numeric initialization
       *
       * Calls init_stepsize and init_time_stepping
       *
       * \returns stepsize
       * chosen timestepsize
       *
       * \note
       * The convection and mass (and stiffness and reaction) matrix must be set and numerically assembled
       * \note
       * If you want to choose your own timestep size and not just use the cfl-condition (e.g. because you use implicit euler)
       * you must call init_stepsize and then init_time_stepping one after another
       */
      DataType init_numeric()
      {
        _stepsize = init_stepsize();
        init_time_stepping(_stepsize);
        return _stepsize;
      }

      /**
       * \brief Calls done_numeric() and done_symbolic()
       */
      void done()
      {
        this->done_numeric();
        this->done_symbolic();
      }

      /**
       * \brief Release any data allocated in the call to init_symbolic()
       */
      void done_symbolic()
      {
        XASSERTM(_init_symbolic_called, "Need to call init_symbolic() first");
        _matrix_hlp.local().clear();
        _matrix_right.local().clear();
        _matrix_art_diff.local().clear();
        if(_theta > 0)
          _matrix_left.local().clear();

        _vec_mass_lump.clear();
        _vec_bc.clear();
        for(auto& vec : _vec_hlp)
          vec.clear();
        _vec_hlp.clear();

        _matrix_freq.clear();

        if(_theta > 0)
          _solver->done_symbolic();

        _init_symbolic_called = false;
        _init_stepsize_called = false;
        _init_time_stepping_called = false;
      }

      /**
       * \brief Does nothing
       */
      void done_numeric()
      {
      }

      /**
       * \brief Applies transport solver to a given rhs
       */
      Solver::Status apply(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        XASSERTM(_init_time_stepping_called, "Need to call init_numeric() or init_time_stepping() first");

        watch_solve.start();

        if(_limiter == 1)
          limit_fct(vec_sol, vec_rhs);
        else if(_limiter == 2)
          limit_mcl_explicit(vec_sol, vec_rhs);

        watch_solve.stop();

        return Solver::Status::success;
      }

      /**
       * \brief Apply Heun
       */
      Solver::Status apply_heun(VectorType& vec_sol, const VectorType& vec_rhs1, const VectorType& vec_rhs2)
      {
        XASSERTM(_theta == 0, "Only use heun when theta = 0 (because Heun makes two explicit Euler steps)");
        watch_solve.start();
        auto vec_hlp = vec_sol.clone();
        watch_solve.stop();

        // Heun step:
        this->apply(vec_sol, vec_rhs1);
        this->apply(vec_sol, vec_rhs2);

        watch_solve.start();
        vec_sol.axpy(vec_hlp);
        vec_sol.scale(vec_sol, 0.5);
        watch_solve.stop();

        return Solver::Status::success;
      }

      /**
       * \brief Print measured time
       *
       * Prints the time for the init*() calls and the apply calls
       */
      void print_times()
      {
        _row_gate->_comm->print("Time initializing solver: " + watch_init.elapsed_string());
        _row_gate->_comm->print("Time solving:             " + watch_solve.elapsed_string());
      }

      /**
       * \brief Makes an FCT step
       *
       * \note Its recommended to just call the apply method
       */
      void limit_fct(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        if(_bc_config == 2)
          _filter.filter_sol(vec_sol);
        // solve for low order solution
        VectorType& vec_rhs_low = _vec_hlp.at(2);
        // if(_bc_config == 2)
        //   _filter.filter_rhs(vec_rhs_low);
        _matrix_right.apply(vec_rhs_low, vec_sol);
        if(_bc_config == 1)
          vec_rhs_low.axpy(_vec_bc, _stepsize);
        vec_rhs_low.axpy(vec_rhs, _stepsize);

        VectorType& vec_sol_low = _vec_hlp.at(3);
        vec_sol_low.format();
        if(_theta == 0)
          vec_sol_low.component_product(vec_rhs_low, _vec_mass_lump_inverse);
        else
        {
          if(_bc_config == 2)
          {
            _filter.filter_rhs(vec_rhs_low);
            _filter.filter_sol(vec_sol_low);
          }
          Solver::Status solving = Solver::solve(*_solver, vec_sol_low, vec_rhs_low, _matrix_left, _filter);

          if (solving != Solver::Status::success)
          {
            std::cerr << std::endl
                      << "ERROR: solving was not successful" << std::endl;
            Runtime::abort();
          }
        }

        // local bounds for each node
        VectorType& vec_min = _vec_hlp.at(4);
        VectorType& vec_max = _vec_hlp.at(5);
        _calc_local_bounds_fct(vec_min.local(), vec_max.local(), vec_sol_low.local());
        // which need a special synchronisation, where we don't add the values up, but look for the min/max of the bordering patches
        _synch_vector_max(vec_min.local(), *_row_gate->_comm, _row_gate->_ranks, _col_gate->_mirrors);
        _synch_vector_max(vec_max.local(), *_row_gate->_comm, _row_gate->_ranks, _col_gate->_mirrors);

        //for the sum of the negative/positive fluxes that go into each node
        VectorType& vec_neg_flux = _vec_hlp.at(6);
        VectorType& vec_pos_flux = _vec_hlp.at(7);
        // calculates above mentioned variables, matrix_hlp has the unlimited fluxes
        _calc_antidiffusive_flux_fct(_matrix_hlp.local(), vec_neg_flux.local(), vec_pos_flux.local(), vec_sol.local(), vec_sol_low.local());
        // the sum of the fluxes has to be synchronised
        vec_neg_flux.sync_0();
        vec_pos_flux.sync_0();

        // limits the antidiffusive fluxes and return vector, that has to be added to low order solution
        VectorType& vec_lim_flux = _vec_hlp.at(2);
        _calc_lim_flux_fct(vec_lim_flux.local(), vec_min.local(), vec_max.local(), vec_neg_flux.local(), vec_pos_flux.local(), _matrix_hlp.local(), _vec_mass_lump.local());
        // the limited flux has to be synced as well
        vec_lim_flux.sync_0();

        //update solution with limited flux
        vec_sol.copy(vec_lim_flux);
        vec_sol.axpy(vec_sol_low);
      }

      /**
       * \brief Makes an MCL step
       *
       * \note Its recommended to just call the apply method
       */
      void limit_mcl_explicit(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        VectorType& vec_lim_flux = _vec_hlp.at(4);
        _calc_lim_flux_mcl(vec_lim_flux, vec_sol);

        VectorType& vec_hlp = _vec_hlp.at(1);
        // add fluxes to rhs of ode form
        vec_hlp.copy(_vec_hlp.at(0));
        vec_hlp.axpy(vec_lim_flux);// in vec_hlp[0] ist rhs_low = R*u gespeichert
        // also add rhs
        vec_hlp.axpy(vec_rhs);
        // solve ode by multipling with inverse lumped mass matrix
        vec_hlp.component_product(vec_hlp, _vec_mass_lump_inverse);
        // scaling with time step
        vec_hlp.scale(vec_hlp, _stepsize);
        // and finally add to solution of the old time step (i.e. the explicit euler step)
        vec_sol.axpy(vec_hlp);
      }

    protected:
      // syncs the matrix from 0 to 1 to 0
      void _matrix_sync_0_1_0(MatrixType& mat)
      {
        auto& mat_local = mat.local();
        _synch_matrix_ticket->exec(mat_local);
        const Memory::TypedView<DataType> vfreq = _matrix_freq.val_view_r();
        Memory::TypedView<DataType> val = mat_local.val_view_rw();
        for(Index i = 0; i < mat_local.num_nzes(); ++i)
          val[i] /= vfreq[i];
      }

      void _setup_lin_solver()
      {
        // we need a preconditioner
        std::shared_ptr<Solver::SolverBase<LocalVectorType>> precond;
        if(_precond == 1)
          precond = Solver::new_ilu_precond(PreferredBackend::generic, _matrix_schwarz, _filter.local());
        else if(_precond == 2)
          precond = Solver::new_jacobi_precond(_matrix_schwarz, _filter.local());
        else if(_precond == 3)
          precond = Solver::new_ssor_precond(PreferredBackend::generic, _matrix_schwarz, _filter.local());
        else
          XABORTM("No valid preconditioner chosen");

        // and the schwarz stuff
        auto schwarz = Solver::new_schwarz_precond(precond, _filter);

        // to make our solver
        if(_lin_solver == 1)
          _solver = Solver::new_bicgstab(_matrix_left, _filter, schwarz);
        else if(_lin_solver == 2)
          _solver = Solver::new_fgmres(_matrix_left, _filter, 20, 0, schwarz);
        else
          XABORTM("No valid linear solver chosen");
      }

      /**
       * Calculated artifical diffusion of konvection matrix k by
       * d_{ij} = max(|k_{ij}|, 0, |k_{ji}|) for j != i;
       *          -sum_{l != i} d_{ij}       for i = j
       */
      void _calc_artificial_diffusion()
      {
        _matrix_art_diff.local().format();

        _matrix_right.local().copy(_matrix_conv.local());
        if(_limit_matrix_stiff)
        {
          XASSERTM(_matrix_stiff_set, "Stiffness matrix needs to be set to limit it");
          _matrix_right.local().axpy(_matrix_stiff.local());
        }
        _matrix_hlp.local().transpose(_matrix_right.local());

        const Memory::TypedView<IndexType> row_ptr_view = _matrix_mass.local().row_ptr_view_r();
        const Memory::TypedView<IndexType> col_idx_view = _matrix_mass.local().col_idx_view_r();
        const Memory::TypedView<DataType> val = _matrix_right.local().val_view_r();
        const Memory::TypedView<DataType> val_t = _matrix_hlp.local().val_view_r();
        Memory::TypedView<DataType> val_art_diff = _matrix_art_diff.local().val_view_rw();

        for(IndexType row = 0; row < _num_rows; ++row)
        {
          const IndexType start = row_ptr_view[row];
          const IndexType end = row_ptr_view[row + 1];

          IndexType diag = end;
          DataType lump = DataType(0);

          for(IndexType k = start; k < end; ++k)
          {
            if(col_idx_view[k] == row)
              diag = k;
            else
            {
              const DataType max_val_val_t = Math::max(Math::abs(val[k]), Math::abs(val_t[k]));
              lump += (val_art_diff[k] = Math::max(DataType(0), max_val_val_t));
            }
          }

          val_art_diff[diag] = -lump;
        }
      }

      // calculates the vector 1/(1 - theta) * 1/(e_ii + l_ii + k_ii - d_ii), needed for the cfl condition (e is react mat, l is stiffnes mat, k is conv mat, d is art diff mat)
      void _get_cfl_vec(VectorType& vec_cfl)
      {
        // was passiert bei theta = 1?
        VectorType& vec_hlp = _vec_hlp.at(0);

        // returns synched diagonal entries
        _matrix_conv.extract_diag(vec_cfl);
        if(_matrix_stiff_set)
        {
          _matrix_stiff.extract_diag(vec_hlp);
          vec_cfl.axpy(vec_hlp);
        }
        if(_matrix_react_set)
        {
          _matrix_react.extract_diag(vec_hlp);
          vec_cfl.axpy(vec_hlp);
        }
        _matrix_art_diff.extract_diag(vec_hlp);
        vec_cfl.axpy(vec_hlp, DataType(-1.0));

        vec_cfl.component_invert(vec_cfl, DataType(1.0/(1.0 - _theta)));
      }

      // ------------------- everything FCT -------------------
      // u_i^min = max_j(0, -(u_i^L - u_j^L)); u_i^max = max_j(0, (u_i^L - u_j^L))
      void _calc_local_bounds_fct(LocalVectorType& vec_u_min, LocalVectorType& vec_u_max, LocalVectorType& vec_sol_low)
      {
        XASSERT(vec_sol_low.size() == vec_u_min.size());
        XASSERT(vec_sol_low.size() == vec_u_max.size());

        vec_u_min.format(); //prob not necessary?
        vec_u_max.format();

        const Memory::TypedView<IndexType> row_ptr_view = _matrix_mass.local().row_ptr_view_r();
        const Memory::TypedView<IndexType> col_idx_view = _matrix_mass.local().col_idx_view_r();
        const Memory::TypedView<DataType> vec_sol_low_ptr = vec_sol_low.elements_view_r();
        Memory::TypedView<DataType> vec_u_min_ptr = vec_u_min.elements_view_w();
        Memory::TypedView<DataType> vec_u_max_ptr = vec_u_max.elements_view_w();

        for(IndexType row = 0; row < _num_rows; ++row)
        {
          const IndexType start = row_ptr_view[row];
          const IndexType end = row_ptr_view[row + 1];
          DataType mi = Math::max(0.0, -(vec_sol_low_ptr[row] - vec_sol_low_ptr[col_idx_view[start]]));
          DataType ma = Math::max(0.0, (vec_sol_low_ptr[row] - vec_sol_low_ptr[col_idx_view[start]]));
          for (IndexType i = start + 1; i < end; ++i)
          {
            mi = Math::max(mi, Math::max(0.0, -(vec_sol_low_ptr[row] - vec_sol_low_ptr[col_idx_view[i]])));
            ma = Math::max(ma, Math::max(0.0, (vec_sol_low_ptr[row] - vec_sol_low_ptr[col_idx_view[i]])));
          }
          vec_u_min_ptr[row] = mi;
          vec_u_max_ptr[row] = ma;
        }
      }

      void _calc_antidiffusive_flux_fct(LocalMatrixType& matrix_flux, LocalVectorType& vec_neg_flux,
        LocalVectorType& vec_pos_flux, LocalVectorType& vec_sol_corr, LocalVectorType& vec_sol_low)
      {
        // calc derivatve
        LocalVectorType& vec_deriv = _vec_hlp.at(0).local();
        vec_deriv.copy(vec_sol_low);
        vec_deriv.axpy(vec_sol_corr, -1.0);
        vec_deriv.scale(vec_deriv, 1.0/_stepsize);

        matrix_flux.format();
        vec_neg_flux.format();
        vec_pos_flux.format();


        const Memory::TypedView<IndexType> row_ptr_view = _matrix_mass.local().row_ptr_view_r();
        const Memory::TypedView<IndexType> col_idx_view = _matrix_mass.local().col_idx_view_r();
        const Memory::TypedView<DataType> val_mass = _matrix_mass.local().val_view_r();
        const Memory::TypedView<DataType> val_art_diff = _matrix_art_diff.local().val_view_r();
        const Memory::TypedView<DataType> val_deriv = vec_deriv.elements_view_r();
        const Memory::TypedView<DataType> val_sol_corr = vec_sol_corr.elements_view_r();
        const Memory::TypedView<DataType> val_sol_low = vec_sol_low.elements_view_r();
        Memory::TypedView<DataType> val_neg_flux = vec_neg_flux.elements_view_rw();
        Memory::TypedView<DataType> val_pos_flux = vec_pos_flux.elements_view_rw();
        Memory::TypedView<DataType> val_flux = matrix_flux.val_view_r();

        for(IndexType i = 0; i < _num_rows; ++i)
        {
          const IndexType start = row_ptr_view[i];
          const IndexType end = row_ptr_view[i+1];

          for(IndexType m = start; m < end; ++m)
          {
            const IndexType j = col_idx_view[m];

            // flux is calculated as in Christophs matlab code
            const DataType flux = _stepsize * val_mass[m] * (val_deriv(j) - val_deriv(i))
              - _stepsize * _theta * val_art_diff[m] * (val_sol_low(i) - val_sol_low(j))
              - _stepsize * (1.0 - _theta) * val_art_diff[m] * (val_sol_corr(i) - val_sol_corr(j));

            // prelimiting step
            if(flux * (val_sol_low(i) - val_sol_low(j)) > 0)
              val_flux[m] = 0;
            else
              val_flux[m] = flux;

            // sum of negative (sum_i!=j max(0, -f_ij)) and positive (sum_i!=j max(0, f_ij)) flux into one node
            if(j != i)
            {
              val_neg_flux[i] += Math::max(0.0, -val_flux[m]);
              val_pos_flux[i] += Math::max(0.0, val_flux[m]);
            }
          }
        }
      }

      // returns the limited fluxes of the fct-scheme, that have to be added to the low order solution
      void _calc_lim_flux_fct(LocalVectorType& vec_lim_flux, LocalVectorType& vec_min, LocalVectorType& vec_max,
        LocalVectorType& vec_neg_flux, LocalVectorType& vec_pos_flux, LocalMatrixType& matrix_flux, LocalVectorType& vec_mass_lump)
      {
        LocalVectorType& vec_bnd_neg = _vec_hlp.at(0).local();
        vec_bnd_neg.format();
        LocalVectorType& vec_bnd_pos = _vec_hlp.at(1).local();
        vec_bnd_pos.format();

        const Memory::TypedView<DataType> val_mass_lump = vec_mass_lump.elements_view_r();
        const Memory::TypedView<DataType> val_min = vec_min.elements_view_r();
        const Memory::TypedView<DataType> val_max = vec_max.elements_view_r();
        const Memory::TypedView<DataType> val_pos_flux = vec_pos_flux.elements_view_r();
        const Memory::TypedView<DataType> val_neg_flux = vec_neg_flux.elements_view_r();

        Memory::TypedView<DataType> vbn = vec_bnd_neg.elements_view_rw();
        Memory::TypedView<DataType> vbp = vec_bnd_pos.elements_view_rw();

        for(IndexType row = 0; row < _num_rows; ++row)
        {
          // distance from low order solution to local min/max, looks like this because of definition of local min/max
          DataType Q_p = val_mass_lump(row) * val_max(row);
          DataType Q_m = val_mass_lump(row) * val_min(row);

          // bounds for correction factors
          vbp[row] = Math::min(1.0, Q_p/val_pos_flux(row));
          vbn[row] = Math::min(1.0, Q_m/val_neg_flux(row));
        }

        if(_bc_config == 2)
        {
          // look for every boundary node and set the pos and neg bound to 1
          // (because then alpha_ij is also one)
          auto& vec_dof_fltr = _filter.local().get_filter_vector();
          const Memory::TypedView<IndexType> idx_fltr = vec_dof_fltr.indices_view_r();
          for(Index row = 0; row < vec_dof_fltr.num_nzes(); ++row)
            vbn[idx_fltr[row]] = DataType(1);
        }

        // calc correction factor and multiply together with unlim fluxes
        const Memory::TypedView<DataType> val_flux = matrix_flux.val_view_r();

        vec_lim_flux.format();
        Memory::TypedView<DataType> val_lim_flux = vec_lim_flux.elements_view_rw();

        const Memory::TypedView<IndexType> row_ptr_view = _matrix_mass.local().row_ptr_view_r();
        const Memory::TypedView<IndexType> col_idx_view = _matrix_mass.local().col_idx_view_r();
        for(IndexType i = 0; i < _num_rows; ++i)
        {
          const IndexType start = row_ptr_view[i];
          const IndexType end = row_ptr_view[i+1];

          // add up to get flux into one node
          DataType sum = 0;
          for(IndexType m = start; m < end; ++m)
          {
            const IndexType j = col_idx_view[m];
            if(i == j) continue;

            if(val_flux[m] > 0)
              sum += Math::min(vbp(i), vbn(j)) * val_flux[m];
            else if(val_flux[m] < 0)
              sum += Math::min(vbn(i), vbp(j)) * val_flux[m];
          }
          val_lim_flux[i] = -sum/val_mass_lump(i);
        }
      }

      // ------------------- everything MCL -------------------
      void _calc_lim_flux_mcl(VectorType& vec_lim_flux, VectorType& vec_sol)
      {
        VectorType& vec_rhs_low = _vec_hlp.at(0);
        // apply matrix to get right hand side of ode (low order)
        _matrix_right.apply(vec_rhs_low, vec_sol);
        // add type-1 bc vector
        if(_bc_config == 1)
          vec_rhs_low.axpy(_vec_bc);

        VectorType& vec_deriv = _vec_hlp.at(1);
        vec_deriv.component_product(vec_rhs_low, _vec_mass_lump_inverse);

        VectorType& vec_min = _vec_hlp.at(2);
        VectorType& vec_max = _vec_hlp.at(3);
        // local min/max for each node
        _calc_local_bounds_mcl(vec_min.local(), vec_max.local(), vec_sol.local());
        // need the same special sync as seen in fct as well
        _synch_vector_min(vec_min.local(), *_row_gate->_comm, _row_gate->_ranks, _col_gate->_mirrors);
        _synch_vector_max(vec_max.local(), *_row_gate->_comm, _row_gate->_ranks, _col_gate->_mirrors);

        // scaled bar states
        _calc_bar_states(_matrix_hlp.local(), vec_sol.local());

        // calculation of the limited antidiffusive flux
        _calc_lim_flux_mcl_from_bar_state(vec_lim_flux.local(), vec_sol.local(), vec_deriv.local(), vec_min.local(), vec_max.local(), _matrix_hlp.local());
        // which has to be synced, since its based on the local matrices (type 0)
        vec_lim_flux.sync_0();
      }

      // u_i^min = min_j u_j; u_i^max = max_j u_j
      void _calc_local_bounds_mcl(LocalVectorType& vec_u_min, LocalVectorType& vec_u_max, LocalVectorType& vec_sol)
      {
        XASSERT(vec_sol.size() == vec_u_min.size());
        XASSERT(vec_sol.size() == vec_u_max.size());

        vec_u_min.format(); //prob not necessary?
        vec_u_max.format();

        const Memory::TypedView<IndexType> row_ptr_view = _matrix_mass.local().row_ptr_view_r();
        const Memory::TypedView<IndexType> col_idx_view = _matrix_mass.local().col_idx_view_r();
        const Memory::TypedView<DataType> vec_sol_ptr = vec_sol.elements_view_r();
        Memory::TypedView<DataType> vec_u_min_ptr = vec_u_min.elements_view_w();
        Memory::TypedView<DataType> vec_u_max_ptr = vec_u_max.elements_view_w();

        for(IndexType row = 0; row < _num_rows; ++row)
        {
          const IndexType start = row_ptr_view[row];
          const IndexType end = row_ptr_view[row + 1];

          DataType mi = vec_sol_ptr[col_idx_view[start]];
          DataType ma = mi;
          for (IndexType i = start + 1; i < end; ++i)
          {
            Math::minimax(vec_sol_ptr[col_idx_view[i]], mi, ma);
          }
          vec_u_min_ptr[row] = mi;
          vec_u_max_ptr[row] = ma;
        }
      }

      // scaled bar states d_ij \bar{u_ij} (convex average of u_i and u_j for no diffusion)
      void _calc_bar_states(LocalMatrixType& matrix, LocalVectorType& vec_sol)
      {
        matrix.format();

        const Memory::TypedView<IndexType> row_ptr_view = _matrix_mass.local().row_ptr_view_r();
        const Memory::TypedView<IndexType> col_idx_view = _matrix_mass.local().col_idx_view_r();
        const Memory::TypedView<DataType> val_sol = vec_sol.elements_view_w();
        Memory::TypedView<DataType> val_art_dif = _matrix_art_diff.local().val_view_r();
        Memory::TypedView<DataType> val_conv = _matrix_conv.local().val_view_r();
        Memory::TypedView<DataType> val_mat = matrix.val_view_w();

        for(IndexType row = 0; row < _num_rows; ++row)
        {
          const IndexType start = row_ptr_view[row];
          const IndexType end = row_ptr_view[row+1];

          for(IndexType i = start; i < end; ++i)
          {
            const IndexType col = col_idx_view[i];
            val_mat[i] = val_art_dif[i] * (val_sol(col) + val_sol(row)) - val_conv[i] * (val_sol(col) - val_sol(row));
          }
        }
      }

      // computes the limited flux that has to be added to the right hand side of the ode for the mcl scheme
      void _calc_lim_flux_mcl_from_bar_state(LocalVectorType& vec_lim_flux, LocalVectorType& vec_sol, LocalVectorType& vec_deriv, LocalVectorType& vec_min, LocalVectorType& vec_max, LocalMatrixType& matrix_bar_state)
      {
        vec_lim_flux.format();
        _matrix_hlp_t.local() = matrix_bar_state.transpose();

        const Memory::TypedView<IndexType> row_ptr_view = _matrix_mass.local().row_ptr_view_r();
        const Memory::TypedView<IndexType> col_idx_view = _matrix_mass.local().col_idx_view_r();
        const Memory::TypedView<DataType> val_sol = vec_sol.elements_view_r();
        const Memory::TypedView<DataType> val_deriv = vec_deriv.elements_view_r();
        const Memory::TypedView<DataType> val_min = vec_min.elements_view_r();
        const Memory::TypedView<DataType> val_max = vec_max.elements_view_r();
        Memory::TypedView<DataType> val_mass = _matrix_mass.local().val_view_r();
        Memory::TypedView<DataType> val_art_diff = _matrix_art_diff.local().val_view_r();
        Memory::TypedView<DataType> val_bar_state = matrix_bar_state.val_view_r();
        Memory::TypedView<DataType> val_bar_state_t = _matrix_hlp_t.local().val_view_r();
        Memory::TypedView<DataType> vec_lim_flux_ptr = vec_lim_flux.elements_view_w();

        for(IndexType row = 0; row < _num_rows; ++row)
        {
          const IndexType start = row_ptr_view[row];
          const IndexType end = row_ptr_view[row+1];

          DataType sum = 0;
          for(IndexType i = start; i < end; ++i)
          {
            const IndexType col = col_idx_view[i];

            if( row == col )
              continue;

            // the flux between node "row" and "col"
            const DataType flux = val_art_diff[i] * (val_sol(row) - val_sol(col))
              + val_mass[i] * (val_deriv(row) - val_deriv(col));
            // add prelimiting?

            // the limited flux
            if(flux > 0)
              sum += Math::min(flux, Math::min(2 * val_art_diff[i] * val_max(row) - val_bar_state[i], val_bar_state_t[i] - 2 * val_art_diff[i] * val_min(col)));
            else
              sum += Math::max(flux, Math::max(2 * val_art_diff[i] * val_min(row) - val_bar_state[i], val_bar_state_t[i] - 2 * val_art_diff[i] * val_max(col)));
          }
          vec_lim_flux_ptr[row] = sum;
        }
      }

      // ------------------- syncs for both limiters -------------------
      // synchronise a vector by not summing up over the values on each patch, but by taking the minimum over all patches
      void _synch_vector_min(FEAT::LAFEM::DenseVector<DataType, IndexType>& vector, const FEAT::Dist::Comm& comm,
        const std::vector<int>& ranks, const std::vector<FEAT::LAFEM::VectorMirror<DataType, IndexType>>& mirrors)
      {
        typedef LAFEM::DenseVector<DataType, IndexType> BufferMain;

        /// send and receive request vectors
        Dist::RequestVector send_reqs, recv_reqs;
        /// send and receive buffers
        std::vector<BufferMain> send_bufs, recv_bufs;
        std::vector<Memory::TypedView<DataType>> send_views, recv_views;

        const std::size_t n = ranks.size();

        XASSERTM(mirrors.size() == n, "invalid vector mirror count");

        // post receives
        recv_reqs.reserve(n);
        recv_bufs.resize(n);
        recv_views.resize(n);
        for(std::size_t i(0); i < n; ++i)
        {
          // create buffer vector
          recv_bufs.at(i) = BufferMain(mirrors.at(i).buffer_size(vector));
          recv_views.at(i) = recv_bufs.at(i).elements_view_w();

          // post receive
          recv_reqs.push_back(comm.irecv(recv_views.at(i).get_w(), recv_bufs.at(i).size(), ranks.at(i)));
        }

        // post sends
        send_reqs.reserve(n);
        send_bufs.resize(n);
        send_views.resize(n);
        for(std::size_t i(0); i < n; ++i)
        {
          // create buffer vector
          send_bufs.at(i) = BufferMain(mirrors.at(i).buffer_size(vector));

          // gather from mirror
          mirrors.at(i).gather(send_bufs.at(i), vector);
          send_views.at(i) = send_bufs.at(i).elements_view_r();

          // post send
          send_reqs.push_back(comm.isend(send_views.at(i).get_r(), send_bufs.at(i).size(), ranks.at(i)));
        }

        Memory::TypedView<DataType> values = vector.elements_view_rw();

        // process all pending receives
        for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
        {
          // get buffer indices and values
          const Index m = recv_bufs.at(idx).size();
          const Memory::TypedView<IndexType> ind = mirrors.at(idx).indices_view_r();
          //const DataType* val = recv_bufs.at(idx).elements();
          const Memory::TypedView<DataType>& val = recv_views.at(idx);
          for(Index k(0); k < m; ++k)
            //Math::mini(values(idx[k]), val[k]);
          values[ind[k]] = Math::min(values[ind[k]], val[k]);
        }

        // wait for all sends to finish
        send_reqs.wait_all();
      }

      // synchronise a vector by not summing up over the values on each patch, but by taking the maximum over all patches
      void _synch_vector_max(FEAT::LAFEM::DenseVector<DataType, IndexType>& vector, const FEAT::Dist::Comm& comm,
        const std::vector<int>& ranks, const std::vector<FEAT::LAFEM::VectorMirror<DataType, IndexType>>& mirrors)
      {
        typedef LAFEM::DenseVector<DataType, IndexType> BufferMain;

        /// send and receive request vectors
        Dist::RequestVector send_reqs, recv_reqs;
        /// send and receive buffers
        std::vector<BufferMain> send_bufs, recv_bufs;
        std::vector<Memory::TypedView<DataType>> send_views, recv_views;

        const std::size_t n = ranks.size();

        XASSERTM(mirrors.size() == n, "invalid vector mirror count");

        // post receives
        recv_reqs.reserve(n);
        recv_bufs.resize(n);
        recv_views.resize(n);
        for(std::size_t i(0); i < n; ++i)
        {
          // create buffer vector
          recv_bufs.at(i) = BufferMain(mirrors.at(i).buffer_size(vector));
          recv_views.at(i) = recv_bufs.at(i).elements_view_w();

          // post receive
          recv_reqs.push_back(comm.irecv(recv_views.at(i).get_w(), recv_bufs.at(i).size(), ranks.at(i)));
        }

        // post sends
        send_reqs.reserve(n);
        send_bufs.resize(n);
        send_views.resize(n);
        for(std::size_t i(0); i < n; ++i)
        {
          // create buffer vector
          send_bufs.at(i) = BufferMain(mirrors.at(i).buffer_size(vector));

          // gather from mirror
          mirrors.at(i).gather(send_bufs.at(i), vector);
          send_views.at(i) = send_bufs.at(i).elements_view_r();

          // post send
          send_reqs.push_back(comm.isend(send_views.at(i).get_r(), send_bufs.at(i).size(), ranks.at(i)));
        }

        Memory::TypedView<DataType> values = vector.elements_view_rw();

        // process all pending receives
        for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
        {
          // get buffer indices and values
          const Index m = recv_bufs.at(idx).size();
          const Memory::TypedView<IndexType> ind = mirrors.at(idx).indices_view_r();
          //const Memory::TypedView<DataType> val = recv_bufs.at(idx).elements_view_r();
          const Memory::TypedView<DataType>& val = recv_views.at(idx);
          for(Index k(0); k < m; ++k)
            values[ind[k]] = Math::max(values[ind[k]], val[k]);
        }

        // wait for all sends to finish
        send_reqs.wait_all();
      }
    };
  } // namespace Control
} // namespace FEAT
