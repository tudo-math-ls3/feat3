#pragma once
#ifndef KERNEL_LAFEM_PROTO_SOLVER_HPP
#define KERNEL_LAFEM_PROTO_SOLVER_HPP 1

// includes, FEAST
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/umfpack.hpp>
#include <memory>
#include <utility>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Solver status return codes enumeration
     *
     * This enumeration defined the solver status return codes, which specify whether
     * the solver application was successful or failed due to some reason.
     */
    enum class SolverStatus
    {
      /// undefined status
      undefined = 0,
      /// continue iteration (internal use only)
      progress,
      /// solving successful (convergence criterion fulfilled)
      success,
      /// premature abort (solver aborted due to internal errors or preconditioner failure)
      aborted,
      /// solver diverged (divergence criterion fulfilled)
      diverged,
      /// solver reached maximum iterations
      max_iter,
      /// solver stagnated (stagnation criterion fulfilled)
      stagnated
    };

    /**
     * \brief Status success check function
     *
     * This function takes a SolverStatus value as input and checks whether it respresents a
     * 'successul' run. A solving run is interpreted as successful, if one of the following
     * status codes was returned:
     *
     *  - SolverStatus::success
     *  - SolverStatus::max_iter
     *  - SolverStatus::stagnated
     *
     * For any other status code, the solving run is interpreted as unsuccessful.
     *
     * \param[in] status
     * A status code returned by a solver.
     *
     * \returns
     * \c true, if the run was successful, otherwise \c fase.
     */
    inline bool status_success(SolverStatus status)
    {
      switch(status)
      {
      case SolverStatus::success:
      case SolverStatus::max_iter:
      case SolverStatus::stagnated:
        return true;

      default:
        return false;
      }
    }

    /**
     * \brief Polymorphic provisional solver interface
     *
     * This class template defines the interface for any solver implementation.
     *
     * \tparam IOVector_
     * The class of the vector that is passed to the solver in the \c solve() method.
     *
     * \author Peter Zajac
     */
    template<typename IOVector_>
    class SolverInterface
    {
    public:
      /// virtual destructor
      virtual ~SolverInterface()
      {
      }

      /**
       * \brief Symbolic initialisation method
       *
       * This method is called to perform symbolic initialisation of the solver.
       *
       * \returns
       * \c true, if the initialisation was successful, otherwise \c false.
       */
      virtual bool init_symbolic()
      {
        return true;
      }

      /**
       * \brief Numeric initialisation method
       *
       * This method is called to perform numeric initialisation of the solver.\n
       * Before this function can be called, the symbolic initialisation must be performed.
       *
       * \returns
       * \c true, if the initialisation was successful, otherwise \c false.
       */
      virtual bool init_numeric()
      {
        return true;
      }

      /**
       * \brief Numeric finalisation method
       *
       * This method is called to release any data allocated in the numeric initialisation step.
       */
      virtual void done_numeric()
      {
      }

      /**
       * \brief Symbolic finalisation method
       *
       * This method is called to release any data allocated in the symbolic initialisation step.
       */
      virtual void done_symbolic()
      {
      }

      /**
       * \brief Initialisation method
       *
       * This function performs both the symbolic and numeric initialisation.
       */
      virtual bool init()
      {
        if(!init_symbolic())
          return false;
        if(!init_numeric())
          return false;
        return true;
      }

      /**
       * \brief Finalisation method
       *
       * This function performs both the symbolic and numeric finalisation.
       */
      virtual void done()
      {
        done_numeric();
        done_symbolic();
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the solver.
       */
      virtual String name()
      {
        return "SolverInterface";
      }

      /**
       * \brief Solver application method
       *
       * This method applies the solver represented by this object onto a given right-hand-side vector
       * and returns the corresponding solution vector.
       *
       * \note Solvers which derive from the IterativeSolver base class also provide a \c correct() method
       * which corrects an initial solution instead of starting with the null vector.
       *
       * \param[in,out] vec_sol
       * The vector that shall receive the solution of the linear system. It is assumed to be allocated, but
       * its numerical contents may be undefined upon calling this method.
       *
       * \param[in] vec_rhs
       * The vector that represents the right-hand-side of the linear system to be solved.
       *
       * \returns
       * A SolverStatus code that represents the status of the solution step.
       */
      virtual SolverStatus solve(IOVector_& vec_sol, const IOVector_& vec_rhs) = 0;
    }; // class SolverInterface<...>

    /**
     * \brief LAFEM Preconditioner Wrapper class template
     *
     * This class template acts as a wrapper around the preconditioners implemented in the <c>preconditioner.hpp</c>
     * header file.
     *
     * \tparam Algo_
     * The algorithm tag; is passed as the first parameter to the preconditioner class template.
     *
     * \tparam Matrix_
     * The matrix class; is passed as the second parameter to the preconditioner class template.
     *
     * \tparam Precond_
     * The preconditioner class template.
     *
     * <b>Example</b>:\n
     * To use the ILUPreconditioner class for CSR-matrices, one would have to use the following class template
     * combination:
     * <c>PreconWrapper<Algo::Generic, SparseMatrixCSR<Mem::Main,double>, ILUPreconditioner></c>.
     *
     * \author Peter Zajac
     */
    template<
      typename Algo_,
      typename Matrix_,
      template<typename,typename,typename> class Precon_>
    class PreconWrapper :
      public SolverInterface<typename Matrix_::VectorTypeR>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef typename MatrixType::VectorTypeR VectorType;

    protected:
      /// the actual preconditioner object
      Precon_<Algo_, MatrixType, VectorType> _precond;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] args
       * The remaining arguments which are passed to the preconditioner. For the required set of arguments,
       * see the documentation of the corresponding preconditioner class template.
       */
      template<typename... Args_>
      explicit PreconWrapper(Args_&&... args) :
        _precond(std::forward<Args_>(args)...)
      {
      }

      /// Applies the preconditioner.
      virtual SolverStatus solve(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        _precond.apply(vec_sol, vec_rhs);
        return SolverStatus::success;
      }
    }; // class PreconWrapper<...>

#ifdef FEAST_HAVE_UMFPACK
    /**
     *\brief Umfpack Solver class
     *
     * This class wraps around the \c Umfpack class implemented in the <c>umfpack.hpp</c> header.
     *
     * \author Peter Zajac
     */
    class UmfpackSolver :
      public SolverInterface<DenseVector<Mem::Main, double, Index>>
    {
    public:
      typedef SparseMatrixCSR<Mem::Main, double, Index> MatrixType;
      typedef DenseVector<Mem::Main, double, Index> VectorType;

    protected:
      const MatrixType& _matrix;
      Umfpack _umfpack;

    public:
      explicit UmfpackSolver(const MatrixType& matrix) :
        _matrix(matrix),
        _umfpack()
        {
        }

      virtual String name() override
      {
        return "Umfpack";
      }

      virtual bool init_symbolic() override
      {
        try
        {
          _umfpack.init_symbolic(&_matrix, true);
        }
        catch(...)
        {
          return false;
        }
        return true;
      }

      virtual bool init_numeric() override
      {
        try
        {
          _umfpack.init_numeric();
        }
        catch(...)
        {
          return false;
        }
        return true;
      }

      virtual void done_numeric() override
      {
        _umfpack.free_numeric();
      }

      virtual void done_symbolic() override
      {
        _umfpack.free_symbolic();
      }

      virtual SolverStatus solve(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        try
        {
          _umfpack.solve(vec_sol, vec_rhs);
        }
        catch(...)
        {
          return SolverStatus::aborted;
        }
        return SolverStatus::success;
      }
    }; // class UmfpackSolver
#endif // FEAST_HAVE_UMFPACK

    /**
     * \brief Abstract base-class for iterative solvers.
     *
     * This class template acts as an abstract base class for iterative solvers.
     * It also implements various auxiliary features for convergence control.
     *
     * \tparam AlgoType_
     * The algorithm tag to be used by the solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Peter Zajac
     */
    template<typename AlgoType_, typename Matrix_, typename Filter_>
    class IterativeSolver :
      public SolverInterface<typename Matrix_::VectorTypeR>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename Matrix_::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef SolverInterface<VectorType> BaseClass;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// name of the solver in plots
      String _plot_name;
      /// relative tolerance parameter
      DataType _tol_rel;
      /// relative tolerance parameter
      DataType _tol_abs;
      /// relative divergence parameter
      DataType _div_rel;
      /// absolute divergence parameter
      DataType _div_abs;
      /// minimum number of iterations
      Index _min_iter;
      /// maximum number of iterations
      Index _max_iter;
      /// number of performed iterations
      Index _num_iter;
      /// initial defect
      DataType _def_init;
      /// current defect
      DataType _def_cur;
      /// iteration count digits for plotting
      Index _iter_digits;
      /// whether to plot something
      bool _plot;

      /**
       * \brief Protected constructor
       *
       * This constructor initialises the following values:
       *
       * - relative tolerance: sqrt(eps) (~1E-8 for double)
       * - absolute tolerance: 1/eps^2 (~1E+32 for double)
       * - relative divergence: 1/eps (~1E+16 for double)
       * - absolute divergence: 1/eps^2 (~1E+32 for double)
       * - minimum iterations: 0
       * - maximum iterations: 100
       * - convergence plot: false
       *
       * \param[in] name
       * Specifies the name of the iterative solver. This is used as a prefix for the convergence plot.
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       */
      explicit IterativeSolver(String plot_name, const MatrixType& matrix, const FilterType& filter) :
        BaseClass(),
        _system_matrix(matrix),
        _system_filter(filter),
        _plot_name(plot_name),
        _tol_rel(Math::sqrt(Math::eps<DataType>())),
        _tol_abs(DataType(1) / Math::sqr(Math::eps<DataType>())),
        _div_rel(DataType(1) / Math::eps<DataType>()),
        _div_abs(DataType(1) / Math::sqr(Math::eps<DataType>())),
        _min_iter(0),
        _max_iter(100),
        _num_iter(0),
        _def_init(0),
        _def_cur(0),
        _iter_digits(4),
        _plot(false)
      {
      }

    public:
      /// Sets the relative tolerance for the solver.
      void set_tol_rel(DataType tol_rel)
      {
        _tol_rel = tol_rel;
      }

      /// Sets the absolute tolerance for the solver.
      void set_tol_abs(DataType tol_abs)
      {
        _tol_abs = tol_abs;
      }

      /// Returns the relative tolerance.
      DataType get_tol_rel() const
      {
        return _tol_rel;
      }

      /// Returns the absolute tolerance.
      DataType get_tol_abs() const
      {
        return _tol_abs;
      }

      /// Sets the relative divergence for the solver.
      void set_div_rel(DataType div_rel)
      {
        _div_rel = div_rel;
      }

      /// Sets the absolute divergence for the solver.
      void set_div_abs(DataType div_abs)
      {
        _div_abs = div_abs;
      }

      /// Returns the relative divergence.
      DataType get_div_rel() const
      {
        return _div_rel;
      }

      /// Returns the absolute divergence.
      DataType get_div_abs() const
      {
        return _div_abs;
      }

      /// Sets the minimum iteration count for the solver.
      void set_min_iter(Index min_iter)
      {
        _min_iter = min_iter;
      }

      /// Sets the maximum iteration count for the solver.
      void set_max_iter(Index max_iter)
      {
        _max_iter = max_iter;
      }

      /// Returns number of performed iterations
      Index get_num_iter() const
      {
        return _num_iter;
      }

      /// Returns the minimal number of iterations
      Index get_min_iter() const
      {
        return _min_iter;
      }

      /// Returns the maximum number of iterations
      Index get_max_iter() const
      {
        return _max_iter;
      }

      /**
       * \brief Sets the plot mode of the solver.
       *
       * \param[in] plot
       * If set to \c true, the solver will print a convergence plot to std::cout.
       */
      void set_plot(bool plot)
      {
        _plot = plot;
      }

      /// checks for convergence
      bool is_converged() const
      {
        return (_def_cur <= _tol_abs) && (_def_cur <= (_tol_rel * _def_init));
      }

      /// checks for divergence
      bool is_diverged() const
      {
        return (_def_cur > _div_abs) || (_def_cur > (_div_rel * _def_init));
      }

      /// Returns the initial defect
      DataType get_def_initial() const
      {
        return _def_init;
      }

      /// Returns the final defect
      DataType get_def_final() const
      {
        return _def_cur;
      }

      /// Returns the overall convergence rate.
      DataType get_conv_rate() const
      {
        // no iterations performed?
        if(_num_iter <= Index(0))
          return DataType(0);
        // initial defect zero?
        if(_def_init < Math::eps<DataType>())
          return DataType(0);

        // compute convergence rate: (def_final / def_initial) ^ (1 / #iter)
        return Math::pow(_def_cur / _def_init, DataType(1) / DataType(_num_iter));
      }

      /**
       * \brief Solver correction method
       *
       * This method applies the solver represented by this object onto a given right-hand-side vector
       * and updates the corresponding solution vector.
       *
       * In contrast to the solve() method of the SolverInterface base class, this method uses the
       * vector \p vec_sol as the initial solution vector for the iterative solution process instead of
       * ignoring its contents upon entry and strating with the null vector.
       *
       * \note The default implementation creates a temporary vector, computes the defect, applies
       * the solve() method to obtain a correction vector and update the solution vector with it.
       *
       * \param[in,out] vec_sol
       * The vector that contains the initial solution upon entry and receives the solution
       * of the linear system upon exit.
       *
       * \param[in] vec_rhs
       * The vector that represents the right-hand-side of the linear system to be solved.
       *
       * \returns
       * A SolverStatus code that represents the status of the solution step.
       */
      virtual SolverStatus correct(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        // create a defect and a correction vector
        auto vec_def = vec_rhs.clone(LAFEM::CloneMode::Layout);
        auto vec_cor = vec_sol.clone(LAFEM::CloneMode::Layout);

        // compute defect
        this->_system_matrix.template apply<AlgoType>(vec_def, vec_sol, vec_rhs, -DataType(1));

        // apply defect filter
        this->_system_filter.template filter_def<AlgoType>(vec_def);

        // apply solver
        SolverStatus status = this->solve(vec_cor, vec_def);

        // apply correction if successful
        if(status_success(status))
        {
          // apply correction filter
          this->_system_filter.template filter_cor<AlgoType>(vec_cor);

          // update solution vector
          vec_sol.template axpy<AlgoType>(vec_cor, vec_sol);
        }

        // return status
        return status;
      }

    protected:
      /**
       * \brief Internal function: sets the initial defect vector
       *
       * \param[in] vector
       * The initial defect vector.
       *
       * \returns
       * A SolverStatus code.
       */
      SolverStatus _set_initial_defect(const VectorType& vector)
      {
        // store new defect
        this->_def_init = this->_def_cur = vector.template norm2<AlgoType>();
        this->_num_iter = Index(0);

        // plot?
        if(this->_plot)
          std::cout << this->_plot_name << ": " << stringify(0).pad_front(_iter_digits)
                    << " : " << scientify(this->_def_init) << std::endl;

        // continue iterating
        return SolverStatus::progress;
      }

      /**
       * \brief Internal function: sets the new (next) defect vector
       *
       * This function computes the defect vector's norm, increments the iteration count,
       * plots an output line to std::cout and checks whether any of the stopping criterions is fulfilled.
       *
       * \param[in] vector
       * The new defect vector.
       *
       * \returns
       * A SolverStatus code.
       */
      SolverStatus _set_new_defect(const VectorType& vector)
      {
        // increase iteration count
        ++this->_num_iter;

        // first, let's see if we have to compute the defect at all
        bool calc_def = false;
        calc_def = calc_def || (this->_min_iter < this->_max_iter);
        calc_def = calc_def || this->_plot;

        // compute new defect
        if(calc_def)
          this->_def_cur = vector.template norm2<AlgoType>();

        // plot?
        if(this->_plot)
          std::cout << this->_plot_name << ": " << stringify(this->_num_iter).pad_front(_iter_digits)
                    << " : " << scientify(this->_def_cur) << std::endl;

        // is diverged?
        if(is_diverged())
          return SolverStatus::diverged;

        // minimum number of iterations performed?
        if(this->_num_iter < this->_min_iter)
          return SolverStatus::progress;

        // is converged?
        if(is_converged())
          return SolverStatus::success;

        // maximum number of iterations performed?
        if(this->_num_iter >= this->_max_iter)
          return SolverStatus::max_iter;

        // continue iterating
        return SolverStatus::progress;
      }

      /**
       * \brief Internal function: sets the new (next) defect vector norm
       *
       * This function stores the defect vector's norm, increments the iteration count,
       * plots an output line to std::cout. This function is used by solvers with intermediate
       * iterations, such as FGMRES.
       *
       * \param[in] def_norm
       * The new defect vector norm.
       */
      void _set_new_defect_norm(DataType def_norm)
      {
        // increase iteration count
        ++this->_num_iter;

        // compute new defect
        this->_def_cur = def_norm;

        // plot?
        if(this->_plot)
          std::cout << this->_plot_name << "* " << stringify(this->_num_iter).pad_front(_iter_digits)
                    << " : " << scientify(this->_def_cur) << std::endl;
      }
    }; // class IterativeSolver

    /**
     * \brief Abstract base-class for preconditioned iterative solvers.
     *
     * This class extends the functionality of the IterativeSolver class template by providing overrides
     * for the initialisation and finalisation methods of the SolverInterface class template, which take
     * care of forwarding these steps to the preconditioner.
     *
     * \tparam AlgoType_
     * The algorithm tag to be used by the solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \authro Peter Zajac
     */
    template<typename AlgoType_, typename Matrix_, typename Filter_>
    class PreconditionedIterativeSolver :
      public IterativeSolver<AlgoType_, Matrix_, Filter_>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef IterativeSolver<AlgoType_, MatrixType, FilterType> BaseClass;

      /// the interface for the preconditioner
      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// the pointer to the preconditioner
      PrecondType* _precond;
      /// specifies whether to delete the preconditioner upon destruction
      bool _del_precond;

      /**
       * \brief Constructor
       *
       * \param[in] name
       * The name of the solver.
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be nullptr.
       *
       * \param[in] del_precond
       * Specifies whether the preconditioner object should be deleted upund
       * destruction of this solver object.
       */
      explicit PreconditionedIterativeSolver(String plot_name, const MatrixType& matrix, const FilterType& filter,
        PrecondType* precond = nullptr, bool del_precond = false) :
        BaseClass(plot_name, matrix, filter),
        _precond(precond),
        _del_precond(del_precond)
      {
      }

      /// virtual destructor
      virtual ~PreconditionedIterativeSolver()
      {
        if((_precond != nullptr) && _del_precond)
          delete _precond;
      }

    public:
      virtual bool init_symbolic() override
      {
        if(!BaseClass::init_symbolic())
          return false;
        if(_precond && !(_precond->init_symbolic()))
          return false;
        return true;
      }

      virtual bool init_numeric() override
      {
        if(!BaseClass::init_numeric())
          return false;
        if(_precond && !(_precond->init_numeric()))
          return false;
        return true;
      }

      virtual void done_numeric() override
      {
        if(_precond)
          _precond->done_numeric();
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        if(_precond)
          _precond->done_symbolic();
        BaseClass::done_symbolic();
      }

    protected:
      /**
       * \brief Applies the preconditioner onto a defect vector.
       *
       * \note
       * If no preconditioner is present, this function will simply copy the input vector's
       * contents into the output vector, therefore emulating an "identity preconditioner".
       *
       * \param[in,out] vec_out
       * A reference to the vector that shall receive the preconditioned defect.
       *
       * \param[in] vec_in
       * A reference to the vector that is to be preconditioned.
       *
       * \returns
       * \c true, if the preconditioner application was successful, otherwise \c false.
       */
      bool _apply_precond(VectorType& vec_out, const VectorType& vec_in)
      {
        if(this->_precond)
        {
          return status_success(this->_precond->solve(vec_out, vec_in));
        }
        else
        {
          vec_out.copy(vec_in);
          return true;
        }
      }
    }; // class PreconditionedIterativeSolver<...>

    /**
     * \brief Fix-Point (Richardson) solver implementation
     *
     * This class implements a simple Fix-Point iteration solver, a.k.a. Richardson.
     *
     * \tparam AlgoType_
     * The algorithm tag to be used by the solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Peter Zajac
     */
    template<typename AlgoType_, typename Matrix_, typename Filter_>
    class FixPointSolver :
      public PreconditionedIterativeSolver<AlgoType_, Matrix_, Filter_>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<AlgoType, MatrixType, FilterType> BaseClass;

      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// defect vector
      VectorType _vec_def;
      /// correction vector
      VectorType _vec_cor;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       *
       * \param[in] del_precond
       * Specifies whether the preconditioner object should be deleted upund
       * destruction of this solver object.
       */
      explicit FixPointSolver(const MatrixType& matrix, const FilterType& filter,
        PrecondType* precond = nullptr, bool del_precond = false) :
        BaseClass("FixPoint", matrix, filter, precond, del_precond)
      {
      }

      virtual bool init_symbolic() override
      {
        if(!BaseClass::init_symbolic())
          return false;
        // create three temporary vectors
        _vec_def = this->_system_matrix.create_vector_r();
        _vec_cor = this->_system_matrix.create_vector_r();
        return true;
      }

      virtual void done_symbolic() override
      {
        this->_vec_cor.clear();
        this->_vec_def.clear();
        BaseClass::done_symbolic();
      }

      virtual SolverStatus solve(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // save defect
        this->_vec_def.copy(vec_rhs);
        this->_system_filter.template filter_def<AlgoType>(this->_vec_def);

        // clear solution vector
        vec_sol.format();

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

      virtual SolverStatus correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.template apply<AlgoType>(this->_vec_def, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.template filter_def<AlgoType>(this->_vec_def);

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

    protected:
      virtual SolverStatus _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        VectorType& vec_def(this->_vec_def);
        VectorType& vec_cor(this->_vec_cor);
        const MatrixType& mat_sys(this->_system_matrix);
        const FilterType& fil_sys(this->_system_filter);
        SolverStatus status(SolverStatus::progress);

        // compute initial defect
        status = this->_set_initial_defect(vec_def);
        if(status != SolverStatus::progress)
          return status;

        // start iterating
        while(status == SolverStatus::progress)
        {
          // apply preconditioner
          if(!this->_apply_precond(vec_cor, vec_def))
            return SolverStatus::aborted;
          fil_sys.template filter_cor<AlgoType>(vec_cor);

          // update solution vector
          vec_sol.template axpy<AlgoType>(vec_cor, vec_sol);

          // compute new defect vector
          mat_sys.template apply<AlgoType>(vec_def, vec_sol, vec_rhs, -DataType(1));
          fil_sys.template filter_def<AlgoType>(vec_def);

          // compute new defect norm
          status = this->_set_new_defect(vec_def);
        }

        // return our status
        return status;
      }
    }; // class FixPointSolver<...>

    /**
     * \brief (Preconditioned) Conjugate-Gradient solver implementation
     *
     * This class implements a simple PCG solver.
     *
     * \tparam AlgoType_
     * The algorithm tag to be used by the solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Peter Zajac
     */
    template<
      typename AlgoType_,
      typename Matrix_,
      typename Filter_>
    class PCGSolver :
      public PreconditionedIterativeSolver<AlgoType_, Matrix_, Filter_>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<AlgoType, MatrixType, FilterType> BaseClass;

      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// defect vector
      VectorType _vec_def;
      /// descend direction vector
      VectorType _vec_dir;
      /// temporary vector
      VectorType _vec_tmp;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       *
       * \param[in] del_precond
       * Specifies whether the preconditioner object should be deleted upund
       * destruction of this solver object.
       */
      explicit PCGSolver(const MatrixType& matrix, const FilterType& filter,
        PrecondType* precond = nullptr, bool del_precond = false) :
        BaseClass("PCG", matrix, filter, precond, del_precond)
      {
      }

      virtual String name() override
      {
        return "PCG";
      }

      virtual bool init_symbolic() override
      {
        if(!BaseClass::init_symbolic())
          return false;
        // create three temporary vectors
        _vec_def = this->_system_matrix.create_vector_r();
        _vec_dir = this->_system_matrix.create_vector_r();
        _vec_tmp = this->_system_matrix.create_vector_r();
        return true;
      }

      virtual void done_symbolic() override
      {
        this->_vec_tmp.clear();
        this->_vec_dir.clear();
        this->_vec_def.clear();
        BaseClass::done_symbolic();
      }

      virtual SolverStatus solve(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // save defect
        this->_vec_def.copy(vec_rhs);
        this->_system_filter.template filter_def<AlgoType>(this->_vec_def);

        // clear solution vector
        vec_sol.format();

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

      virtual SolverStatus correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.template apply<AlgoType>(this->_vec_def, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.template filter_def<AlgoType>(this->_vec_def);

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

    protected:
      virtual SolverStatus _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
      {
        VectorType& vec_def(this->_vec_def);
        VectorType& vec_dir(this->_vec_dir);
        VectorType& vec_tmp(this->_vec_tmp);
        const MatrixType& mat_sys(this->_system_matrix);
        const FilterType& fil_sys(this->_system_filter);
        SolverStatus status(SolverStatus::progress);

        // compute initial defect
        status = this->_set_initial_defect(vec_def);
        if(status != SolverStatus::progress)
          return status;

        // apply preconditioner to defect vector
        if(!this->_apply_precond(vec_dir, vec_def))
          return SolverStatus::aborted;
        fil_sys.template filter_cor<AlgoType>(vec_dir);

        // compute initial gamma
        DataType gamma = vec_def.template dot<AlgoType>(vec_dir);

        // start iterating
        while(status == SolverStatus::progress)
        {
          // compute A*d
          mat_sys.template apply<AlgoType>(vec_tmp, vec_dir);
          fil_sys.template filter_def<AlgoType>(vec_tmp);

          // compute alpha
          DataType alpha = gamma / vec_tmp.template dot<AlgoType>(vec_dir);

          // update solution vector
          vec_sol.template axpy<AlgoType>(vec_dir, vec_sol, alpha);

          // update defect vector
          vec_def.template axpy<AlgoType>(vec_tmp, vec_def, -alpha);

          // compute defect norm
          status = this->_set_new_defect(vec_def);
          if(status != SolverStatus::progress)
            return status;

          // apply preconditioner
          if(!this->_apply_precond(vec_tmp, vec_def))
            return SolverStatus::aborted;
          fil_sys.template filter_cor<AlgoType>(vec_tmp);

          // compute new gamma
          DataType gamma2 = gamma;
          gamma = vec_def.template dot<AlgoType>(vec_tmp);

          // compute beta
          DataType beta = gamma / gamma2;

          // update direction vector
          vec_dir.template axpy<AlgoType>(vec_dir, vec_tmp, beta);
        }

        // we should never reach this point...
        return SolverStatus::undefined;
      }
    }; // class PCGSolver<...>

    /**
     * \brief FGMRES(k) solver implementation
     *
     * This class implements the Restarted Flexible Generalised Minimal Residual solver (FGMRES(k)).
     *
     * \tparam AlgoType_
     * The algorithm tag to be used by the solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Peter Zajac
     */
    template<
      typename AlgoType_,
      typename Matrix_,
      typename Filter_>
    class FGMRESSolver :
      public PreconditionedIterativeSolver<AlgoType_, Matrix_, Filter_>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<AlgoType, MatrixType, FilterType> BaseClass;

      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// krylov dimension
      Index _krylov_dim;
      /// krylov basis vectors
      std::vector<VectorType> _vec_v, _vec_z;
      /// Givens rotation coefficients
      std::vector<DataType> _c, _s, _q;
      /// Hessenberg matrix
      std::vector<std::vector<DataType>> _h;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] krylov_dim
       * The maximum Krylov subspace dimension. Must be > 0.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       *
       * \param[in] del_precond
       * Specifies whether the preconditioner object should be deleted upund
       * destruction of this solver object.
       */
      explicit FGMRESSolver(const MatrixType& matrix, const FilterType& filter, Index krylov_dim,
        PrecondType* precond = nullptr, bool del_precond = false) :
        BaseClass("FGMRES(" + stringify(krylov_dim) + ")", matrix, filter, precond, del_precond),
        _krylov_dim(krylov_dim)
      {
        _c.reserve(krylov_dim);
        _s.reserve(krylov_dim);
        _q.reserve(krylov_dim);
        _h.resize(krylov_dim);
        for(Index i(0); i < krylov_dim; ++i)
          _h.at(i).resize(i+1);
      }

      virtual String name() override
      {
        return "FGMRES";
      }

      virtual bool init_symbolic() override
      {
        if(!BaseClass::init_symbolic())
          return false;
        _vec_v.push_back(this->_system_matrix.create_vector_r());
        for(Index i(0); i < _krylov_dim; ++i)
        {
          _vec_v.push_back(this->_vec_v.front().clone(CloneMode::Layout));
          _vec_z.push_back(this->_vec_v.front().clone(CloneMode::Layout));
        }
        return true;
      }

      virtual void done_symbolic() override
      {
        _vec_v.clear();
        _vec_z.clear();
        BaseClass::done_symbolic();
      }

      virtual SolverStatus solve(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // save input rhs vector as initial defect
        this->_vec_v.at(0).copy(vec_rhs);
        this->_system_filter.template filter_def<AlgoType>(this->_vec_v.at(0));

        // clear solution vector
        vec_sol.format();

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

      virtual SolverStatus correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute initial defect
        this->_system_matrix.template apply<AlgoType>(this->_vec_v.at(0), vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.template filter_def<AlgoType>(this->_vec_v.at(0));

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

    protected:
      virtual SolverStatus _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        // compute initial defect
        SolverStatus status = this->_set_initial_defect(this->_vec_v.at(0));

        // outer GMRES loop
        while(status == SolverStatus::progress)
        {
          _q.clear();
          _s.clear();
          _c.clear();
          _q.push_back(this->_def_cur);

          // normalise v[0]
          this->_vec_v.at(0).template scale<AlgoType>(this->_vec_v.at(0), DataType(1) / _q.back());

          // inner GMRES loop
          Index i(0);
          for(; i < this->_krylov_dim; ++i)
          {
            // apply preconditioner
            if(!this->_apply_precond(this->_vec_z.at(i), this->_vec_v.at(i)))
              return SolverStatus::aborted;
            this->_system_filter.template filter_cor<AlgoType>(this->_vec_z.at(i));

            // v[i+1] := A*z[i]
            this->_system_matrix.template apply<AlgoType>(this->_vec_v.at(i+1), this->_vec_z.at(i));
            this->_system_filter.template filter_def<AlgoType>(this->_vec_v.at(i+1));

            // Gram-Schmidt process
            for(Index k(0); k <= i; ++k)
            {
              this->_h.at(i).at(k) = this->_vec_v.at(i+1).template dot<AlgoType>(this->_vec_v.at(k));
              this->_vec_v.at(i+1).template axpy<AlgoType>(this->_vec_v.at(k), this->_vec_v.at(i+1), -this->_h.at(i).at(k));
            }

            // normalise v[i+1]
            DataType alpha = this->_vec_v.at(i+1).template norm2<AlgoType>();
            this->_vec_v.at(i+1).template scale<AlgoType>(this->_vec_v.at(i+1), DataType(1) / alpha);

            // apply Givens rotations
            for(Index k(0); k < i; ++k)
            {
              DataType t(this->_h.at(i).at(k));
              this->_h.at(i).at(k  ) = this->_c.at(k) * t + this->_s.at(k) * this->_h.at(i).at(k+1);
              this->_h.at(i).at(k+1) = this->_s.at(k) * t - this->_c.at(k) * this->_h.at(i).at(k+1);
            }

            // compute beta
            DataType beta = Math::sqrt(Math::sqr(this->_h.at(i).at(i)) + Math::sqr(alpha));

            // compute next plane rotation
            _s.push_back(alpha / beta);
            _c.push_back(this->_h.at(i).at(i) / beta);

            this->_h.at(i).at(i) = beta;
            this->_q.push_back(this->_s.back() * this->_q.at(i));
            this->_q.at(i) *= this->_c.back();

            // push our new defect
            if((i+1) < this->_krylov_dim)
              this->_set_new_defect_norm(this->_q.back());
          }

          Index n = Math::min(i, this->_krylov_dim);

          // solve H*q = q
          for(Index k(n); k > 0;)
          {
            --k;
            this->_q.at(k) /= this->_h.at(k).at(k);
            for(Index j(k); j > 0;)
            {
              --j;
              this->_q.at(j) -= this->_h.at(k).at(j) * this->_q.at(k);
            }
          }

          // update solution
          for(Index k(0); k < n; ++k)
            vec_sol.template axpy<AlgoType>(this->_vec_z.at(k), vec_sol, this->_q.at(k));

          // compute "real" residual
          this->_system_matrix.template apply<AlgoType>(this->_vec_v.at(0), vec_sol, vec_rhs, -DataType(1));
          this->_system_filter.template filter_def<AlgoType>(this->_vec_v.at(0));

          // set the current defect
          status = this->_set_new_defect(this->_vec_v.at(0));
        }

        // finished
        return status;
      }
    }; // class FGMRESSolver<...>

    /**
     * \brief (Preconditioned) Biconjugate gradient stabilized solver implementation
     *
     * This class implements a simple BiCGStab solver.
     *
     * \tparam AlgoType_
     * The algorithm tag to be used by the solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Christoph Lohmann
     */
    template<
      typename AlgoType_,
      typename Matrix_,
      typename Filter_>
    class BiCGStabSolver :
      public PreconditionedIterativeSolver<AlgoType_, Matrix_, Filter_>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<AlgoType, MatrixType, FilterType> BaseClass;

      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// temporary vectors
      VectorType _vec_r;
      VectorType _vec_r_tilde;
      VectorType _vec_r_tilde_0;
      VectorType _vec_p_tilde;
      VectorType _vec_v;
      VectorType _vec_v_tilde;
      VectorType _vec_s;
      VectorType _vec_s_tilde;
      VectorType _vec_t;
      VectorType _vec_t_tilde;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       *
       * \param[in] del_precond
       * Specifies whether the preconditioner object should be deleted upund
       * destruction of this solver object.
       */
      explicit BiCGStabSolver(const MatrixType& matrix, const FilterType& filter,
        PrecondType* precond = nullptr, bool del_precond = false) :
        BaseClass("BiCGStab", matrix, filter, precond, del_precond)
      {
      }

      virtual String name() override
      {
        return "BiCGStab";
      }

      virtual bool init_symbolic() override
      {
        if(!BaseClass::init_symbolic())
          return false;
        // create all temporary vectors
        // Each *_tilde vector is a correction, all others are defects.
        _vec_r         = this->_system_matrix.create_vector_r();
        _vec_r_tilde   = this->_system_matrix.create_vector_r();
        _vec_r_tilde_0 = this->_system_matrix.create_vector_r();
        _vec_p_tilde   = this->_system_matrix.create_vector_r();
        _vec_v         = this->_system_matrix.create_vector_r();
        _vec_v_tilde   = this->_system_matrix.create_vector_r();
        _vec_s         = this->_system_matrix.create_vector_r();
        _vec_s_tilde   = this->_system_matrix.create_vector_r();
        _vec_t         = this->_system_matrix.create_vector_r();
        _vec_t_tilde   = this->_system_matrix.create_vector_r();
        return true;
      }

      virtual void done_symbolic() override
      {
        this->_vec_r.clear();
        this->_vec_r_tilde.clear();
        this->_vec_r_tilde_0.clear();
        this->_vec_p_tilde.clear();
        this->_vec_v.clear();
        this->_vec_v_tilde.clear();
        this->_vec_s.clear();
        this->_vec_s_tilde.clear();
        this->_vec_t.clear();
        this->_vec_t_tilde.clear();
        BaseClass::done_symbolic();
      }

      virtual SolverStatus correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute initial defect
        this->_system_matrix.template apply<AlgoType_>(this->_vec_r, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.template filter_def<AlgoType_>(this->_vec_r);

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

      virtual SolverStatus solve(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // save rhs vector as initial defect
        this->_vec_r.copy(vec_rhs);
        this->_system_filter.template filter_def<AlgoType_>(this->_vec_r);

        // format solution vector
        vec_sol.format();
        return _apply_intern(vec_sol, vec_rhs);
      }

    protected:
      virtual SolverStatus _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        VectorType& vec_r        (this->_vec_r);
        VectorType& vec_r_tilde  (this->_vec_r_tilde);
        VectorType& vec_r_tilde_0(this->_vec_r_tilde_0);
        VectorType& vec_p_tilde  (this->_vec_p_tilde);
        VectorType& vec_v        (this->_vec_v);
        VectorType& vec_v_tilde  (this->_vec_v_tilde);
        VectorType& vec_s        (this->_vec_s);
        VectorType& vec_s_tilde  (this->_vec_s_tilde);
        VectorType& vec_t        (this->_vec_t);
        VectorType& vec_t_tilde  (this->_vec_t_tilde);
        const MatrixType& mat_sys(this->_system_matrix);
        const FilterType& fil_sys(this->_system_filter);
        SolverStatus status(SolverStatus::progress);

        DataType rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
        //bool early_exit = 0;
        bool restarted = false;

        while(status == SolverStatus::progress)
        {
          if (restarted == false)
          {
            // initial defect is already computed
            status = this->_set_initial_defect(vec_r);
          }
          else
          {
            mat_sys.template apply<AlgoType_>(vec_r, vec_sol, vec_rhs, -DataType(1));
            fil_sys.template filter_def<AlgoType_>(vec_r);
          }

          // apply preconditioner
          if(!this->_apply_precond(vec_r_tilde_0, vec_r))
            return SolverStatus::aborted;
          fil_sys.template filter_cor<AlgoType_>(vec_r_tilde_0);

          vec_r_tilde.copy(vec_r_tilde_0);
          vec_p_tilde.copy(vec_r_tilde_0);

          rho_tilde = vec_r_tilde_0.template dot<AlgoType_>(vec_r_tilde_0);

          // main BiCGStab loop
          while(status == SolverStatus::progress)
          {
            mat_sys.template apply<AlgoType_>(vec_v, vec_p_tilde);
            fil_sys.template filter_def<AlgoType_>(vec_v);
            // apply preconditioner
            if(!this->_apply_precond(vec_v_tilde, vec_v))
              return SolverStatus::aborted;
            fil_sys.template filter_cor<AlgoType_>(vec_v_tilde);

            gamma_tilde = vec_v_tilde.template dot<AlgoType_>(vec_r_tilde_0);

            if (Math::abs(gamma_tilde) < Math::abs(rho_tilde)*1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 1" << std::endl;
              break;
            }

            alpha_tilde = rho_tilde / gamma_tilde;

            if ((Math::abs(alpha_tilde) * vec_v_tilde.template norm2<AlgoType_>()) / this->_def_cur < 1e-5)
            {
              restarted = true;
              //std::cout << "Breakpoint 2" << std::endl;
              // \TODO warum ist das break hier nicht aktiv?
              //break;
            }

            DataType malpha_tilde(-alpha_tilde);
            vec_s.template axpy<AlgoType_>(vec_v, vec_r, malpha_tilde);

            // compute new defect norm
            /*status = this->_set_new_defect(vec_s);
            if (status == SolverStatus::success)
            {
              vec_sol.template axpy<AlgoType_>(vec_p_tilde, vec_sol, alpha_tilde);

              //early_exit = 1;
              //std::cout << "Breakpoint 3 (converged)" << std::endl;
              return status;
            }*/
            vec_s_tilde.template axpy<AlgoType_>(vec_v_tilde, vec_r_tilde, malpha_tilde);

            mat_sys.template apply<AlgoType_>(vec_t, vec_s_tilde);
            fil_sys.template filter_def<AlgoType_>(vec_t);

            // apply preconditioner
            if(!this->_apply_precond(vec_t_tilde, vec_t))
              return SolverStatus::aborted;
            fil_sys.template filter_cor<AlgoType_>(vec_t_tilde);

            gamma_tilde = vec_t_tilde.template dot<AlgoType_>(vec_t_tilde);
            omega_tilde = vec_t_tilde.template dot<AlgoType_>(vec_s_tilde);

            if (Math::abs(gamma_tilde) < Math::abs(omega_tilde) * 1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 4" << std::endl;
              break;
            }
            omega_tilde = omega_tilde / gamma_tilde;

            vec_sol.template axpy<AlgoType_>(vec_s_tilde, vec_sol, omega_tilde);
            vec_sol.template axpy<AlgoType_>(vec_p_tilde, vec_sol, alpha_tilde);

            DataType momega_tilde(-omega_tilde);
            vec_r.template axpy<AlgoType_>(vec_t, vec_s, momega_tilde);

            // compute new defect norm
            status = this->_set_new_defect(vec_r);
            if (status == SolverStatus::success)
            {
              //std::cout << "Breakpoint 5 (converged)" << std::endl;
              return status;
            }

            vec_r_tilde.template axpy<AlgoType_>(vec_t_tilde, vec_s_tilde, momega_tilde);

            rho_tilde_old = rho_tilde;
            rho_tilde = vec_r_tilde.template dot<AlgoType_>(vec_r_tilde_0);

            beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

            vec_p_tilde.template axpy<AlgoType_>(vec_v_tilde, vec_p_tilde, momega_tilde);
            vec_p_tilde.template scale<AlgoType_>(vec_p_tilde, beta_tilde);
            vec_p_tilde.template axpy<AlgoType_>(vec_p_tilde, vec_r_tilde);
          }
        }

        // finished
        return status;
      }
    }; // class BiCGStabSolver<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_PROTO_SOLVER_HPP
