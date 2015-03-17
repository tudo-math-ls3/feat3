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
     * The class of the vector that is passed to the solver in the \c apply method.
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
       * \brief Solver application method
       *
       * This method applies the solver represented by this object onto a given right-hand-side vector
       * and returns the corresponding solution vector.
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
      virtual SolverStatus apply(IOVector_& vec_sol, const IOVector_& vec_rhs) = 0;
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
       * \param[in] matrix
       * A reference to the matrix for the preconditioner.
       *
       * \param[in] args
       * The remaining arguments which are passed to the preconditioner. For the required set of arguments,
       * see the documentation of the corresponding preconditioner class template.
       */
      template<typename... Args_>
      explicit PreconWrapper(const Matrix_& matrix, Args_&&... args) :
        _precond(matrix, std::forward<Args_>(args)...)
      {
      }

      /// Applies the preconditioner.
      virtual SolverStatus apply(VectorType& vec_sol, const VectorType& vec_rhs) override
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

      virtual SolverStatus apply(VectorType& vec_sol, const VectorType& vec_rhs) override
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
     * \tparam Vector_
     * The vector class used by the solver.
     *
     * \author Peter Zajac
     */
    template<typename AlgoType_, typename Vector_>
    class IterativeSolver :
      public SolverInterface<Vector_>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Vector_ VectorType;
      typedef typename VectorType::DataType DataType;
      typedef SolverInterface<VectorType> BaseClass;

    protected:
      /// name of the solver
      String _name;
      /// relative tolerance parameter
      DataType _tol_rel;
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
       * \param[in] name
       * Specifies the name of the iterative solver. This is used as a prefix for the convergence plot.
       *
       * This constructor initialises the following values:
       *
       * - relative tolerance: 1E-8
       * - minimum iterations: 0
       * - maximum iterations: 100
       * - convergence plot: false
       */
      explicit IterativeSolver(String name) :
        BaseClass(),
        _name(name),
        _tol_rel(1E-8),
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
        return ((_def_cur / _def_init) <= _tol_rel);
      }

      /// checks for divergence
      bool is_diverged() const
      {
        return _def_cur > DataType(1E+99);
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
          std::cout << this->_name << ": " << stringify(0).pad_front(_iter_digits)
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
          std::cout << this->_name << ": " << stringify(this->_num_iter).pad_front(_iter_digits)
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
          std::cout << this->_name << "* " << stringify(this->_num_iter).pad_front(_iter_digits)
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
     * \tparam Vector_
     * The vector class used by the solver.
     *
     * \authro Peter Zajac
     */
    template<typename AlgoType_, typename Vector_>
    class PreconditionedIterativeSolver :
      public IterativeSolver<AlgoType_, Vector_>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Vector_ VectorType;
      typedef typename VectorType::DataType DataType;
      typedef IterativeSolver<AlgoType_, VectorType> BaseClass;

      /// the interface for the preconditioner
      typedef SolverInterface<VectorType> PrecondType;

    protected:
      // a pointer to the preconditioner
      PrecondType* _precond;

      /**
       * \brief Constructor
       *
       * \param[in] name
       * The name of the solver.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be nullptr.
       *
       * \attention
       * This class does \b not delete the preconditioner object upon destruction!
       */
      explicit PreconditionedIterativeSolver(String name, PrecondType* precond = nullptr) :
        BaseClass(name),
        _precond(precond)
      {
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
          return status_success(this->_precond->apply(vec_out, vec_in));
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
      public PreconditionedIterativeSolver<AlgoType_, typename Matrix_::VectorTypeR>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<AlgoType, VectorType> BaseClass;

      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// the system matrix
      const MatrixType& _system_matrix;
      /// the system filter
      const FilterType& _system_filter;
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
       */
      explicit FixPointSolver(const MatrixType& matrix, const FilterType& filter, PrecondType* precond = nullptr) :
        BaseClass("FixPoint", precond),
        _system_matrix(matrix),
        _system_filter(filter)
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
      }

      virtual SolverStatus apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        VectorType& vec_def(this->_vec_def);
        VectorType& vec_cor(this->_vec_cor);
        const MatrixType& mat_sys(this->_system_matrix);
        const FilterType& fil_sys(this->_system_filter);
        SolverStatus status(SolverStatus::progress);

        // save defect
        vec_def.copy(vec_rhs);
        fil_sys.template filter_def<AlgoType>(vec_def);

        // clear solution vector
        vec_sol.format();

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

        // we should never reach this point...
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
      public PreconditionedIterativeSolver<AlgoType_, typename Matrix_::VectorTypeR>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<AlgoType, VectorType> BaseClass;

      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// the system matrix
      const MatrixType& _system_matrix;
      /// the system filter
      const FilterType& _system_filter;
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
       */
      explicit PCGSolver(const MatrixType& matrix, const FilterType& filter, PrecondType* precond = nullptr) :
        BaseClass("PCG", precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
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
      }

      virtual SolverStatus apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        VectorType& vec_def(this->_vec_def);
        VectorType& vec_dir(this->_vec_dir);
        VectorType& vec_tmp(this->_vec_tmp);
        const MatrixType& mat_sys(this->_system_matrix);
        const FilterType& fil_sys(this->_system_filter);
        SolverStatus status(SolverStatus::progress);

        // save defect
        vec_def.copy(vec_rhs);
        fil_sys.template filter_def<AlgoType>(vec_def);

        // clear solution vector
        vec_sol.format();

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
      public PreconditionedIterativeSolver<AlgoType_, typename Matrix_::VectorTypeR>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<AlgoType, VectorType> BaseClass;

      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// the system matrix
      const MatrixType& _system_matrix;
      /// the system filter
      const FilterType& _system_filter;
      /// krylov dimension
      Index _krylov_dim;
      /// right-hand-side vector
      VectorType _vec_b;
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
       */
      explicit FGMRESSolver(const MatrixType& matrix, const FilterType& filter, Index krylov_dim, PrecondType* precond = nullptr) :
        BaseClass("FGMRES(" + stringify(krylov_dim) + ")", precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _krylov_dim(krylov_dim)
      {
        _c.reserve(krylov_dim);
        _s.reserve(krylov_dim);
        _q.reserve(krylov_dim);
        _h.resize(krylov_dim);
        for(Index i(0); i < krylov_dim; ++i)
          _h.at(i).resize(i+1);
      }

      virtual bool init_symbolic() override
      {
        if(!BaseClass::init_symbolic())
          return false;
        _vec_b = this->_system_matrix.create_vector_r();
        _vec_v.push_back(this->_vec_b.clone(CloneMode::Layout));
        for(Index i(0); i < _krylov_dim; ++i)
        {
          _vec_v.push_back(this->_vec_b.clone(CloneMode::Layout));
          _vec_z.push_back(this->_vec_b.clone(CloneMode::Layout));
        }
        return true;
      }

      virtual void done_symbolic() override
      {
        _vec_v.clear();
        _vec_z.clear();
      }

      virtual SolverStatus apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // save input defect vector
        this->_vec_b.copy(vec_rhs);
        this->_system_filter.template filter_def<AlgoType>(this->_vec_b);
        this->_vec_v.at(0).copy(this->_vec_b);

        // clear solution vector
        vec_sol.format();

        // compute initial defect
        SolverStatus status = this->_set_initial_defect(this->_vec_b);

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
          this->_system_matrix.template apply<AlgoType>(this->_vec_v.at(0), vec_sol, this->_vec_b, -DataType(1));
          this->_system_filter.template filter_def<AlgoType>(this->_vec_v.at(0));

          // set the current defect
          status = this->_set_new_defect(this->_vec_v.at(0));
        }

        // finished
        return status;
      }
    }; // class FGMRESSolver<...>

    /// a simple SSOR-preconditioner
    template<typename Matrix_>
    class SSORPrecond;

    /// SSOR-specialisation for CSR matrices stored in Main memory
    template<
      typename DataType_,
      typename IndexType_>
    class SSORPrecond<LAFEM::SparseMatrixCSR<Mem::Main, DataType_, IndexType_>>
      : public SolverInterface<LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>>
    {
    public:
      typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_, IndexType_> MatrixType;
      typedef LAFEM::DenseVector<Mem::Main, DataType_, IndexType_> VectorType;
      typedef SolverInterface<LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>> BaseClass;
      typedef Algo::Generic AlgoType;

    protected:
      const MatrixType& _matrix;
      DataType_ _relax;

    public:
      explicit SSORPrecond(const MatrixType& matrix, DataType_ relax = DataType_(1)) :
        _matrix(matrix),
        _relax(relax)
        {
        }

      void set_relax(DataType_ relax)
      {
        _relax = relax;
      }

      virtual SolverStatus apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        const IndexType_* row_ptr(this->_matrix.row_ptr());
        const IndexType_* col_idx(this->_matrix.col_ind());
        const DataType_* data(this->_matrix.val());
        DataType_* x(vec_sol.elements());
        const DataType_* y(vec_rhs.elements());

        const Index n = this->_matrix.rows();

        // forward loop
        for(Index i(0); i < n; ++i)
        {
          DataType_ d(0);
          Index j(0);
          for(j = row_ptr[i]; col_idx[j] < i; ++j)
            d += data[j] * x[col_idx[j]];
          x[i] = (y[i] - _relax * d) / data[j];
        }

        // backward loop
        for(Index i(n); i > 0; )
        {
          --i;
          DataType_ d(0);
          Index j(0);
          for(j = row_ptr[i+1]-1; col_idx[j] > i; --j)
            d += data[j] * x[col_idx[j]];
          x[i] -= _relax * d / data[j];
        }

        // okay
        return SolverStatus::success;
      }
    };

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
      public PreconditionedIterativeSolver<AlgoType_, typename Matrix_::VectorTypeR>
    {
    public:
      typedef AlgoType_ AlgoType;
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<AlgoType, VectorType> BaseClass;

      typedef SolverInterface<VectorType> PrecondType;

    protected:
      /// the system matrix
      const MatrixType& _system_matrix;
      /// the system filter
      const FilterType& _system_filter;

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
       */
      explicit BiCGStabSolver(const MatrixType& matrix, const FilterType& filter, PrecondType* precond = nullptr) :
        BaseClass("BiCGStab", precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
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
      }

      virtual SolverStatus apply(VectorType& vec_sol, const VectorType& vec_rhs) override
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
          mat_sys.template apply<AlgoType_>(vec_r, vec_sol, vec_rhs, -DataType(1));
          fil_sys.template filter_def<AlgoType_>(vec_r);

          if (restarted == false)
          {
            status = this->_set_initial_defect(vec_r);
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
