#pragma once
#ifndef KERNEL_SOLVER_BASE_HPP
#define KERNEL_SOLVER_BASE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <memory>

namespace FEAST
{
  /**
   * \brief Solver namespace
   */
  namespace Solver
  {
    /**
     * \brief Solver status return codes enumeration
     *
     * This enumeration defined the solver status return codes, which specify whether
     * the solver application was successful or failed due to some reason.
     */
    enum class Status
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

    /// \cond internal
    inline std::ostream& operator<<(std::ostream& os, Status status)
    {
      switch(status)
      {
      case Status::undefined:
        return os << "undefined";
      case Status::progress:
        return os << "progress";
      case Status::success:
        return os << "success";
      case Status::aborted:
        return os << "aborted";
      case Status::diverged:
        return os << "diverged";
      case Status::max_iter:
        return os << "max-iter";
      case Status::stagnated:
        return os << "stagnated";
      default:
        return os << "-unknown-";
      }
    }
    /// \endcond

    /**
     * \brief Status success check function
     *
     * This function takes a Status value as input and checks whether it respresents a
     * 'successul' run. A solving run is interpreted as successful, if one of the following
     * status codes was returned:
     *
     *  - Status::success
     *  - Status::max_iter
     *  - Status::stagnated
     *
     * For any other status code, the solving run is interpreted as unsuccessful.
     *
     * \param[in] status
     * A status code returned by a solver.
     *
     * \returns
     * \c true, if the run was successful, otherwise \c fase.
     */
    inline bool status_success(Status status)
    {
      switch(status)
      {
      case Status::success:
      case Status::max_iter:
      case Status::stagnated:
        return true;

      default:
        return false;
      }
    }

    /**
     * \brief Base-class for solver generated exceptions.
     */
    class SolverException :
      public Exception
    {
    public:
      explicit SolverException(const String& msg) :
        Exception(msg)
      {
      }
    };

    /**
     * \brief Invalid Matrix structure exception
     *
     * This exception may be thrown by a solver object if the matrix structure
     * is not valid.
     */
    class InvalidMatrixStructureException :
      public SolverException
    {
    public:
      explicit InvalidMatrixStructureException(const String& msg = "invalid matrix structure") :
        SolverException(msg)
      {
      }
    };

    /**
     * \brief Singular Matrix exception
     *
     * This exception may be thrown by a solver object if it detects that
     * the matrix to be factorised is singular.
     */
    class SingularMatrixException :
      public SolverException
    {
    public:
      explicit SingularMatrixException(const String& msg = "singular matrix") :
        SolverException(msg)
      {
      }
    };

    /**
     * \brief Polymorphic solver interface
     *
     * This class template defines the interface for any solver implementation.
     *
     * \tparam Vector_
     * The class of the vector that is passed to the solver in the \c solve() method.
     *
     * \author Peter Zajac
     */
    template<typename Vector_>
    class SolverBase
    {
    private:

      /// Description of this solvers branch, with itself at the leaf.
      String _branch;

    public:

      /// virtual destructor
      virtual ~SolverBase()
      {
      }

      /**
       * \brief Symbolic initialisation method
       *
       * This method is called to perform symbolic initialisation of the solver.
       */
      virtual void init_symbolic()
      {
      }

      /**
       * \brief Numeric initialisation method
       *
       * This method is called to perform numeric initialisation of the solver.\n
       * Before this function can be called, the symbolic initialisation must be performed.
       */
      virtual void init_numeric()
      {
      }

      /// Initialise solver branch description
      virtual void init_branch(String root = "")
      {
        _branch = root + "::" + name();
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
       * This function performs both the symbolic and numeric initialisation, i.e. it simply performs
         \verbatim
         this->init_symbolic();
         this->init_numeric();
         this->init_branch();
         \endverbatim
       */
      virtual void init()
      {
        init_symbolic();
        init_numeric();
        init_branch();
      }

      /**
       * \brief Finalisation method
       *
       * This function performs both the symbolic and numeric finalisation, i.e. it simply performs
         \verbatim
         this->done_numeric();
         this->done_symbolic();
         \endverbatim
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
      virtual String name() const = 0;

      String get_solver_branch()
      {
        return _branch;
      }

      /**
       * \brief Solver application method
       *
       * This method applies the solver represented by this object onto a given defect vector
       * and returns the corresponding correction vector.
       *
       * \note Solvers which derive from the IterativeSolver base class also provide a \c correct() method
       * which corrects an initial solution instead of starting with the null vector.
       *
       * \param[out] vec_cor
       * The vector that shall receive the solution of the linear system. It is assumed to be allocated, but
       * its numerical contents may be undefined upon calling this method.
       *
       * \param[in] vec_def
       * The vector that represents the right-hand-side of the linear system to be solved.
       *
       * \returns
       * A Status code that represents the status of the solution step.
       */
      virtual Status apply(Vector_& vec_cor, const Vector_& vec_def) = 0;
    }; // class SolverBase<...>

    /**
     * \brief Solve linear system with initial solution guess.
     *
     * \param[in] solver
     * The solver object that is to be used for solving the system.
     *
     * \param[in,out] vec_sol
     * The solution vector containing an initial guess.
     *
     * \param[in] vec_rhs
     * The right-hand-side vector of the linear system.
     *
     * \param[in] matrix
     * The matrix of the linear system.
     *
     * \param[in] filter
     * The filter of the linear system.
     *
     * \returns
     * A status code representing the status of the solving step.
     */
    template<
      typename Vector_,
      typename Matrix_,
      typename Filter_>
    inline Status solve(
      SolverBase<Vector_>& solver,
      Vector_& vec_sol,
      const Vector_& vec_rhs,
      const Matrix_& matrix,
      const Filter_& filter)
    {
      typedef typename Vector_::DataType DataType;

      // create two temporary vectors
      Vector_ vec_def(vec_rhs.clone(LAFEM::CloneMode::Layout));
      Vector_ vec_cor(vec_sol.clone(LAFEM::CloneMode::Layout));

      // compute defect
      matrix.apply(vec_def, vec_sol, vec_rhs, -DataType(1));

      // apply defect filter
      filter.filter_def(vec_def);

      // apply solver
      Status status = solver.apply(vec_cor, vec_def);
      if(status_success(status))
      {
        // apply correction filter
        filter.filter_cor(vec_cor);

        // update solution
        vec_sol.axpy(vec_sol, vec_cor);
      }

      // return status
      return status;
    }
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_BASE_HPP
