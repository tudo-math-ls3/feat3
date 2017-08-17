#pragma once
#ifndef KERNEL_SOLVER_BASE_HPP
#define KERNEL_SOLVER_BASE_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <memory>

namespace FEAT
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
     * If the solver was created by calling the PropertyMap based constructor, it can return the name of the
     * corresponding section by \c get_section_name(). Otherwise, it returns its name through the virtual \c name()
     * function.
     *
     * \author Peter Zajac
     */
    template<typename Vector_>
    class SolverBase
    {
    public:
      /// The type of vector this solver can be applied to
      typedef Vector_ VectorType;

    private:
      /// The name of the PropertyMap section this solver was configured by (if any)
      String _section_name;

    public:
      /**
       * \brief Empty standard constructor
       */
      SolverBase()
        : _section_name("")
      {
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * This configures the solver using the key/value pairs in the PropertyMap. This base class has nothing to
       * configure, but this constructor is called from the derived classes as well.
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       */
      explicit SolverBase(const String& config_section_name, PropertyMap* config_section) :
        _section_name(config_section_name)
      {
        XASSERT(config_section != nullptr);
      }

      /**
       * \brief Virtual destructor
       */
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
         \endverbatim
       */
      virtual void init()
      {
        init_symbolic();
        init_numeric();
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

      /**
       * \brief Returns the name of the PropertyMap section this solver was configured by
       *
       * If the standard constructor was used, _section_name is the empty String and name() is returned instead.
       * This is important for derived classes' write_config() functionality, so that the written section can be
       * assigned a name.
       */
      virtual String get_section_name() const
      {
        if(_section_name == "")
        {
          return name();
        }
        else
        {
          return _section_name;
        }
      }

      /**
       * \brief Writes the solver configuration to a PropertyMap
       *
       * \param[in] parent
       * The PropertyMap to add a new section to
       *
       * \param[in] new_section_name
       * The name the new section will have
       *
       * \returns A pointer to the new section created in the PropertyMap
       */
      virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name = "") const
      {
        XASSERT(parent != nullptr);

        PropertyMap* config_section(nullptr);

        if(new_section_name == "")
        {
          if(_section_name == "")
          {
            config_section = parent->add_section(name());
          }
          else
          {
            config_section = parent->add_section(_section_name);
          }
        }
        else
        {
          config_section = parent->add_section(new_section_name);
        }

        config_section->add_entry("type",name());

        return config_section;
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

    /**
     * \brief Helper class for iteration statistics collection
     *
     * This helper class collects timings for iterative solvers.
     * Create an instance (object) of this class at the beginning of
     * each iteration.
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    class IterationStats
    {
    private:
      const String _solver_name;
      TimeStamp _at;
      double _mpi_execute_reduction_start;
      double _mpi_execute_reduction_stop;
      double _mpi_execute_spmv_start;
      double _mpi_execute_spmv_stop;
      double _mpi_execute_collective_start;
      double _mpi_execute_collective_stop;
      double _mpi_wait_start_reduction;
      double _mpi_wait_start_spmv;
      double _mpi_wait_start_collective;
      double _mpi_wait_stop_reduction;
      double _mpi_wait_stop_spmv;
      double _mpi_wait_stop_collective;
      bool _destroyed;

    public:
      /**
       * \brief Constructor
       *
       * This constructor is called at the beginning of each iteration of
       * an iterative solver.
       *
       * \param[in] solver
       * The iterative solver object that is to benchmarked.
       */
      template<typename Vector_>
      explicit IterationStats(const SolverBase<Vector_>& solver) :
        _solver_name(solver.name()),
        _destroyed(false)
      {
        _mpi_execute_reduction_start = Statistics::get_time_mpi_execute_reduction();
        _mpi_execute_spmv_start = Statistics::get_time_mpi_execute_spmv();
        _mpi_execute_collective_start = Statistics::get_time_mpi_execute_collective();
        _mpi_wait_start_reduction    = Statistics::get_time_mpi_wait_reduction();
        _mpi_wait_start_spmv    = Statistics::get_time_mpi_wait_spmv();
        _mpi_wait_start_collective    = Statistics::get_time_mpi_wait_collective();
        _mpi_execute_reduction_stop = _mpi_execute_reduction_start;
        _mpi_execute_spmv_stop = _mpi_execute_spmv_start;
        _mpi_execute_collective_stop = _mpi_execute_collective_start;
        _mpi_wait_stop_reduction    = _mpi_wait_start_reduction;
        _mpi_wait_stop_spmv    = _mpi_wait_start_spmv;
        _mpi_wait_stop_collective    = _mpi_wait_start_collective;
      }

      // delete copy-ctor and assign operator
      IterationStats(const IterationStats&) = delete;
      IterationStats& operator=(const IterationStats&) = delete;

      /**
       * \brief Destructor
       *
       * This is called at the end of each iteration. This function computes the
       * iteration runtimes and commits them to the Statistics collection system.
       */
      ~IterationStats()
      {
        if (!_destroyed)
          destroy();
      }

      ///destroy the objects contents (and generate Statistics::expression) before the actual destructor call
      void destroy()
      {
        if (_destroyed)
          throw InternalError(__func__, __FILE__, __LINE__, "IterationStats destroy method was already called before!");

        _mpi_execute_reduction_stop = Statistics::get_time_mpi_execute_reduction();
        _mpi_execute_spmv_stop = Statistics::get_time_mpi_execute_spmv();
        _mpi_execute_collective_stop = Statistics::get_time_mpi_execute_collective();
        _mpi_wait_stop_reduction    = Statistics::get_time_mpi_wait_reduction();
        _mpi_wait_stop_spmv    = Statistics::get_time_mpi_wait_spmv();
        _mpi_wait_stop_collective    = Statistics::get_time_mpi_wait_collective();
        Statistics::add_solver_expression(std::make_shared<ExpressionTimings>(_solver_name, _at.elapsed_now(),
              _mpi_execute_reduction_stop - _mpi_execute_reduction_start, _mpi_execute_spmv_stop - _mpi_execute_spmv_start, _mpi_execute_collective_stop - _mpi_execute_collective_start,
          _mpi_wait_stop_reduction - _mpi_wait_start_reduction, _mpi_wait_stop_spmv - _mpi_wait_start_spmv, _mpi_wait_stop_collective - _mpi_wait_start_collective));

        _destroyed = true;
      }
    }; // class IterationStats
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_BASE_HPP
