// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/statistics.hpp>

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    namespace Intern
    {
      /// returns the underlying comm, if Matrix_ is a Global::Matrix - SFINAE at its best
      template<typename Matrix_>
      const Dist::Comm* get_mat_comm(const Matrix_& mat, typename Matrix_::GateRowType*)
      {
        return mat.get_comm();
      }

      template<typename Matrix_>
      const Dist::Comm* get_mat_comm(const Matrix_&, ...)
      {
        return nullptr;
      }

      /// returns the underlying comm, if Vector_ is a Global::Vector - SFINAE at its best
      template<typename Vector_>
      const Dist::Comm* get_vec_comm(const Vector_& vec, typename Vector_::GateType*)
      {
        return vec.get_comm();
      }

      template<typename Vector_>
      const Dist::Comm* get_vec_comm(const Vector_&, ...)
      {
        return nullptr;
      }
    } // namespace Intern
    /// \endcond

    /**
     * \brief Solver plot modes enumeration
     *
     */
    enum class PlotMode
    {
      /// No plotting whatsoever
      none = 0,
      /// Plot every iteration (if applicable)
      iter,
      /// Plot a summary after each solver run
      summary,
      /// Plot every iteration (if applicable) and a summary
      all
    };

    /// \cond internal
    inline std::ostream& operator<<(std::ostream& os, PlotMode mode)
    {
      switch(mode)
      {
        case PlotMode::none:
          return os << "none";
        case PlotMode::iter:
          return os << "iter";
        case PlotMode::summary:
          return os << "summary";
        case PlotMode::all:
          return os << "all";
        default:
          return os << "-unknown-";
      }
    }

    inline std::istream& operator>>(std::istream& is, PlotMode& mode)
    {
      String s;
      if((is >> s).fail())
        return is;

      if(s.compare_no_case("none") == 0)
        mode = PlotMode::none;
      else if(s.compare_no_case("iter") == 0)
        mode = PlotMode::iter;
      else if(s.compare_no_case("summary") == 0)
        mode = PlotMode::summary;
      else if(s.compare_no_case("all") == 0)
        mode = PlotMode::all;
      else
        is.setstate(std::ios_base::failbit);

      return is;
    }
    /// \endcond

    /**
     * \brief Abstract base-class for iterative solvers.
     *
     * This class template acts as an abstract base class for iterative solvers and it implements various auxiliary
     * features for convergence control, output formatting and a posteriori convergence analysis.
     *
     * \tparam Vector_
     * The class of the vector that is passed to the solver in the \c solve() method.
     *
     * <u><b>Details about the convergence control mechanisms offered by this class:</b></u>\n
     * This class offers a wide variety of convergence control mechanisms and stopping criterions, which are summarized
     * in this section.
     *
     * In the following, let \f$k\geq 0\f$ denote the current iteration and let \f$r_0 := b - Ax_0\f$ denote the initial
     * defect/residual vector and let \f$r_k := b - Ax_k\f$ denote the defect/residual vector of the current iteration.
     *
     * <b>1. Divergence Criterion:</b>
     * The first stopping criterion checks whether the residual vector norm has risen above a relative or absolute
     * divergence tolerance and, if so, the solver stops and returns the status code Status::diverged.\n
     * The divergence criterion is defined as follows:
     * \f[ (\vert \vert r_k\vert\vert \geq \text{div\_abs}) \lor (\vert \vert r_k\vert\vert \geq \text{div\_rel}\cdot\vert \vert r_0\vert\vert) \f]
     * By default, the divergence tolerances are chosen as
     * - div_abs = 1/eps^2
     * - dil_rel = 1/eps
     *
     * where eps is the machine precision constant for the solvers underlying data type, which evaluates to approximately
     * div_abs = 1E+32 and div_rel = 1E+16 for double precision.
     *
     * <b>2. Maximum Iteration Number Criterion:</b>
     * The second stopping criterion checks whether the solver has reached the maximum number of allowed
     * iterations, i.e. whether k >= max_iter, and, if so, the solver returns the status code Status::max_iter.
     *
     * <b>3. Convergence Criterion:</b>
     * The third stopping criterion checks whether the residual vector has fallen below a relative and/or absolute
     * convergence tolerance and, if so, the solver returns the status code Status::converged.
     * The convergence criterion consists of a relative tolerance and two absolute tolerances and is defined as follows:
     *
     * \f[ (\vert \vert r_k\vert\vert \leq \text{tol\_abs}) \land ((\vert \vert r_k\vert\vert \leq \text{tol\_rel}\cdot\vert \vert r_0\vert\vert) \lor (\vert \vert r_k\vert\vert \leq \text{tol\_abs\_low})) \f]
     * By default, the convergence tolerances are chosen as
     * - tol_rel = sqrt(eps)
     * - tol_abs = 1/eps^2
     * - tol_abs_low = 0
     *
     * where eps is the machine precision constant for the solvers underlying data type, which evaluates to approximately
     * tol_abs = 1E-8 and tol_rel = 1E+32 for double precision, which effectively "disables" both absolute tolerances
     * therefore resulting in a convergence criterion that only depends on the relative tolerance.
     *
     * The tol_abs_low parameter can be used to "disarm" an overly strict relative tolerance tol_rel in the case where
     * the initial defect is already quite small and therefore fulfilling the relative tolerance may become impossible
     * due to finite precision arithmetic combined with a large condition number. A typical use case would be to set
     * tol_abs to 1E-8 and tol_rel to 1E-5 (thus asking that the solver should gain 5 digits) and tol_abs_low to 1E-12,
     * which will "disable" the tol_rel criterion if (and only if) the initial defect norm is less than 1E-7. If the
     * initial defect norm is larger than 1E-7, then the tol_abs_low parameter will become irrelevant and the tol_abs
     * and tol_rel criterions have to be fulfilled for convergence.
     *
     * - If you want the solver to fulfill <b>only a relative</b> tolerance, then set tol_rel to the corresponding
     *   tolerance and leave tol_abs_low set to 0 and tol_abs to set a large value (1/eps^2 by default).
     *
     * - If you want the solver to fulfill <b>only an absolute</b> tolerance, then set both tol_abs and tol_abs_low
     *   to the corresponding tolerance and set tol_rel to any value you like, e.g. 0.
     *
     * - If you want the solver to fulfill <b>either a relative or an absolute</b> tolerance, then set tol_rel and
     *   tol_abs_low to the corresponding tolerances and leave tol_abs to set a large value (1/eps^2 by default).
     *
     * - If you want the solver to fulfill <b>both a relative and an absolute</b> tolerance, then set tol_rel and
     *   tol_abs to the corresponding tolerances and leave tol_abs_low set to 0.
     *
     * <b>5. Stagnation Criterion:</b>
     * The fifth stopping criterion checks whether the residual vector norms do not decrease by a minimum factor over
     * the last few iterations, thus indicating that the solver seems to be stagnating and that further iterations
     * will probably not improve the solution any further and, if so, the solver stops and returns the status code
     * Status::stagnated.\n
     * The stagnation criterion is fulfilled if min_stag_iter is set to a value > 0 and if
     * \f[ \vert \vert r_k\vert\vert \geq \text{stag\_rate}\cdot\vert \vert r_{k-1}\vert\vert \f]
     * for at least min_stag_iter consecutive iterations.
     *
     * By default, the stagnation criterion is disabled and its parameters are chosen as
     * - stag_rate = 0.95
     * - min_stag_iter = 0
     *
     * It is important to keep in mind that the stagnation criterion can usually only be helpful for simple fix-point
     * type iterative solvers, which have more or less monotonous convergence behavior. Other types of solvers,
     * especially Krylov subspace solvers, often show a somewhat "ballistic" convergence curve, where the defect norms
     * seem to indicate divergence for the first few dozen/hundred iterations before the defect norms start falling
     * down very quickly -- in this case, the stagnation criterion will trigger a false positive during the initial
     * climb of the defect norms. Also, some solvers can run into a stagnation scenario where the defect norms
     * oscillate, which will not be detected by the stagnation criterion because the defect norm drops every other
     * iteration before rising back up again in the next iteration.
     *
     * \author Peter Zajac
     */
    template<typename Vector_>
    class IterativeSolver :
      public SolverBase<Vector_>
    {
    public:
      /// The vector type this solver can be applied to
      typedef Vector_ VectorType;
      /// Floating point type
      typedef typename VectorType::DataType DataType;
      /// The base class
      typedef SolverBase<VectorType> BaseClass;

    protected:
      /// Communicator of the solver
      const Dist::Comm* _comm;
      /// name of the solver in plots
      String _plot_name;
      /// current status of the solver
      Status _status;
      /// relative tolerance parameter
      DataType _tol_rel;
      /// absolute tolerance parameter
      DataType _tol_abs;
      /// absolute low tolerance parameter
      DataType _tol_abs_low;
      /// relative divergence parameter
      DataType _div_rel;
      /// absolute divergence parameter
      DataType _div_abs;
      /// stagnation rate
      DataType _stag_rate;
      /// minimum number of iterations
      Index _min_iter;
      /// maximum number of iterations
      Index _max_iter;
      /// number of performed iterations
      Index _num_iter;
      /// minimum number of stagnation iterations
      Index _min_stag_iter;
      /// number of consecutive stagnated iterations
      Index _num_stag_iter;
      /// initial defect
      DataType _def_init;
      /// current defect
      DataType _def_cur;
      /// previous iteration defect
      DataType _def_prev;
      /// iteration count digits for plotting
      Index _iter_digits;
      /// whether to plot something
      PlotMode _plot_mode;
      /// plot output interval
      Index _plot_interval;
      /// whether skipping of defect computation is allowed or not
      bool _force_def_norm_calc;

      /**
       * \brief Protected constructor
       *
       * This constructor initializes the following values:
       *
       * - relative tolerance: sqrt(eps) (~1E-8 for double)
       * - absolute tolerance: 1/eps^2 (~1E+32 for double)
       * - lower absolute tolerance: 0
       * - relative divergence: 1/eps (~1E+16 for double)
       * - absolute divergence: 1/eps^2 (~1E+32 for double)
       * - stagnation rate: 0.95
       * - minimum iterations: 0
       * - maximum iterations: 100
       * - minimum stagnation iterations: 0
       * - convergence plot: false
       *
       * \param[in] plot_name
       * Specifies the name of the iterative solver. This is used as a prefix for the convergence plot.
       */
      explicit IterativeSolver(const String& plot_name) :
        BaseClass(),
        _comm(nullptr),
        _plot_name(plot_name),
        _status(Status::undefined),
        _tol_rel(Math::sqrt(Math::eps<DataType>())),
        _tol_abs(DataType(1) / Math::sqr(Math::eps<DataType>())),
        _tol_abs_low(DataType(0)),
        _div_rel(DataType(1) / Math::eps<DataType>()),
        _div_abs(DataType(1) / Math::sqr(Math::eps<DataType>())),
        _stag_rate(DataType(0.95)),
        _min_iter(0),
        _max_iter(100),
        _num_iter(0),
        _min_stag_iter(0),
        _num_stag_iter(0),
        _def_init(0),
        _def_cur(0),
        _def_prev(0),
        _iter_digits(Math::ilog10(_max_iter)),
        _plot_mode(PlotMode::none),
        _plot_interval(1),
        _force_def_norm_calc(false)
      {
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \sa \ref solver_configuration
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] plot_name
       * The name of the solver.
       *
       */
      explicit IterativeSolver(const String& plot_name, const String& section_name, const PropertyMap* section) :
        IterativeSolver(plot_name)
      {
        auto plot_mode_p = section->get_entry("plot_mode");
        if (plot_mode_p.second && !plot_mode_p.first.parse(this->_plot_mode))
          throw ParseError(section_name + ".plot_mode", plot_mode_p.first, "one of: none, iter, summary, all");

        auto plot_name_p = section->get_entry("plot_name");
        if (plot_name_p.second)
          this->_plot_name = plot_name_p.first;

        auto plot_interval_p = section->get_entry("plot_interval");
        if (plot_interval_p.second && !plot_interval_p.first.parse(this->_plot_interval))
          throw ParseError(section_name + ".plot_interval", plot_interval_p.first, "an integer >= 0");

        auto tol_abs_p = section->get_entry("tol_abs");
        if (tol_abs_p.second && !tol_abs_p.first.parse(this->_tol_abs))
          throw ParseError(section_name + ".tol_abs", tol_abs_p.first, "a floating point number");

        auto tol_rel_p = section->get_entry("tol_rel");
        if (tol_rel_p.second && !tol_rel_p.first.parse(this->_tol_rel))
          throw ParseError(section_name + ".tol_rel", tol_rel_p.first, "a floating point number");

        auto tol_abs_low_p = section->get_entry("tol_abs_low");
        if (tol_abs_low_p.second && !tol_abs_low_p.first.parse(this->_tol_abs_low))
          throw ParseError(section_name + ".tol_abs_low", tol_abs_low_p.first, "a floating point number");

        auto div_abs_p = section->get_entry("div_abs");
        if (div_abs_p.second && !div_abs_p.first.parse(this->_div_abs))
          throw ParseError(section_name + ".div_abs", div_abs_p.first, "a floating point number");

        auto div_rel_p = section->get_entry("div_rel");
        if (div_rel_p.second && !div_rel_p.first.parse(this->_div_rel))
          throw ParseError(section_name + ".div_rel", div_rel_p.first, "a floating point number");

        auto stag_rate_p = section->get_entry("stag_rate");
        if (stag_rate_p.second && !stag_rate_p.first.parse(this->_stag_rate))
          throw ParseError(section_name + ".stag_rate", stag_rate_p.first, "a floating point number in range [0,1)");

        auto max_iter_p = section->get_entry("max_iter");
        if (max_iter_p.second && !max_iter_p.first.parse(this->_max_iter))
          throw ParseError(section_name + ".max_iter", max_iter_p.first, "an integer >= 0");

        auto min_iter_p = section->get_entry("min_iter");
        if (min_iter_p.second && !min_iter_p.first.parse(this->_min_iter))
          throw ParseError(section_name + ".min_iter", min_iter_p.first, "an integer >= 0");

        auto min_stag_iter_p = section->get_entry("min_stag_iter");
        if (min_stag_iter_p.second && !min_stag_iter_p.first.parse(this->_min_stag_iter))
          throw ParseError(section_name + ".min_stag_iter", min_stag_iter_p.first, "an integer >= 0");
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

      /// Sets the lower absolute tolerance for the solver.
      void set_tol_abs_low(DataType tol_abs_low)
      {
        _tol_abs_low = tol_abs_low;
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

      /// Returns the lower absolute tolerance.
      DataType get_tol_abs_low() const
      {
        return _tol_abs_low;
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

      /// Sets the stagnation rate for the solver.
      void set_stag_rate(DataType rate)
      {
        _stag_rate = rate;
      }

      /// Returns the stagnation rate
      DataType get_stag_rate() const
      {
        return _stag_rate;
      }

      /// Sets the minimum stagnate iteration count for the solver
      void set_min_stag_iter(Index min_iter)
      {
        _min_stag_iter = min_iter;
      }

      /// Returns the minimum stagnation iteration count.
      Index get_min_stag_iter() const
      {
        return _min_stag_iter;
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
        _iter_digits = Math::ilog10(_max_iter);
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
       * \brief Forces the calculation of defect norms in every iteration
       *
       * \note Please note that allowing the skipping of defect norm calculations is merely a hint and a derived class
       * may override this function therefore always forcing the calculation of defect norms. One possible reason is
       * that the derived solver class may require defect norms as part of its solver algorithm and therefore skipping
       * the computation of the defect norms may break the solver; two examples are the BiCGStab and GMRES solvers.
       *
       * \param[in] force
       * Specifies whether defect norm computation is to be enforced in every iteration
       */
      virtual void force_defect_norm_calc(bool force)
      {
        _force_def_norm_calc = force;
      }

      /**
       * \brief Sets the plot mode of the solver.
       *
       * \param[in] plot_mode
       * The desired plot mode
       */
      void set_plot_mode(const PlotMode plot_mode)
      {
        _plot_mode = plot_mode;
      }

      /**
       * \brief Sets the interval between two plot outputs, if any.
       *
       * \param[in] plot_interval The desired interval of iteration plots
       */
      void set_plot_interval(const Index plot_interval)
      {
        _plot_interval = plot_interval;
      }

      /// Sets the plot name of the solver.
      void set_plot_name(const String& plot_name)
      {
        _plot_name = plot_name;
      }

      /// Returns the plot name of the solver.
      String get_plot_name() const
      {
        return _plot_name;
      }

      /**
       * \brief checks for convergence
       *
       * \returns \c true, if converged, else \c false
       */
      bool is_converged() const
      {
        return is_converged(_def_cur);
      }

      /**
       * \brief checks for convergence
       *
       * \param[in] def_cur
       * The defect norm that is to be analyzed
       *
       * \returns \c true, if converged, else \c false
       */
      bool is_converged(const DataType def_cur) const
      {
        return (def_cur <= _tol_abs) && ((def_cur <= (_tol_rel * _def_init)) || (def_cur <= _tol_abs_low));
      }

      /// checks for divergence
      bool is_diverged() const
      {
        return is_diverged(_def_cur);
      }

      /// checks for divergence
      bool is_diverged(const DataType def_cur) const
      {
        return (def_cur > _div_abs) || (def_cur > (_div_rel * _def_init));
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

      /// Returns the status
      Status get_status() const
      {
        return _status;
      }

      /// Computes the overall convergence rate: (defect_final / defect_initial) ^ (1 / number of iterations)
      DataType calc_convergence_rate() const
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

      /// Computes the overall defect reduction factor: (defect_final / defect_inital)
      DataType calc_defect_reduction() const
      {
        // avoid division by zero
        if(this->_def_init <= Math::abs(this->_def_cur * Math::eps<DataType>()))
          return DataType(0);

        // compute defect reduction (def_final / def_inital)
        return this->_def_cur / this->_def_init;
      }

      /**
       * \brief Returns a summary string
       */
      virtual String get_summary() const
      {
        String msg;
        msg += String(  "Name............: ") + this->get_plot_name();
        msg += String("\nStatus..........: ") + stringify(this->get_status());
        msg += String("\nIterations......: ") + stringify(this->get_num_iter());
        msg += String("\nInitial Defect..: ") + stringify_fp_sci(this->get_def_initial());
        msg += String("\nFinal Defect....: ") + stringify_fp_sci(this->get_def_final());
        msg += String("\nDefect Reduction: ") + stringify_fp_sci(this->calc_defect_reduction());
        msg += String("\nConvergence Rate: ") + stringify_fp_fix(this->calc_convergence_rate());

        return msg;
      }

      /**
       * \brief Plot a summary of the last solver run
       */
      virtual void plot_summary() const
      {
        // Print solver summary
        if(!_plot_summary())
          return;

        this->_print_line(this->get_summary());
      }

      /**
       * \brief Solver correction method
       *
       * This method applies the solver represented by this object onto a given right-hand-side vector
       * and updates the corresponding solution vector.
       *
       * In contrast to the apply() method of the SolverBase base class, this method uses the
       * vector \p vec_sol as the initial solution vector for the iterative solution process instead of
       * ignoring its contents upon entry and starting with the null vector.
       *
       * \param[in,out] vec_sol
       * The vector that contains the initial solution upon entry and receives the solution
       * of the linear system upon exit.
       *
       * \param[in] vec_rhs
       * The vector that represents the right-hand-side of the linear system to be solved.
       *
       * \attention vec_sol and vec_rhs must \b not refer to the same vector object!
       *
       * \returns
       * A Status code that represents the status of the solution step.
       */
      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) = 0;

    protected:
      /**
       * \brief Sets the communicator for the solver directly.
       *
       * \param[in] comm
       * A pointer to the communicator that is to be used by the solver.
       */
      void _set_comm(const Dist::Comm* comm)
      {
        this->_comm = comm;
      }

      /**
       * \brief Sets the communicator for the solver from a matrix.
       *
       * \param[in] matrix
       * A reference to a matrix. If 'Matrix_' is a 'Global::Matrix',
       * then the communicator of the matrix gate is taken.
       */
      template<typename Matrix_>
      void _set_comm_by_matrix(const Matrix_& matrix)
      {
        this->_comm = Intern::get_mat_comm(matrix, nullptr);
      }

      /**
       * \brief Sets the communicator for the solver from a vector.
       *
       * \param[in] vector
       * A reference to a vector. If 'Vector_' is a 'Global::Vector',
       * then the communicator of the vector gate is taken.
       */
      void _set_comm_by_vector(const Vector_& vector)
      {
        this->_comm = Intern::get_vec_comm(vector, nullptr);
      }

      /**
       * \brief Computes the defect norm.
       *
       * \param[in] vec_def
       * The current defect vector.
       *
       * \param[in] vec_sol
       * The current solution vector approximation.
       */
      virtual DataType _calc_def_norm(const VectorType& vec_def, const VectorType& DOXY(vec_sol))
      {
        return vec_def.norm2();
      }

      /**
       * \brief Plot the current iteration?
       *
       * \param[in] st
       * The status of the current iteration.
       *
       * \returns \c true if the plot mode is set to \c iter or \c all and the plot interval matches
       * or the solver does not continue to the next iteration
       */
      bool _plot_iter(Status st = Status::progress) const
      {
        return ((_num_iter % _plot_interval == 0) || (st != Status::progress))
          && ((_plot_mode == PlotMode::iter) || (_plot_mode == PlotMode::all));
      }

      /**
       * \brief Plot summary?
       *
       * \returns \c true if the plot mode is set to \c summary or \c all.
       */
      bool _plot_summary() const
      {
        return _plot_mode == PlotMode::summary || _plot_mode == PlotMode::all;
      }

      /**
       * \brief Progress iteration?
       *
       * \returns \c true if the solver should process, otherwise \c false.
       */
      bool _progress() const
      {
        return this->_status == Status::progress;
      }

      /**
       * \brief Prints a line.
       *
       * \param[in] line
       * The line to be printed.
       */
      void _print_line(const String& line) const
      {
        // print message line via comm (if available)
        if(_comm != nullptr)
          _comm->print(line);
        else
          std::cout << line << "\n";
      }

      /**
       * \brief Plots an iteration line.
       *
       * \param[in] num_iter
       * Current number of iterations; usually = this->_num_iter.
       *
       * \param[in] def_cur
       * Current defect norm; usually = this->_def_cur.
       *
       * \param[in] def_prev
       * Previous defect norm; usually = this->_def_prev.
       */
      virtual void _plot_iter_line(Index num_iter, DataType def_cur, DataType def_prev)
      {
        // compose message line
        String msg = this->_plot_name
          +  ": " + stringify(num_iter).pad_front(this->_iter_digits)
          + " : " + stringify_fp_sci(def_cur);

        // not first iteration?
        if(num_iter > Index(0))
        {
          msg += " / " + stringify_fp_sci(def_cur / this->_def_init);
          msg += " / " + stringify_fp_fix(def_cur / def_prev);
        }

        this->_print_line(msg);
      }

      /**
       * \brief Internal function: sets the initial defect vector
       *
       * \param[in] vec_def
       * The initial defect vector.
       *
       * \param[in] vec_sol
       * The current solution vector approximation.
       *
       * \returns
       * The updated Status code.
       */
      virtual Status _set_initial_defect(const VectorType& vec_def, const VectorType& vec_sol)
      {
        // reset iteration variables
        this->_def_init = this->_def_prev = this->_def_cur = DataType(0);
        this->_num_iter = Index(0);
        this->_num_stag_iter = Index(0);

        // first, let's see if we have to compute the defect at all
        bool calc_def = _force_def_norm_calc
          || (this->_min_iter < this->_max_iter)
          || (this->_plot_iter() || this->_plot_summary())
          || (this->_min_stag_iter > Index(0));

        // if we don't have to compute the initial defect, we can exit here
        if(!calc_def)
          return Status::progress;

        // compute new initial defect
        this->_def_init = this->_def_cur = this->_def_prev = this->_calc_def_norm(vec_def, vec_sol);
        Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_init, this->get_num_iter()));

        // plot iteration line?
        if(this->_plot_iter())
          this->_plot_iter_line(Index(0), this->_def_init, this->_def_init);

        // ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(this->_def_init))
          return Status::aborted;

        // check against low absolute tolerance; this may be fulfilled if we have only a handful of dofs
        if(this->_def_init < this->_tol_abs_low)
          return Status::success;

        // check if the initial defect is zero; we test against eps^2 here
        if(this->_def_init <= Math::sqr(Math::eps<DataType>()))
          return Status::success;

        // continue iterating
        return Status::progress;
      }

      /**
       * \brief Internal function: analyze the current defect
       *
       * \note
       * This function is called by _set_new_defect() and _update_defect().
       *
       * \param[in] num_iter
       * Current number of iterations; usually = this->_num_iter.
       *
       * \param[in] def_cur
       * Current defect norm; usually = this->_def_cur.
       *
       * \param[in] def_prev
       * Previous defect norm; usually = this->_def_prev.
       *
       * \param[in] check_stag
       * Specifies whether to check (and update) the stagnation criterion.
       * This is typically set to \c false if one wants to check anything
       * else than the 'true' next defect norm.
       *
       * \returns
       * The updated status code.
       */
      virtual Status _analyze_defect(Index num_iter, DataType def_cur, DataType def_prev, bool check_stag)
      {
        // ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(def_cur))
          return Status::aborted;

        // is diverged?
        if(this->is_diverged(def_cur))
          return Status::diverged;

        // maximum number of iterations performed?
        if(num_iter >= this->_max_iter)
          return Status::max_iter;

        // minimum number of iterations performed?
        if(num_iter < this->_min_iter)
          return Status::progress;

        // is converged?
        if(this->is_converged(def_cur))
          return Status::success;

        // check for stagnation?
        if(check_stag && (this->_min_stag_iter > Index(0)))
        {
          // did this iteration stagnate?
          if(def_cur >= this->_stag_rate * def_prev)
          {
            // increment stagnation count
            if(++this->_num_stag_iter >= this->_min_stag_iter)
              return Status::stagnated;
          }
          else
          {
            // this iteration did not stagnate
            this->_num_stag_iter = Index(0);
          }
        }

        // continue iterating
        return Status::progress;
      }

      /**
       * \brief Internal function: sets the new (next) defect vector
       *
       * This function computes the defect vector's norm, increments the iteration count,
       * plots an output line to std::cout and checks whether any of the stopping criterions is fulfilled.
       *
       * \param[in] vec_def
       * The new defect vector.
       *
       * \param[in] vec_sol
       * The current solution vector approximation.
       *
       * \returns
       * The updated Status code.
       */
      virtual Status _set_new_defect(const VectorType& vec_def, const VectorType& vec_sol)
      {
        // increase iteration count
        ++this->_num_iter;

        // store previous defect
        this->_def_prev = this->_def_cur;

        // first, let's see if we have to compute the defect at all
        bool calc_def = _force_def_norm_calc
          || (this->_min_iter < this->_max_iter)
          || (this->_plot_iter() || this->_plot_summary())
          || (this->_min_stag_iter > Index(0));

        // compute new defect
        if(calc_def)
        {
          this->_def_cur = this->_calc_def_norm(vec_def, vec_sol);
          Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_cur, this->get_num_iter()));
        }

        // analyze defect
        Status status = this->_analyze_defect(this->_num_iter, this->_def_cur, this->_def_prev, true);

        // plot defect?
        if(this->_plot_iter(status))
          this->_plot_iter_line(this->_num_iter, this->_def_cur, this->_def_prev);

        // return status
        return status;
      }

      /**
       * \brief Internal function: sets the new (next) defect norm
       *
       * This function takes a pre-calculated defect vector norm, increments the iteration count,
       * plots an output line to std::cout and checks whether any of the stopping criterions is fulfilled.
       *
       * \param[in] def_cur_norm
       * The new defect norm.
       *
       * \returns
       * The updated Status code.
       *
       * \note This function is preferred over _set_new_defect when using asynchronous mpi operations.
       */
      virtual Status _update_defect(const DataType def_cur_norm)
      {
        // increase iteration count
        ++this->_num_iter;

        // store previous defect
        this->_def_prev = this->_def_cur;

        // update current defect
        this->_def_cur = def_cur_norm;
        Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_cur, this->get_num_iter()));

        // analyze defect
        Status status = this->_analyze_defect(this->_num_iter, this->_def_cur, this->_def_prev, true);

        // plot defect?
        if(this->_plot_iter(status))
          this->_plot_iter_line(this->_num_iter, this->_def_cur, this->_def_prev);

        // return status
        return status;
      }
    }; // class IterativeSolver

    /**
     * \brief Overload of solve() function for IterativeSolver solvers
     */
    template<
      typename Vector_,
      typename Matrix_,
      typename Filter_>
    inline Status solve(
      IterativeSolver<Vector_>& solver,
      Vector_& vec_sol,
      const Vector_& vec_rhs,
      const Matrix_& DOXY(matrix),
      const Filter_& DOXY(filter))
    {
      // simply call the 'correct' method
      return solver.correct(vec_sol, vec_rhs);
    }

    /**
     * \brief Abstract base-class for preconditioned iterative solvers.
     *
     * This class extends the functionality of the IterativeSolver class template by providing overrides
     * for the initialization and finalization methods of the SolverBase class template, which take
     * care of forwarding these steps to the preconditioner.
     *
     * \tparam Vector_
     * The vector type this solver can be applied to
     *
     * \author Peter Zajac
     */
    template<typename Vector_>
    class PreconditionedIterativeSolver :
      public IterativeSolver<Vector_>
    {
    public:
      /// The vector type this solver can be applied to
      typedef Vector_ VectorType;
      /// Floating point data type
      typedef typename VectorType::DataType DataType;
      /// Our base class
      typedef IterativeSolver<VectorType> BaseClass;
      /// the interface for the preconditioner
      typedef SolverBase<VectorType> PrecondType;

    protected:
      /// the pointer to the preconditioner
      std::shared_ptr<PrecondType> _precond;

      /**
       * \brief Constructor
       *
       * \param[in] plot_name
       * The name of the solver.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be nullptr.
       */
      explicit PreconditionedIterativeSolver(const String& plot_name,std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass(plot_name),
        _precond(precond)
      {
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] plot_name
       * The name of the solver.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be nullptr.
       */
      explicit PreconditionedIterativeSolver(const String& plot_name, const String& section_name,
        const PropertyMap* section, std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass(plot_name, section_name, section),
        _precond(precond)
      {
      }

    public:
      /// virtual destructor
      virtual ~PreconditionedIterativeSolver()
      {
      }

      /// \copydoc SolverBase::init_symbolic()
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        if(_precond)
          _precond->init_symbolic();
      }

      /// \copydoc SolverBase::init_numeric()
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        if(_precond)
          _precond->init_numeric();
      }

      /// \copydoc SolverBase::done_numeric()
      virtual void done_numeric() override
      {
        if(_precond)
          _precond->done_numeric();
        BaseClass::done_numeric();
      }

      /// \copydoc SolverBase::done_symbolic()
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
       * \attention vec_cor and vec_def must \b not refer to the same vector object!
       *
       * \param[in,out] vec_cor
       * A reference to the vector that shall receive the preconditioned defect.
       *
       * \param[in] vec_def
       * A reference to the vector that is to be preconditioned.
       *
       * \param[in] filter
       * A reference to the system filter. This filter is only used if no preconditioner is present.
       *
       * \returns
       * \c true, if the preconditioner application was successful, otherwise \c false.
       */
      template<typename Filter_>
      bool _apply_precond(VectorType& vec_cor, const VectorType& vec_def, const Filter_& filter)
      {
        XASSERTM(&vec_cor != &vec_def, "vec_cor and vec_def must not refer to the same object");
        if(this->_precond)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), this->_precond->name()));
          return status_success(this->_precond->apply(vec_cor, vec_def));
        }
        else
        {
          vec_cor.copy(vec_def);
          filter.filter_cor(vec_cor);
          return true;
        }
      }
    }; // class PreconditionedIterativeSolver<...>
  } // namespace Solver
} // namespace FEAT
