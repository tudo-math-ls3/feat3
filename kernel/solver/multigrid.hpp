// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_MULTIGRID_HPP
#define KERNEL_SOLVER_MULTIGRID_HPP 1

// includes, FEAT
#include <kernel/solver/base.hpp>

// includes, system
#include <deque>
#include <vector>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Multigrid Cycle enumeration
     *
     * This enumeration specifies the various cycles supported by the MultiGrid solver.
     */
    enum class MultiGridCycle
    {
      /// Multigrid V-Cycle
      V,
      /// Multigrid F-Cycle
      F,
      /// Multigrid W-Cycle
      W
    };

    /**
     * \brief Multigrid adaptive Coarse-Grid-Correction enumeration
     */
    enum class MultiGridAdaptCGC
    {
      /// fixed coarse grid correction damping
      Fixed,
      /// Energy-Minimization
      MinEnergy,
      /// Defect-Minimization
      MinDefect
    };

    /// \cond internal
    inline std::ostream& operator<<(std::ostream& os, MultiGridCycle cycle)
    {
      switch(cycle)
      {
      case MultiGridCycle::V:
        return os << "V";
      case MultiGridCycle::F:
        return os << "F";
      case MultiGridCycle::W:
        return os << "W";
      default:
        return os << "?";
      }
    }

    inline std::istream& operator>>(std::istream& is, MultiGridCycle& cycle)
    {
      // read character
      char c(' ');
      if((is >> c).fail())
        return is;

      switch(c)
      {
      case 'v':
      case 'V':
        cycle = MultiGridCycle::V;
        break;

      case 'f':
      case 'F':
        cycle = MultiGridCycle::F;
        break;

      case 'w':
      case 'W':
        cycle = MultiGridCycle::W;
        break;

      default:
        // invalid character
        is.putback(c);
        is.setstate(std::ios_base::failbit);
        break;
      }

      return is;
    }
    /// \endcond

    /**
     * \brief Multigrid level base class
     *
     * This class acts as an interface to the MultiGridHierarchy class and is
     * responsible for providing the necessary data for each level.
     *
     * On the coarse mesh level, one has to provide:
     * - the system matrix (mandatory)
     * - the system filter (mandatory)
     * - the coarse-grid solver (optional)
     *
     * On each refined level, one has to provide:
     * - the system matrix (mandatory)
     * - the system filter (mandatory)
     * - prolongation and restriction operators/matrices (mandatory)
     * - pre-/post- and peak-smoothers (optional)
     *
     *
     * \tparam SystemMatrix_
     * The class representing the system matrix.
     *
     * \tparam SystemFilter_
     * The class representing the system filter.
     *
     * \tparam TransferOperator_
     * The class representing the grid-transfer operator.
     *
     * \author Peter Zajac
     */
    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename TransferOperator_>
    class MultiGridLevelBase
    {
    public:
      /// the system matrix type
      typedef SystemMatrix_ SystemMatrixType;
      /// the system filter type
      typedef SystemFilter_ SystemFilterType;
      /// the transfer operator type
      typedef TransferOperator_ TransferOperatorType;
      /// the system vector type
      typedef typename SystemMatrix_::VectorTypeR SystemVectorType;
      /// the coarse-grid solver/smoother type
      typedef SolverBase<SystemVectorType> SolverType;

    public:
      /// virtual destructor
      virtual ~MultiGridLevelBase()
      {
      }

      /**
       * \brief Returns a const reference to the system matrix.
       */
      virtual const SystemMatrixType& get_system_matrix() const = 0;

      /**
       * \brief Returns a const reference to the system filter.
       */
      virtual const SystemFilterType& get_system_filter() const = 0;

      /**
       * \brief Returns a const pointer to the transfer operator.
       *
       * \note
       * This function may return \c nullptr on the coarse-grid level.
       */
      virtual const TransferOperatorType* get_transfer_operator() const = 0;

      /**
       * \brief Returns a shared pointer to the coarse grid solver.
       *
       * \note
       * This function may return \c nullptr on each refined level.
       */
      virtual std::shared_ptr<SolverType> get_coarse_solver() = 0;

      /**
       * \brief Returns a shared pointer to the pre-smoother.
       *
       * \note
       * This function may return \c nullptr if no pre-smoother
       * is to be used on this level.
       */
      virtual std::shared_ptr<SolverType> get_smoother_pre() = 0;

      /**
       * \brief Returns a shared pointer to the post-smoother.
       *
       * \note
       * This function may return \c nullptr if no post-smoother
       * is to be used on this level.
       */
      virtual std::shared_ptr<SolverType> get_smoother_post() = 0;

      /**
       * \brief Returns a shared pointer to the peak-smoother.
       *
       * \note
       * This function may return \c nullptr if both the pre-
       * and post-smoothers are to be used as the peak-smoother.
       */
      virtual std::shared_ptr<SolverType> get_smoother_peak() = 0;
    }; // class MultiGridLevelBase<...>

    /**
     * \brief Standard Multigrid level class implementation
     *
     * This class template implements the MultiGridLevelBase interface by storing
     * references and pointers to the required objects.\n
     * This class has two constructors:
     * - One that accepts all required object references/pointers for the coarse-grid level
     * - One that accepts all required object references/pointers for a refined level
     *
     * \note
     * This class template is used by the MultiGridHierarchy::push_level() functions.
     *
     * \tparam SystemMatrix_
     * The class representing the system matrix.
     *
     * \tparam SystemFilter_
     * The class representing the system filter.
     *
     * \tparam TransferOperator_
     * The class representing the grid-transfer operator.
     *
     * \author Peter Zajac
     */
    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename TransferOperator_>
    class MultiGridLevelStd :
      public MultiGridLevelBase<SystemMatrix_, SystemFilter_, TransferOperator_>
    {
    public:
      /// the base-class type
      typedef MultiGridLevelBase<SystemMatrix_, SystemFilter_, TransferOperator_> BaseClass;
      /// the system-matrix type
      typedef typename BaseClass::SystemMatrixType SystemMatrixType;
      /// the system-filter type
      typedef typename BaseClass::SystemFilterType SystemFilterType;
      /// the system-vector type
      typedef typename BaseClass::SystemVectorType SystemVectorType;
      /// the transfer operator type
      typedef typename BaseClass::TransferOperatorType TransferOperatorType;
      /// the coarse-grid solver/smoother type
      typedef typename BaseClass::SolverType SolverType;

    protected:
      /// the system matrix
      const SystemMatrixType& system_matrix;
      /// the system filter
      const SystemFilterType& system_filter;
      /// the transfer matrix
      const TransferOperatorType* transfer_operator;
      /// the coarse-grid solver
      std::shared_ptr<SolverType> coarse_solver;
      /// the pre-smoother
      std::shared_ptr<SolverType> smoother_pre;
      /// the post-smoother
      std::shared_ptr<SolverType> smoother_post;
      /// the peak-smoother
      std::shared_ptr<SolverType> smoother_peak;

    public:
      /**
       * \brief Coarse-Grid level constructor
       *
       * This constructor is used for the construction of a coarse-grid level.
       *
       * \param[in] sys_matrix
       * A reference to the system matrix.
       *
       * \param[in] sys_filter
       * A reference to the system filter.
       *
       * \param[in] crs_solver
       * A shared pointer to the coarse-grid solver.
       */
      explicit MultiGridLevelStd(
        const SystemMatrixType& sys_matrix,
        const SystemFilterType& sys_filter,
        std::shared_ptr<SolverType> crs_solver)
         :
        system_matrix(sys_matrix),
        system_filter(sys_filter),
        transfer_operator(nullptr),
        coarse_solver(crs_solver),
        smoother_pre(),
        smoother_post(),
        smoother_peak()
      {
      }

      /**
       * \brief Refined level constructor
       *
       * This constructor is used for the construction of a refined level.
       *
       * \param[in] sys_matrix
       * A reference to the system matrix.
       *
       * \param[in] sys_filter
       * A reference to the system filter.
       *
       * \param[in] trans_operat
       * A reference to the transfer operator.
       *
       * \param[in] smooth_pre
       * A shared pointer to the pre-smoother.
       *
       * \param[in] smooth_post
       * A shared pointer to the post-smoother.
       *
       * \param[in] smooth_peak
       * A shared pointer to the peak-smoother.
       * Is set to \c nullptr, the multigrid will use both the pre- and post-smoothers
       * as the peak-smoother.
       *
       * \param[in] crs_solver
       * A shared pointer to the coarse-grid solver, usually a \c nullptr.
       */
      explicit MultiGridLevelStd(
        const SystemMatrixType& sys_matrix,
        const SystemFilterType& sys_filter,
        const TransferOperatorType& trans_operat,
        std::shared_ptr<SolverType> smooth_pre,
        std::shared_ptr<SolverType> smooth_post,
        std::shared_ptr<SolverType> smooth_peak,
        std::shared_ptr<SolverType> crs_solver = nullptr)
         :
        system_matrix(sys_matrix),
        system_filter(sys_filter),
        transfer_operator(&trans_operat),
        coarse_solver(crs_solver),
        smoother_pre(smooth_pre),
        smoother_post(smooth_post),
        smoother_peak(smooth_peak)
      {
      }

      /// virtual destructor
      virtual ~MultiGridLevelStd()
      {
      }

      /// \copydoc MultiGridLevelBase::get_system_matrix()
      virtual const SystemMatrixType& get_system_matrix() const override
      {
        return system_matrix;
      }

      /// \copydoc MultiGridLevelBase::get_system_filter()
      virtual const SystemFilterType& get_system_filter() const override
      {
        return system_filter;
      }

      /// \copydoc MultiGridLevelBase::get_transfer_operator()
      virtual const TransferOperatorType* get_transfer_operator() const override
      {
        return transfer_operator;
      }

      /// \copydoc MultiGridLevelBase::get_coarse_solver()
      virtual std::shared_ptr<SolverType> get_coarse_solver() override
      {
        return coarse_solver;
      }

      /// \copydoc MultiGridLevelBase::get_smoother_pre()
      virtual std::shared_ptr<SolverType> get_smoother_pre() override
      {
        return smoother_pre;
      }

      /// \copydoc MultiGridLevelBase::get_smoother_post()
      virtual std::shared_ptr<SolverType> get_smoother_post() override
      {
        return smoother_post;
      }

      /// \copydoc MultiGridLevelBase::get_smoother_peak()
      virtual std::shared_ptr<SolverType> get_smoother_peak() override
      {
        return smoother_peak;
      }
    }; // class MultiGridLevelStd<...>

    // forward declaration
    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename TransferOperator_>
    class MultiGrid;

    /**
     * \brief Multigrid hierarchy management class template
     *
     * \todo document this
     *
     * \tparam SystemMatrix_
     * The class representing the system matrix.
     *
     * \tparam SystemFilter_
     * The class representing the system filter.
     *
     * \tparam TransferOperator_
     * The class representing the transfer operator.
     *
     * \author Peter Zajac
     */
    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename TransferOperator_>
    class MultiGridHierarchy
    {
    public:
      /// MultiGrid is our friend
      friend class MultiGrid<SystemMatrix_, SystemFilter_, TransferOperator_>;

      /// the level base class type
      typedef MultiGridLevelBase<SystemMatrix_, SystemFilter_, TransferOperator_> LevelType;
      /// the standard level class type
      typedef MultiGridLevelStd<SystemMatrix_, SystemFilter_, TransferOperator_> StdLevelType;

      /// the system matrix type
      typedef typename LevelType::SystemMatrixType SystemMatrixType;
      /// the system filter type
      typedef typename LevelType::SystemFilterType SystemFilterType;
      /// the system vector type
      typedef typename LevelType::SystemVectorType SystemVectorType;
      /// the transfer operator type
      typedef typename LevelType::TransferOperatorType TransferOperatorType;
      /// the coarse-grid solver/smoother type
      typedef typename LevelType::SolverType SolverType;

    protected:
      /**
       * \brief Multigrid level info class
       *
       * This helper class is used internally by the MultiGridHierarchy and MultiGrid classes
       * to manage the various sub-solvers and temporary vectors.
       */
      class LevelInfo
      {
      public:
        /// the actual multigrid level
        std::shared_ptr<LevelType> level;
        /// deque of unique solver object pointers
        std::deque<SolverType*> unique_solvers;
        /// four temporary vectors
        SystemVectorType vec_rhs, vec_sol, vec_def, vec_cor, vec_tmp;

        /// total smoother application time
        double time_smooth;

        /// coarse grid solver application time
        double time_coarse;

        /// defect computation time
        double time_defect;

        /// prolongation/restriction times
        double time_transfer;

        explicit LevelInfo(std::shared_ptr<LevelType> lvl) :
          level(lvl),
          unique_solvers()
        {
          reset_timings();
        }

        void reset_timings()
        {
          time_smooth = time_coarse = time_defect = time_transfer = 0.0;
        }

        void init_symbolic()
        {
          // build unique solver list
          unique_solvers.clear();
          _push_solver(level->get_coarse_solver().get());
          _push_solver(level->get_smoother_pre().get());
          _push_solver(level->get_smoother_post().get());
          _push_solver(level->get_smoother_peak().get());

          // initialize all unique solvers in forward order
          for(auto it = unique_solvers.begin(); it != unique_solvers.end(); ++it)
            (*it)->init_symbolic();

          // get level system matrix
          const SystemMatrixType& matrix = level->get_system_matrix();

          // create our vectors
          vec_rhs = matrix.create_vector_r();
          vec_sol = matrix.create_vector_r();
          vec_def = matrix.create_vector_r();
          vec_cor = matrix.create_vector_r();
          vec_tmp = matrix.create_vector_r();
        }

        void done_symbolic()
        {
          // release our vectors
          vec_tmp.clear();
          vec_cor.clear();
          vec_def.clear();
          vec_sol.clear();
          vec_rhs.clear();

          // release all unique solvers in reverse order
          for(auto it = unique_solvers.rbegin(); it != unique_solvers.rend(); ++it)
            (*it)->done_symbolic();
        }

        void init_numeric()
        {
          // initialize all unique solvers in forward order
          for(auto it = unique_solvers.begin(); it != unique_solvers.end(); ++it)
            (*it)->init_numeric();
        }

        void done_numeric()
        {
          // release all unique solvers in reverse order
          for(auto it = unique_solvers.rbegin(); it != unique_solvers.rend(); ++it)
            (*it)->done_numeric();
        }

      protected:
        void _push_solver(SolverType* solver)
        {
          if(solver == nullptr)
            return;

          // loop over all existing solvers to determine uniqueness
          for(auto it = unique_solvers.begin(); it != unique_solvers.end(); ++it)
          {
            // do we already have this solver object?
            if(*it == solver)
              return;
          }

          // alright, that's a new unique solver object, so add it
          unique_solvers.push_back(solver);
        }
      }; // class LevelInfo

      /**
       * \brief Returns a level info in the hierarchy.
       *
       * \note
       * This function is used by the MultiGrid solver class.
       *
       * \param[in] lvl
       * The index of the level info to be returned.
       *
       * \returns
       * A reference to the corresponding level info object.
       */
      LevelInfo& _get_level_info(Index lvl)
      {
        XASSERTM(lvl < size_physical(), "invalid level index");
        return _levels.at(std::size_t(lvl));
      }

    protected:
      /// the deque of all level info objects
      std::deque<LevelInfo> _levels;
      /// the virtual multigrid hierarchy sizes
      std::size_t _virt_size;
      /// specifies whether init_symbolic() has been called
      bool _have_init_symbolic;
      /// specifies whether init_numeric() has been called
      bool _have_init_numeric;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] size_virt
       * The virtual multigrid hierarchy size on all processes.
       */
      explicit MultiGridHierarchy(std::size_t size_virt) :
        _levels(),
        _virt_size(size_virt),
        _have_init_symbolic(false),
        _have_init_numeric(false)
      {
      }

      /**
       * \brief Destructor
       */
      virtual ~MultiGridHierarchy()
      {
        // delete levels in reverse order
        while(!_levels.empty())
          _levels.pop_back();
      }

      /**
       * \brief Pushes a new (custom) level into the hierarchy.
       *
       * \param[in] level
       * A shared pointer to the new level.
       */
      void push_level(std::shared_ptr<LevelType> level)
      {
        XASSERTM(!_have_init_symbolic, "cannot push new level, init_symbolic() was already called");
        XASSERTM(!_have_init_numeric,  "cannot push new level, init_numeric() was already called");
        XASSERTM(_levels.size() < _virt_size, "cannot push new level, already have maximum number of levels");
        _levels.push_back(LevelInfo(level));
      }

      /**
       * \brief Pushes a standard coarse-grid level into the hierarchy.
       *
       * \param[in] system_matrix
       * A reference to the system matrix.
       *
       * \param[in] system_filter
       * A reference to the system filter.
       *
       * \param[in] coarse_solver
       * A shared pointer to the coarse grid solver.
       * May be \c nullptr if no coarse-grid solver is to be used.
       */
      void push_level(
        const SystemMatrixType& system_matrix,
        const SystemFilterType& system_filter,
        std::shared_ptr<SolverType> coarse_solver)
      {
        push_level(std::make_shared<StdLevelType>(system_matrix, system_filter, coarse_solver));
      }

      /**
       * \brief Pushes a standard refined level into the hierarchy.
       *
       * \param[in] system_matrix
       * A reference to the system matrix.
       *
       * \param[in] system_filter
       * A reference to the system filter.
       *
       * \param[in] transfer_operator
       * A reference to the transfer operator onto this level.
       *
       * \param[in] smoother_pre
       * A shared pointer to the pre-smoother.
       * May be \c nullptr if no pre-smoother is to be used.
       *
       * \param[in] smoother_post
       * A shared pointer to the post-smoother.
       * May be \c nullptr if no post-smoother is to be used.
       *
       * \param[in] smoother_peak
       * A shared pointer to the peak-smoother.
       * Is set to \c nullptr, the multigrid will use both the pre- and post-smoothers
       * as the peak-smoother.
       *
       * \param[in] coarse_solver
       * A shared pointer to the coarse grid solver.
       * This one is usually \c nullptr, unless one wants to use this refined
       * level also as a coarse-grid level.
       *
       * \note
       * The \p smoother_pre, \p smoother_post, \p smoother_peak and \p coarse_solver
       * pointers need not point to different objects. This class will ensure the correct
       * handling of unique solver objects independent of whether more that one pointer
       * references the same unique object on one level.
       */
      void push_level(
        const SystemMatrixType& system_matrix,
        const SystemFilterType& system_filter,
        const TransferOperatorType& transfer_operator,
        std::shared_ptr<SolverType> smoother_pre,
        std::shared_ptr<SolverType> smoother_post,
        std::shared_ptr<SolverType> smoother_peak,
        std::shared_ptr<SolverType> coarse_solver = nullptr
        )
      {
        push_level(std::make_shared<StdLevelType>(system_matrix, system_filter, transfer_operator,
          smoother_pre, smoother_post, smoother_peak, coarse_solver));
      }

      /**
       * \brief Returns the physical hierarchy size.
       *
       * \returns The physical hierarchy size.
       */
      Index size_physical() const
      {
        return Index(_levels.size());
      }

      /**
       * \brief Returns the virtual hierarchy size.
       *
       * \returns The virtual hierarchy size.
       */
      Index size_virtual() const
      {
        return Index(_virt_size);
      }

      /**
       * \brief Returns a level in the hierarchy.
       *
       * \param[in] lvl
       * The index of the level to be returned.
       *
       * \returns
       * A (const) reference to the corresponding level object.
       */
      LevelType& get_level(Index lvl)
      {
        XASSERTM(lvl < size_physical(), "invalid level index");
        return *(_levels.at(std::size_t(lvl)).level);
      }

      /** \copydoc get_level() */
      const LevelType& get_level(Index lvl) const
      {
        XASSERTM(lvl < size_physical(), "invalid level index");
        return *(_levels.at(std::size_t(lvl)).level);
      }

      /**
       * \brief Symbolic initialization method
       *
       * This method is called to perform symbolic initialization of the sub-solvers
       * in the multigrid hierarchy.
       */
      void init_symbolic()
      {
        XASSERTM(!_have_init_numeric, "cannot call init_symbolic(): init_numeric() has already been called");
        XASSERTM(!_have_init_symbolic, "cannot call init_symbolic(): init_symbolic() has already been called");

        // initialize all levels in forward order
        for(auto it = _levels.begin(); it != _levels.end(); ++it)
          it->init_symbolic();

        _have_init_symbolic = true;
      }

      /**
       * \brief Symbolic finalization method
       *
       * This method is called to release any data allocated in the symbolic initialization step
       * of the sub-solvers in the multigrid hierarchy.
       */
      void done_symbolic()
      {
        XASSERTM(!_have_init_numeric, "cannot call done_symbolic(): done_numeric() has not been called yet");
        XASSERTM(_have_init_symbolic, "cannot call done_symbolic(): init_symbolic() has not been called yet");

        // release all levels in reverse order
        for(auto it = _levels.rbegin(); it != _levels.rend(); ++it)
          it->done_symbolic();

        _have_init_symbolic = false;
      }

      /**
       * \brief Numeric initialization method
       *
       * This method is called to perform numeric initialization of the sub-solvers
       * in the multigrid hierarchy.\n
       * Before this function can be called, the symbolic initialization must be performed.
       */
      void init_numeric()
      {
        XASSERTM(!_have_init_numeric, "cannot call init_numeric(): init_numeric() has already been called");
        XASSERTM(_have_init_symbolic, "cannot call init_numeric(): init_symbolic() has not been called yet");

        // initialize all levels in forward order
        for(auto it = _levels.begin(); it != _levels.end(); ++it)
          it->init_numeric();

        _have_init_numeric = true;
      }

      /**
       * \brief Numeric finalization method
       *
       * This method is called to release any data allocated in the numeric initialization step
       * of the sub-solvers in the multigrid hierarchy.
       */
      void done_numeric()
      {
        XASSERTM(_have_init_numeric, "cannot call done_numeric(): init_numeric() has not been called yet");

        // release all levels in reverse order
        for(auto it = _levels.rbegin(); it != _levels.rend(); ++it)
          it->done_numeric();

        _have_init_numeric = false;
      }

      /**
       * \brief Initialization method
       *
       * This function performs both the symbolic and numeric initialization, i.e. it simply performs
         \verbatim
         this->init_symbolic();
         this->init_numeric();
         \endverbatim
       */
      void init()
      {
        init_symbolic();
        init_numeric();
      }

      /**
       * \brief Finalization method
       *
       * This function performs both the symbolic and numeric finalization, i.e. it simply performs
         \verbatim
         this->done_numeric();
         this->done_symbolic();
         \endverbatim
       */
      void done()
      {
        done_numeric();
        done_symbolic();
      }

      /**
       * \brief Specifies whether the symbolic initialization was performed.
       *
       * \returns
       * \c true between a call of init_symbolic() and done_symbolic(), otherwise \c false.
       */
      bool have_symbolic() const
      {
        return _have_init_symbolic;
      }

      /**
       * \brief Specifies whether the symbolic initialization was performed.
       *
       * \returns
       * \c true between a call of init_symbolic() and done_symbolic(), otherwise \c false.
       */
      bool have_numeric() const
      {
        return _have_init_numeric;
      }

      /**
       * \brief Resets the internal timing statistics.
       */
      void reset_timings()
      {
        for(auto& li : _levels)
          li.reset_timings();
      }

      /**
       * \brief Returns the time for smoother application.
       *
       * \param[in] lvl
       * The level for which the timing is to be retrieved.
       * If set to -1, the sum of all level timings is returned.
       *
       * \returns The time for smoother application.
       */
      double get_time_smooth(int lvl = -1) const
      {
        if(lvl < 0)
        {
          double t(0.0);
          for(const auto& li : _levels)
            t += li.time_smooth;
          return t;
        }

        XASSERT(lvl < int(_levels.size()));
        return _levels.at(std::size_t(lvl)).time_smooth;
      }

      /**
       * \brief Returns the time for coarse grid solver application.
       *
       * \param[in] lvl
       * The level for which the timing is to be retrieved.
       * If set to -1, the sum of all level timings is returned.
       *
       * \returns The time for coarse grid solver application.
       */
      double get_time_coarse(int lvl = -1) const
      {
        if(lvl < 0)
        {
          double t(0.0);
          for(const auto& li : _levels)
            t += li.time_coarse;
          return t;
        }

        XASSERT(lvl < int(_levels.size()));
        return _levels.at(std::size_t(lvl)).time_coarse;
      }

      /**
       * \brief Returns the time for defect computation.
       *
       * \param[in] lvl
       * The level for which the timing is to be retrieved.
       * If set to -1, the sum of all level timings is returned.
       *
       * \returns The time for defect computation.
       */
      double get_time_defect(int lvl = -1) const
      {
        if(lvl < 0)
        {
          double t(0.0);
          for(const auto& li : _levels)
            t += li.time_defect;
          return t;
        }

        XASSERT(lvl < int(_levels.size()));
        return _levels.at(std::size_t(lvl)).time_defect;
      }

      /**
       * \brief Returns the time for grid transfer application.
       *
       * \param[in] lvl
       * The level for which the timing is to be retrieved.
       * If set to -1, the sum of all level timings is returned.
       *
       * \returns The time for grid transfer application.
       */
      double get_time_transfer(int lvl = -1) const
      {
        if(lvl < 0)
        {
          double t(0.0);
          for(const auto& li : _levels)
            t += li.time_transfer;
          return t;
        }

        XASSERT(lvl < int(_levels.size()));
        return _levels.at(std::size_t(lvl)).time_transfer;
      }
    }; // class MultiGridHierarchy<...>

    /**
     * \brief Multigrid preconditioner implementation
     *
     * \sa For details in the design see \ref multigrid_design
     *
     * \todo document this
     *
     * \tparam SystemMatrix_
     * The class representing the system matrix.
     *
     * \tparam SystemFilter_
     * The class representing the system filter.
     *
     * \tparam ProlOperator_
     * The class representing the prolongation operator (usually a matrix type).
     *
     * \tparam RestOperator_
     * The class representing the restriction operator (usually a matrix type).
     *
     * \author Peter Zajac
     */
    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename TransferOperator_>
    class MultiGrid :
      public SolverBase<typename SystemMatrix_::VectorTypeR>
    {
    public:
      /// our base-class
      typedef SolverBase<typename SystemMatrix_::VectorTypeR> BaseClass;

      /// our compatible multigrid hierarchy class
      typedef MultiGridHierarchy<SystemMatrix_, SystemFilter_, TransferOperator_> HierarchyType;

      /// the level type
      typedef typename HierarchyType::LevelType LevelType;
      /// the level info type
      typedef typename HierarchyType::LevelInfo LevelInfo;

      /// the system matrix type
      typedef typename LevelType::SystemMatrixType MatrixType;
      /// the system filter type
      typedef typename LevelType::SystemFilterType FilterType;
      /// the system vector type
      typedef typename LevelType::SystemVectorType VectorType;
      /// the transfer operator type
      typedef typename LevelType::TransferOperatorType TransferOperatorType;
      /// the sub-solver type
      typedef typename LevelType::SolverType SolverType;

      /// our data type
      typedef typename MatrixType::DataType DataType;

    protected:
      /// the multigrid hierarchy object
      std::shared_ptr<HierarchyType> _hierarchy;
      /// the multigrid cycle
      MultiGridCycle _cycle;
      /// the coarse grid correction type
      MultiGridAdaptCGC _adapt_cgc;
      /// the top-level of this multigrid
      Index _top_level;
      /// the coarse-level of this multigrid
      Index _crs_level;

      /// W-cycle level counter vector
      std::vector<int> _counters;

      /// array containing toe for each processed level
      std::vector<double> _toes;
      /// array containing toe of mpi execution reduction for each processed level
      std::vector<double> _mpi_execs_reduction;
      /// array containing toe of mpi execution blas2 for each processed level
      std::vector<double> _mpi_execs_blas2;
      /// array containing toe of mpi execution blas3 for each processed level
      std::vector<double> _mpi_execs_blas3;
      /// array containing toe of mpi execution collective for each processed level
      std::vector<double> _mpi_execs_collective;
      /// array containing toe of mpi wait reduction for each processed level
      std::vector<double> _mpi_waits_reduction;
      /// array containing toe of mpi wait blas2 for each processed level
      std::vector<double> _mpi_waits_blas2;
      /// array containing toe of mpi wait blas3 for each processed level
      std::vector<double> _mpi_waits_blas3;
      /// array containing toe of mpi wait collective for each processed level
      std::vector<double> _mpi_waits_collective;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] hierarchy
       * A pointer to the multigrid hierarchy object.
       *
       * \param[in] cycle
       * The desired multigrid cycle.
       *
       * \param[in] top_level
       * The desired top-level for this multigrid solver.\n
       * Set to 0 to use the finest level in the multigrid hierarchy.
       *
       * \param[in] crs_level
       * The desired coarse-level for this multigrid solver.\n
       * Set to -1 to use the coarsest level in the multigrid hierarchy.
       *
       * \note
       * The cycle may be changed anytime by using the set_cycle() function.
       */
      explicit MultiGrid(std::shared_ptr<HierarchyType> hierarchy, MultiGridCycle cycle, int top_level = 0, int crs_level = -1) :
        BaseClass(),
        _hierarchy(hierarchy),
        _cycle(cycle),
        _adapt_cgc(MultiGridAdaptCGC::Fixed),
        _top_level(Index(top_level)),
        _crs_level(crs_level >= 0 ? Index(crs_level) : Index(int(_hierarchy->size_virtual()) + crs_level))
      {
        XASSERTM(_crs_level < _hierarchy->size_virtual(), "invalid coarse level");
        XASSERTM(_top_level <= _crs_level, "invalid top/coarse level combination");
      }

      /// virtual destructor
      virtual ~MultiGrid()
      {
      }

      /**
       * \brief Sets the adaption mode for the coarse grid correction.
       *
       * \param[in] adapt_cgc
       * The adaption mode for the coarse grid correction.
       */
      void set_adapt_cgc(MultiGridAdaptCGC adapt_cgc)
      {
        _adapt_cgc = adapt_cgc;
      }

      /**
       * \brief Sets a new cycle.
       *
       * \param[in] cycle
       * The new cycle to be used.
       */
      void set_cycle(MultiGridCycle cycle)
      {
        _cycle = cycle;
      }

      /**
       * \brief Returns the currently selected cycle.
       *
       * \returns The currently selected cycle.
       */
      MultiGridCycle get_cycle() const
      {
        return _cycle;
      }

      /**
       * \brief Sets the level range for this multigrid.
       *
       * \param[in] top_level
       * The desired top-level for this multigrid solver.\n
       * Set to 0 to use the finest level in the multigrid hierarchy.
       *
       * \param[in] crs_level
       * The desired coarse-level for this multigrid solver.\n
       * Set to -1 to use the coarsest level in the multigrid hierarchy.
       */
      void set_levels(int top_level, int crs_level)
      {
        _top_level = Index(top_level);
        _crs_level = (crs_level >= 0 ? Index(crs_level) : Index(int(_hierarchy->size_virtual()) + crs_level));

        XASSERTM(_crs_level < _hierarchy->size_virtual(), "invalid coarse level");
        XASSERTM(_top_level <= _crs_level, "invalid top/coarse level combination");
      }

      /**
       * \brief Returns the currently selected top-level index.
       *
       * \returns The currently selected top-level index.
       */
      Index get_top_level() const
      {
        return _top_level;
      }

      /**
       * \brief Returns the currently selected coarse-level index.
       *
       * \returns The currently selected coarse-level index.
       */
      Index get_crs_level() const
      {
        return _crs_level;
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the solver.
       */
      virtual String name() const override
      {
        return "MultiGrid-" + stringify(_cycle);
      }

      /**
       * \brief Symbolic initialization function.
       *
       * \attention
       * This function does \b not perform the actual symbolic initialization
       * of the multigrid hierarchy. You have to perform the symbolic initialization
       * explicitly by calling the init_symbolic function of the corresponding
       * MultiGridHierarchy object \b before initialising the solver tree!
       */
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // ensure that the hierarchy was initialized
        XASSERTM(_hierarchy->have_symbolic(), "init_symbolic() of multigrid hierarchy was not called yet");

        // get the virtual hierarchy size
        const std::size_t num_lvls = _hierarchy->size_virtual();

        // allocate counter vector
        _counters.resize(num_lvls, 0);

        // allocate statistics vectors
        _toes.resize(num_lvls, 0.0);
        _mpi_execs_reduction.resize(num_lvls, 0.0);
        _mpi_execs_blas2.resize(num_lvls, 0.0);
        _mpi_execs_blas3.resize(num_lvls, 0.0);
        _mpi_execs_collective.resize(num_lvls, 0.0);
        _mpi_waits_reduction.resize(num_lvls, 0.0);
        _mpi_waits_blas2.resize(num_lvls, 0.0);
        _mpi_waits_blas3.resize(num_lvls, 0.0);
        _mpi_waits_collective.resize(num_lvls, 0.0);
      }

      /**
       * \brief Symbolic finalization function.
       *
       * \attention
       * This function does \b not perform the actual symbolic finalization
       * of the multigrid hierarchy. You have to perform the symbolic finalization
       * explicitly by calling the done_symbolic function of the corresponding
       * MultiGridHierarchy object \b after finalising the solver tree!
       */
      virtual void done_symbolic() override
      {
        _toes.clear();
        _mpi_execs_reduction.clear();
        _mpi_execs_blas2.clear();
        _mpi_execs_blas3.clear();
        _mpi_execs_collective.clear();
        _mpi_waits_reduction.clear();
        _mpi_waits_blas2.clear();
        _mpi_waits_blas3.clear();
        _mpi_waits_collective.clear();
        _counters.clear();

        BaseClass::done_symbolic();
      }

      /**
       * \brief Numeric initialization function.
       *
       * \attention
       * This function does \b not perform the actual numeric initialization
       * of the multigrid hierarchy. You have to perform the numeric initialization
       * explicitly by calling the init_numeric function of the corresponding
       * MultiGridHierarchy object \b before initialising the solver tree!
       */
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // ensure that the hierarchy was initialized
        XASSERTM(_hierarchy->have_numeric(), "init_numeric() of multigrid hierarchy was not called yet");
      }

      /**
       * \brief Numeric finalization function.
       *
       * \attention
       * This function does \b not perform the actual numeric finalization
       * of the multigrid hierarchy. You have to perform the numeric finalization
       * explicitly by calling the done_numeric function of the corresponding
       * MultiGridHierarchy object \b after finalising the solver tree!
       */
      virtual void done_numeric() override
      {
        BaseClass::done_numeric();
      }

      /**
       * \brief Applies the multigrid preconditioner.
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
      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        // reset statistics counters
        for(std::size_t i(0); i <  std::size_t(_hierarchy->size_virtual()); ++i)
        {
          _toes.at(i) = 0.0;
          _mpi_execs_reduction.at(i) = 0.0;
          _mpi_execs_blas2.at(i) = 0.0;
          _mpi_execs_blas3.at(i) = 0.0;
          _mpi_execs_collective.at(i) = 0.0;
          _mpi_waits_reduction.at(i) = 0.0;
          _mpi_waits_blas2.at(i) = 0.0;
          _mpi_waits_blas3.at(i) = 0.0;
          _mpi_waits_collective.at(i) = 0.0;
        }

        // get a reference to the finest level
        LevelInfo& lvl_top = _hierarchy->_get_level_info(_top_level);

        // copy defect to RHS vector
        lvl_top.vec_rhs.copy(vec_def);

        // apply the corresponding cycle
        Status status(Status::undefined);
        switch(_cycle)
        {
        case MultiGridCycle::V:
          status = this->_apply_cycle_v();
          break;

        case MultiGridCycle::F:
          status = this->_apply_cycle_f();
          break;

        case MultiGridCycle::W:
          status = this->_apply_cycle_w();
          break;

        default:
          // whoops...
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, 1));
          status = Status::aborted;
          break;
        }

        // copy solution to correction vector
        vec_cor.copy(lvl_top.vec_sol);

        // propagate solver statistics
        for(std::size_t i(0); i <  std::size_t(_hierarchy->size_virtual()); ++i)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionLevelTimings>(this->name(), Index(i),
            _toes.at(i), _mpi_execs_reduction.at(i), _mpi_execs_blas2.at(i), _mpi_execs_blas3.at(i), _mpi_execs_collective.at(i), _mpi_waits_reduction.at(i), _mpi_waits_blas2.at(i),
            _mpi_waits_blas3.at(i), _mpi_waits_collective.at(i)));
        }

        // okay
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, 1));
        return status;
      }

    protected:
      /// Applies a single V-cycle.
      Status _apply_cycle_v()
      {
        // The V-cycle is pretty simple:
        // 1) Restrict from the top level to the coarse level
        // 2) Solve coarse level system
        // 3) Prolongate from the coarse level to the top level

        // restrict (apply top-level pre-smoother if it exists)
        this->_apply_rest(_top_level, true);

        // solve coarse level
        this->_apply_coarse();

        // prolongate (apply top-level post-smoother if it exists)
        this->_apply_prol(_top_level, true);

        // okay
        return Status::success;
      }

      /// Applies a single F-cycle.
      Status _apply_cycle_f()
      {
        // The F-cycle is slightly more complex than the V-cycle:
        // 1) Restrict from the top level to the coarse level
        //
        // 2) For each intermediate level (that's everything between the
        //    coarse and top level) in ascending order:
        //    2.1) solve the coarse system
        //    2.2) prolongate to the intermediate level without applying
        //         the post-smoother on the intermediate level
        //    2.3) apply peak-smoother
        //    2.4) restrict to the coarse level without applying
        //         the pre-smoother on the intermediate level
        //
        // 3) Solve the coarse level system
        //
        // 4) Prolongate from the coarse level to the top level

        // restrict (apply top-level pre-smoother if it exists)
        this->_apply_rest(_top_level, true);

        // get the coarsest level on this process; that's either the
        // coarse level or the last level for this process:
        const Index last_level = Math::min(_crs_level, _hierarchy->size_physical());

        // ensure that we do not apply this if there is just one level
        if(last_level > Index(0))
        {
          // F-cycle intermediate peak-level loop
          for(Index peak_lvl(last_level-1); peak_lvl > _top_level; --peak_lvl)
          {
            // solve coarse level
            this->_apply_coarse();

            // prolongate to current peak-level without post-smoothing
            this->_apply_prol(peak_lvl, false);

            // apply peak-smoother
            this->_apply_smooth_peak(peak_lvl);

            // restrict from current peak-level without pre-smoothing
            this->_apply_rest(peak_lvl, false);
          }
        }

        // solve coarse level
        this->_apply_coarse();

        // prolongate (apply top-level post-smoother if it exists)
        this->_apply_prol(_top_level, true);

        // okay
        return Status::success;
      }

      /// Applies a single W-cycle.
      Status _apply_cycle_w()
      {
        // The W-cycle is a tough nut to implement/understand:
        // As always, we start by restricting from the top level to the
        // coarse level and solving the coarse level system afterwards.
        //
        // The next part is the tricky bit:
        // We need to implement the "inner W-cycle part". In theory,
        // this is usually defined as a "recursive 2x V-cycle application",
        // but we do not want to use recursion here explicitly.
        // We know that the W-cycle performs 2^(top_level - crs_level)
        // coarse grid solves, so we use that for our for-loop.
        // Moreover, we manage a counter for each level above the coarse
        // level and count how often we have applied the peak-smoother.
        // Each iteration we choose the lowest level that has a peak-smooth
        // count equal to 0 as the next intermediate level, increment
        // its counter by 1 and reset all level counters below this level
        // back to 0. This may sound perfectly confusing, but it actually
        // does the job :) See the description of the multigrid cycles
        // in the documentation for further explanations and a fancy image
        // of the W-cycle level counters in action.
        //
        // Oh yes, and finally we also need to prolongate from the
        // coarse level to the top level, of course.

        // restrict (apply top-level pre-smoother if it exists)
        this->_apply_rest(_top_level, true);

        // solve coarse level
        this->_apply_coarse();

        // get the coarsest level on this process; that's either the
        // coarse level or the last level for this process:
        const Index last_level = Math::min(_crs_level, _hierarchy->size_physical());

        // reset all level counters to 0
        for(std::size_t i(_top_level); i <= std::size_t(last_level); ++i)
        {
          _counters[i] = 0;
        }

        // compute the total number of coarse-grid solves to perform:
        // that's = 2^(top_level - crs_level)
        const int num_cgs = 1 << (last_level - _top_level);

        // W-cycle loop: count the number of coarse-grid solves
        for(int cgs(1); cgs < num_cgs; ++cgs)
        {
          // find the next peak level; that's the lowest level
          // above the coarse level whose counter is 0
          std::size_t peak_lvl = std::size_t(last_level);
          while(peak_lvl > std::size_t(_top_level))
          {
            // if the counter is 0, we have our next peak level
            if(_counters[--peak_lvl] == 0)
              break;
          }

          // >>>>> DEBUG >>>>>
          //std::cout << stringify(cgs).pad_front(3) << " ";
          //std::cout << stringify(peak_lvl).pad_front(2) << ":";
          //for(std::size_t i(0); i <= std::size_t(top_lvl); ++i)
          //  std::cout << " " << _counters[i];
          // <<<<< DEBUG <<<<<

          // sanity check: do not go beyond our top level
          XASSERTM(peak_lvl >= _top_level, "W-cycle sanity check failed: invalid peak level");

          // reset all counters below our peak level to 0
          for(std::size_t i = std::size_t(last_level-1); i > peak_lvl; --i)
            _counters[i] = 0;

          // increment the counter of our peak level by 1,
          // indicating that we are going to perform an
          // peak-smoother step on this level
          ++_counters[peak_lvl];

          // >>>>> DEBUG >>>>>
          //std::cout << " >";
          //for(std::size_t i(0); i <= std::size_t(top_lvl); ++i)
          //  std::cout << " " << _counters[i];
          //std::cout << std::endl;
          // <<<<< DEBUG <<<<<

          // prolongate to current peak-level without post-smoothing
          this->_apply_prol(Index(peak_lvl), false);

          // apply peak-smoother
          this->_apply_smooth_peak(Index(peak_lvl));

          // restrict from current peak-level without pre-smoothing
          this->_apply_rest(Index(peak_lvl), false);

          // solve coarse level
          this->_apply_coarse();
        }

        // sanity check: all level counters above the coarse level must be 1
        for(std::size_t i(last_level); i > std::size_t(_top_level); )
        {
          --i;
          XASSERTM(this->_counters[i] == 1, "W-cycle sanity check failed: invalid counters");
        }

        // prolongate (apply top-level post-smoother if it exists)
        this->_apply_prol(_top_level, true);

        // okay
        return Status::success;
      }

      /**
       * \brief Applies the coarse grid solver on the coarse level.
       */
      Status _apply_coarse()
      {
        // no coarse level on this process?
        if(_crs_level >= _hierarchy->size_physical())
          return Status::success;

        // process the coarse grid level
        TimeStamp at;
        double mpi_exec_reduction_start(Statistics::get_time_mpi_execute_reduction());
        double mpi_exec_blas2_start(Statistics::get_time_mpi_execute_blas2());
        double mpi_exec_blas3_start(Statistics::get_time_mpi_execute_blas3());
        double mpi_exec_collective_start(Statistics::get_time_mpi_execute_collective());
        double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
        double mpi_wait_start_blas2(Statistics::get_time_mpi_wait_blas2());
        double mpi_wait_start_blas3(Statistics::get_time_mpi_wait_blas3());
        double mpi_wait_start_collective(Statistics::get_time_mpi_wait_collective());

        // get the coarse level info
        LevelInfo& lvl_crs = _hierarchy->_get_level_info(_crs_level);

        // get the system filter
        const FilterType& system_filter = lvl_crs.level->get_system_filter();

        // get the coarse grid solver
        std::shared_ptr<SolverType> coarse_solver = lvl_crs.level->get_coarse_solver();

        // if the have a coarse grid solver, apply it
        TimeStamp stamp_coarse;
        if(coarse_solver)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionCallCoarseSolver>(this->name(), coarse_solver->name()));
          if(!status_success(coarse_solver->apply(lvl_crs.vec_sol, lvl_crs.vec_rhs)))
            return Status::aborted;
        }
        else
        {
          // simply copy the RHS thus emulating an identity solver
          lvl_crs.vec_sol.copy(lvl_crs.vec_rhs);

          // apply the correction filter
          system_filter.filter_cor(lvl_crs.vec_sol);
        }
        lvl_crs.time_coarse += stamp_coarse.elapsed_now();

        _toes.at(std::size_t(_crs_level)) += at.elapsed_now();
        double mpi_exec_reduction_stop(Statistics::get_time_mpi_execute_reduction());
        double mpi_exec_blas2_stop(Statistics::get_time_mpi_execute_blas2());
        double mpi_exec_blas3_stop(Statistics::get_time_mpi_execute_blas3());
        double mpi_exec_collective_stop(Statistics::get_time_mpi_execute_collective());
        double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
        double mpi_wait_stop_blas2(Statistics::get_time_mpi_wait_blas2());
        double mpi_wait_stop_blas3(Statistics::get_time_mpi_wait_blas3());
        double mpi_wait_stop_collective(Statistics::get_time_mpi_wait_collective());
        _mpi_execs_reduction.at(std::size_t(_crs_level)) += mpi_exec_reduction_stop - mpi_exec_reduction_start;
        _mpi_execs_blas2.at(std::size_t(_crs_level)) += mpi_exec_blas2_stop - mpi_exec_blas2_start;
        _mpi_execs_blas3.at(std::size_t(_crs_level)) += mpi_exec_blas3_stop - mpi_exec_blas3_start;
        _mpi_execs_collective.at(std::size_t(_crs_level)) += mpi_exec_collective_stop - mpi_exec_collective_start;
        _mpi_waits_reduction.at(std::size_t(_crs_level)) += mpi_wait_stop_reduction - mpi_wait_start_reduction;
        _mpi_waits_blas2.at(std::size_t(_crs_level)) += mpi_wait_stop_blas2 - mpi_wait_start_blas2;
        _mpi_waits_blas3.at(std::size_t(_crs_level)) += mpi_wait_stop_blas3 - mpi_wait_start_blas3;
        _mpi_waits_collective.at(std::size_t(_crs_level)) += mpi_wait_stop_collective - mpi_wait_start_collective;

        return Status::success;
      }

      bool _apply_smooth_def(Index cur_lvl, SolverType& smoother)
      {
        // get the level info
        LevelInfo& lvl = _hierarchy->_get_level_info(cur_lvl);

        // get the system matrix and filter
        const MatrixType& system_matrix = lvl.level->get_system_matrix();
        const FilterType& system_filter = lvl.level->get_system_filter();

        // apply peak-smoother
        TimeStamp stamp_smooth;
        Statistics::add_solver_expression(std::make_shared<ExpressionCallSmoother>(this->name(), smoother.name()));
        smoother.apply(lvl.vec_cor, lvl.vec_def);
        //if(!status_success(smoother.apply(lvl.vec_cor, lvl.vec_def)))
          //return false;
        lvl.time_smooth += stamp_smooth.elapsed_now();

        // apply correction filter
        system_filter.filter_cor(lvl.vec_cor);

        // update solution vector
        lvl.vec_sol.axpy(lvl.vec_cor, lvl.vec_sol);

        // re-compute defect
        TimeStamp stamp_defect;
        system_matrix.apply(lvl.vec_def, lvl.vec_sol, lvl.vec_rhs, -DataType(1));
        lvl.time_defect += stamp_defect.elapsed_now();

        // apply defect filter
        system_filter.filter_def(lvl.vec_def);

        // okay
        return true;
      }

      /**
       * \brief Applies the peak-smoother on the current level.
       *
       * \param[in] cur_lvl
       * The current level.
       */
      Status _apply_smooth_peak(Index cur_lvl)
      {
        TimeStamp at;
        double mpi_exec_reduction_start(Statistics::get_time_mpi_execute_reduction());
        double mpi_exec_blas2_start(Statistics::get_time_mpi_execute_blas2());
        double mpi_exec_blas3_start(Statistics::get_time_mpi_execute_blas3());
        double mpi_exec_collective_start(Statistics::get_time_mpi_execute_collective());
        double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
        double mpi_wait_start_blas2(Statistics::get_time_mpi_wait_blas2());
        double mpi_wait_start_blas3(Statistics::get_time_mpi_wait_blas3());
        double mpi_wait_start_collective(Statistics::get_time_mpi_wait_collective());

        // get the level info
        LevelInfo& lvl = _hierarchy->_get_level_info(cur_lvl);

        // get the system matrix and filter
        const MatrixType& system_matrix = lvl.level->get_system_matrix();
        const FilterType& system_filter = lvl.level->get_system_filter();

        // compute defect
        TimeStamp stamp_defect;
        system_matrix.apply(lvl.vec_def, lvl.vec_sol, lvl.vec_rhs, -DataType(1));
        lvl.time_defect += stamp_defect.elapsed_now();

        // apply defect filter
        system_filter.filter_def(lvl.vec_def);

        // get the peak-smoother
        std::shared_ptr<SolverType> smoother_peak = lvl.level->get_smoother_peak();

        // do we have a peak-smoother?
        if(smoother_peak)
        {
          // apply peak-smoother
          if(!this->_apply_smooth_def(cur_lvl, *smoother_peak))
            return Status::aborted;
        }
        else
        {
          // no peak-smoother given
          // apply the pre- and post-smoothers instead
          std::shared_ptr<SolverType> smoother_pre = lvl.level->get_smoother_pre();
          std::shared_ptr<SolverType> smoother_post = lvl.level->get_smoother_post();

          if((smoother_pre)  && (!this->_apply_smooth_def(cur_lvl, *smoother_pre)))
            return Status::aborted;
          if((smoother_post) && (!this->_apply_smooth_def(cur_lvl, *smoother_post)))
            return Status::aborted;
        }

        TimeStamp bt;
        _toes.at(std::size_t(cur_lvl)) += bt.elapsed(at);
        double mpi_exec_reduction_stop(Statistics::get_time_mpi_execute_reduction());
        double mpi_exec_blas2_stop(Statistics::get_time_mpi_execute_blas2());
        double mpi_exec_blas3_stop(Statistics::get_time_mpi_execute_blas3());
        double mpi_exec_collective_stop(Statistics::get_time_mpi_execute_collective());
        double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
        double mpi_wait_stop_blas2(Statistics::get_time_mpi_wait_blas2());
        double mpi_wait_stop_blas3(Statistics::get_time_mpi_wait_blas3());
        double mpi_wait_stop_collective(Statistics::get_time_mpi_wait_collective());
        _mpi_execs_reduction.at(std::size_t(cur_lvl)) += mpi_exec_reduction_stop - mpi_exec_reduction_start;
        _mpi_execs_blas2.at(std::size_t(cur_lvl)) += mpi_exec_blas2_stop - mpi_exec_blas2_start;
        _mpi_execs_blas3.at(std::size_t(cur_lvl)) += mpi_exec_blas3_stop - mpi_exec_blas3_start;
        _mpi_execs_collective.at(std::size_t(cur_lvl)) += mpi_exec_collective_stop - mpi_exec_collective_start;
        _mpi_waits_reduction.at(std::size_t(cur_lvl)) += mpi_wait_stop_reduction - mpi_wait_start_reduction;
        _mpi_waits_blas2.at(std::size_t(cur_lvl)) += mpi_wait_stop_blas2 - mpi_wait_start_blas2;
        _mpi_waits_blas3.at(std::size_t(cur_lvl)) += mpi_wait_stop_blas3 - mpi_wait_start_blas3;
        _mpi_waits_collective.at(std::size_t(cur_lvl)) += mpi_wait_stop_collective - mpi_wait_start_collective;

        return Status::success;
      }

      /**
       * \brief Restricts from the current level onto the coarse level
       *
       * \param[in] cur_lvl
       * The level from which to restrict.
       *
       * \param[in] cur_smooth
       * Specifies whether to apply the pre-smoother on the current level.
       */
      Status _apply_rest(const Index cur_lvl, bool cur_smooth)
      {
        const Index last_level = Math::min(_crs_level, _hierarchy->size_physical());

        // restriction loop: from current level to last level
        for(Index i(cur_lvl); i < last_level; ++i)
        {
          TimeStamp at;
          double mpi_exec_reduction_start(Statistics::get_time_mpi_execute_reduction());
          double mpi_exec_blas2_start(Statistics::get_time_mpi_execute_blas2());
          double mpi_exec_blas3_start(Statistics::get_time_mpi_execute_blas3());
          double mpi_exec_collective_start(Statistics::get_time_mpi_execute_collective());
          double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_start_blas2(Statistics::get_time_mpi_wait_blas2());
          double mpi_wait_start_blas3(Statistics::get_time_mpi_wait_blas3());
          double mpi_wait_start_collective(Statistics::get_time_mpi_wait_collective());

          // get our fine and coarse levels
          LevelInfo& lvl_f = _hierarchy->_get_level_info(i);

          // get system matrix and filters
          const MatrixType& system_matrix   = lvl_f.level->get_system_matrix();
          const FilterType& system_filter_f = lvl_f.level->get_system_filter();

          // get our pre-smoother
          std::shared_ptr<SolverType> smoother = lvl_f.level->get_smoother_pre();

          // The following if-block ensures that we never apply the pre-smoother
          // when restricting from an intermediate peak-level in the F- or W-cycle,
          // as otherwise the following statements would discard the solution vector
          // approximation that has already been computed by prior restriction and
          // peak-smoothing operations.
          if(cur_smooth || (i > cur_lvl))
          {
            // apply pre-smoother if we have one
            if(smoother)
            {
              // apply pre-smoother
              TimeStamp stamp_smooth;
              Statistics::add_solver_expression(std::make_shared<ExpressionCallSmoother>(this->name(), smoother->name()));
              smoother->apply(lvl_f.vec_sol, lvl_f.vec_rhs);
              //if(!status_success(smoother->apply(lvl_f.vec_sol, lvl_f.vec_rhs)))
                //return Status::aborted;
              lvl_f.time_smooth += stamp_smooth.elapsed_now();

              // compute defect
              TimeStamp stamp_defect;
              system_matrix.apply(lvl_f.vec_def, lvl_f.vec_sol, lvl_f.vec_rhs, -DataType(1));
              lvl_f.time_defect += stamp_defect.elapsed_now();
            }
            else
            {
              // format solution
              lvl_f.vec_sol.format();

              // set our rhs as defect
              lvl_f.vec_def.copy(lvl_f.vec_rhs);
            }
          }

          // filter our defect
          system_filter_f.filter_def(lvl_f.vec_def);

          // get our transfer operator
          const TransferOperatorType* transfer_operator = lvl_f.level->get_transfer_operator();
          XASSERTM(transfer_operator != nullptr, "transfer operator is missing");

          // break loop?
          bool break_loop = false;

          // ghost transfer?
          if(transfer_operator->is_ghost())
          {
            // send restriction to parent processes and return
            TimeStamp stamp_rest;
            transfer_operator->rest_send(lvl_f.vec_def);
            lvl_f.time_transfer += stamp_rest.elapsed_now();
            break_loop = true;
          }
          else
          {
            // okay, get coarse level
            LevelInfo& lvl_c = _hierarchy->_get_level_info(i+1);
            const FilterType& system_filter_c = lvl_c.level->get_system_filter();

            // restrict onto coarse level
            //Statistics::add_solver_expression(std::make_shared<ExpressionRestriction>(this->name(), i));
            TimeStamp stamp_rest;
            transfer_operator->rest(lvl_f.vec_def, lvl_c.vec_rhs);
            lvl_f.time_transfer += stamp_rest.elapsed_now();

            // filter coarse defect
            system_filter_c.filter_def(lvl_c.vec_rhs);
          }

          // collect timings
          _toes.at((size_t)i) += at.elapsed_now();
          double mpi_exec_reduction_stop(Statistics::get_time_mpi_execute_reduction());
          double mpi_exec_blas2_stop(Statistics::get_time_mpi_execute_blas2());
          double mpi_exec_blas3_stop(Statistics::get_time_mpi_execute_blas3());
          double mpi_exec_collective_stop(Statistics::get_time_mpi_execute_collective());
          double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_stop_blas2(Statistics::get_time_mpi_wait_blas2());
          double mpi_wait_stop_blas3(Statistics::get_time_mpi_wait_blas3());
          double mpi_wait_stop_collective(Statistics::get_time_mpi_wait_collective());
          _mpi_execs_reduction.at((size_t)i) += mpi_exec_reduction_stop - mpi_exec_reduction_start;
          _mpi_execs_blas2.at((size_t)i) += mpi_exec_blas2_stop - mpi_exec_blas2_start;
          _mpi_execs_blas3.at((size_t)i) += mpi_exec_blas3_stop - mpi_exec_blas3_start;
          _mpi_execs_collective.at((size_t)i) += mpi_exec_collective_stop - mpi_exec_collective_start;
          _mpi_waits_reduction.at((size_t)i) += mpi_wait_stop_reduction - mpi_wait_start_reduction;
          _mpi_waits_blas2.at((size_t)i) += mpi_wait_stop_blas2 - mpi_wait_start_blas2;
          _mpi_waits_blas3.at((size_t)i) += mpi_wait_stop_blas3 - mpi_wait_start_blas3;
          _mpi_waits_collective.at((size_t)i) += mpi_wait_stop_collective - mpi_wait_start_collective;

          if(break_loop)
            break;

          // descent to prior level
        }

        return Status::success;
      }

      /**
       * \brief Prolongates from the coarse level onto the current level
       *
       * \param[in] cur_lvl
       * The level onto which to prolongate.
       *
       * \param[in] cur_smooth
       * Specifies whether to apply the post-smoother on the current level.
       */
      virtual Status _apply_prol(const Index cur_lvl, bool cur_smooth)
      {
        const Index last_level = Math::min(_crs_level, _hierarchy->size_physical());

        // prolongation loop: from last level to current level
        for(Index i(last_level); i > cur_lvl;)
        {
          --i;
          TimeStamp at;
          double mpi_exec_reduction_start(Statistics::get_time_mpi_execute_reduction());
          double mpi_exec_blas2_start(Statistics::get_time_mpi_execute_blas2());
          double mpi_exec_blas3_start(Statistics::get_time_mpi_execute_blas3());
          double mpi_exec_collective_start(Statistics::get_time_mpi_execute_collective());
          double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_start_blas2(Statistics::get_time_mpi_wait_blas2());
          double mpi_wait_start_blas3(Statistics::get_time_mpi_wait_blas3());
          double mpi_wait_start_collective(Statistics::get_time_mpi_wait_collective());

          // get our fine level
          LevelInfo& lvl_f = _hierarchy->_get_level_info(i);

          // get our transfer operator
          const TransferOperatorType* transfer_operator = lvl_f.level->get_transfer_operator();
          XASSERTM(transfer_operator != nullptr, "transfer operator is missing");

          // ghost transfer?
          if(transfer_operator->is_ghost())
          {
            // receive prolongation
            TimeStamp stamp_prol;
            transfer_operator->prol_recv(lvl_f.vec_cor);
            lvl_f.time_transfer += stamp_prol.elapsed_now();
          }
          else
          {
            // get coarse level
            LevelInfo& lvl_c = _hierarchy->_get_level_info(i+1);

            // prolongate
            TimeStamp stamp_prol;
            transfer_operator->prol(lvl_f.vec_cor, lvl_c.vec_sol);
            lvl_f.time_transfer += stamp_prol.elapsed_now();
          }

          // get system matrix and filters
          const MatrixType& system_matrix   = lvl_f.level->get_system_matrix();
          const FilterType& system_filter_f = lvl_f.level->get_system_filter();

          // apply correction filter
          system_filter_f.filter_cor(lvl_f.vec_cor);

          // initialize omega for (adaptive) coarse grid correction
          DataType omega_cgc = DataType(1);

          // use some sort of adaptive cgc?
          if(_adapt_cgc != MultiGridAdaptCGC::Fixed)
          {
            // multiply coarse grid correction by system matrix
            TimeStamp stamp_defect;
            system_matrix.apply(lvl_f.vec_tmp, lvl_f.vec_cor);
            lvl_f.time_defect += stamp_defect.elapsed_now();

            // apply filter
            system_filter_f.filter_def(lvl_f.vec_tmp);

            // compute adaptive omega
            switch(_adapt_cgc)
            {
            case MultiGridAdaptCGC::MinEnergy:
              omega_cgc = lvl_f.vec_def.dot(lvl_f.vec_cor)
                        / lvl_f.vec_tmp.dot(lvl_f.vec_cor);
              break;

            case MultiGridAdaptCGC::MinDefect:
              omega_cgc = lvl_f.vec_def.dot(lvl_f.vec_tmp)
                        / lvl_f.vec_tmp.dot(lvl_f.vec_tmp);
              break;

            default:
              // This default block is required to make the GCC happy...
              break;
            }
          }

          // update our solution vector
          lvl_f.vec_sol.axpy(lvl_f.vec_cor, lvl_f.vec_sol, omega_cgc);

          // get our post-smoother
          std::shared_ptr<SolverType> smoother = lvl_f.level->get_smoother_post();

          // apply post-smoother if we have one
          if(smoother && (cur_smooth || (i > cur_lvl)))
          {
            // compute new defect
            if(_adapt_cgc != MultiGridAdaptCGC::Fixed)
            {
              // We have used some sort of adaptive coarse grid correction.
              // In this case, the current solution vector is
              //    sol_new = sol_old + omega * cor
              //
              // Moreover, we have the following vectors at hand:
              //    1) tmp = A*cor
              //    2) def_old = rhs - A*sol_old
              //
              // So we can apply the same "trick" as in the PCG method to update our
              // defect vector instead of recomputing it by a matrix-vector mult, i.e.:
              //   def_new = rhs - A * sol_new
              //           = rhs - A * (sol_old + omega * cor)
              //           = rhs - A *  sol_old - omega * A * cor
              //           = def_old            - omega * tmp
              lvl_f.vec_def.axpy(lvl_f.vec_tmp, lvl_f.vec_def, -omega_cgc);
            }
            else
            {
              // we have to recompute the defect by (rhs - A*sol_new)
              TimeStamp stamp_defect;
              system_matrix.apply(lvl_f.vec_def, lvl_f.vec_sol, lvl_f.vec_rhs, -DataType(1));
              lvl_f.time_defect += stamp_defect.elapsed_now();

              // apply defect filter
              system_filter_f.filter_def(lvl_f.vec_def);
            }

            // apply post-smoother
            Statistics::add_solver_expression(std::make_shared<ExpressionCallSmoother>(this->name(), smoother->name()));
            TimeStamp stamp_smooth;
            smoother->apply(lvl_f.vec_cor, lvl_f.vec_def);
            //if(!status_success(smoother->apply(lvl_f.vec_cor, lvl_f.vec_def)))
              //return Status::aborted;
            lvl_f.time_smooth += stamp_smooth.elapsed_now();

            // update solution vector
            lvl_f.vec_sol.axpy(lvl_f.vec_cor, lvl_f.vec_sol);
          }

          _toes.at((size_t)i) += at.elapsed_now();
          double mpi_exec_reduction_stop(Statistics::get_time_mpi_execute_reduction());
          double mpi_exec_blas2_stop(Statistics::get_time_mpi_execute_blas2());
          double mpi_exec_blas3_stop(Statistics::get_time_mpi_execute_blas3());
          double mpi_exec_collective_stop(Statistics::get_time_mpi_execute_collective());
          double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_stop_blas2(Statistics::get_time_mpi_wait_blas2());
          double mpi_wait_stop_blas3(Statistics::get_time_mpi_wait_blas3());
          double mpi_wait_stop_collective(Statistics::get_time_mpi_wait_collective());
          _mpi_execs_reduction.at((size_t)i) += mpi_exec_reduction_stop - mpi_exec_reduction_start;
          _mpi_execs_blas2.at((size_t)i) += mpi_exec_blas2_stop - mpi_exec_blas2_start;
          _mpi_execs_blas3.at((size_t)i) += mpi_exec_blas3_stop - mpi_exec_blas3_start;
          _mpi_execs_collective.at((size_t)i) += mpi_exec_collective_stop - mpi_exec_collective_start;
          _mpi_waits_reduction.at((size_t)i) += mpi_wait_stop_reduction - mpi_wait_start_reduction;
          _mpi_waits_blas2.at((size_t)i) += mpi_wait_stop_blas2 - mpi_wait_start_blas2;
          _mpi_waits_blas3.at((size_t)i) += mpi_wait_stop_blas3 - mpi_wait_start_blas3;
          _mpi_waits_collective.at((size_t)i) += mpi_wait_stop_collective - mpi_wait_start_collective;

          // ascend to next level
        }

        return Status::success;
      }
    }; // class MultiGrid<...>

    /**
     * \brief Creates a new Multigrid preconditioner object
     *
     * \param[in] hierarchy
     * A pointer to the multigrid hierarchy object.
     *
     * \param[in] cycle
     * The desired multigrid cycle.
     *
     * \param[in] top_level
     * The desired top-level for this multigrid solver.\n
     * Set to 0 to use the finest level in the multigrid hierarchy.
     *
     * \param[in] crs_level
     * The desired coarse-level for this multigrid solver.\n
     * Set to -1 to use the coarsest level in the multigrid hierarchy.
     *
     * \returns
     * A shared pointer to a new MultiGrid object.
     */
    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename TransferOperator_>
    std::shared_ptr<MultiGrid<SystemMatrix_, SystemFilter_, TransferOperator_>> new_multigrid(
      std::shared_ptr<MultiGridHierarchy<SystemMatrix_, SystemFilter_, TransferOperator_>> hierarchy,
      MultiGridCycle cycle = MultiGridCycle::V,
      int top_level = 0,
      int crs_level = -1)
    {
      return std::make_shared<MultiGrid<SystemMatrix_, SystemFilter_, TransferOperator_>>
        (hierarchy, cycle, top_level, crs_level);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_MULTIGRID_HPP
