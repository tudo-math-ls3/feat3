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
      typename ProlOperator_,
      typename RestOperator_>
    class MultiGridLevelBase
    {
    public:
      /// the system matrix type
      typedef SystemMatrix_ SystemMatrixType;
      /// the system filter type
      typedef SystemFilter_ SystemFilterType;
      /// the prolongation operator type
      typedef ProlOperator_ ProlOperatorType;
      /// the restriction operator type
      typedef RestOperator_ RestOperatorType;
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
       * \brief Returns a const pointer to the prolongation operator.
       *
       * \note
       * This function may return \c nullptr on the coarse-grid level.
       */
      virtual const ProlOperatorType* get_prol_operator() const = 0;

      /**
       * \brief Returns a const pointer to the restriction operator.
       *
       * \note
       * This function may return \c nullptr on the coarse-grid level.
       */
      virtual const RestOperatorType* get_rest_operator() const = 0;

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

      /**
       * \brief Returns the scaling for adaptive coarse correction.
       *
       */
     virtual typename SystemMatrixType::DataType get_alpha_adaptive_coarse_correction() = 0;
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
      typename ProlOperator_,
      typename RestOperator_>
    class MultiGridLevelStd :
      public MultiGridLevelBase<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>
    {
    public:
      /// the base-class type
      typedef MultiGridLevelBase<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_> BaseClass;
      /// the system-matrix type
      typedef typename BaseClass::SystemMatrixType SystemMatrixType;
      /// the system-filter type
      typedef typename BaseClass::SystemFilterType SystemFilterType;
      /// the system-vector type
      typedef typename BaseClass::SystemVectorType SystemVectorType;
      /// the prolongation operator type
      typedef typename BaseClass::ProlOperatorType ProlOperatorType;
      /// the restriction operator type
      typedef typename BaseClass::RestOperatorType RestOperatorType;
      /// the coarse-grid solver/smoother type
      typedef typename BaseClass::SolverType SolverType;

    protected:
      /// the system matrix
      const SystemMatrixType& system_matrix;
      /// the system filter
      const SystemFilterType& system_filter;
      /// the prolongation matrix
      const ProlOperatorType* prol_operator;
      /// the restriction matrix
      const RestOperatorType* rest_operator;
      /// the coarse-grid solver
      std::shared_ptr<SolverType> coarse_solver;
      /// the pre-smoother
      std::shared_ptr<SolverType> smoother_pre;
      /// the post-smoother
      std::shared_ptr<SolverType> smoother_post;
      /// the peak-smoother
      std::shared_ptr<SolverType> smoother_peak;
      /// the adaptive coarse correction parameter
      typename BaseClass::SystemMatrixType::DataType alpha_adaptive_coarse_correction;

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
        prol_operator(nullptr),
        rest_operator(nullptr),
        coarse_solver(crs_solver),
        smoother_pre(),
        smoother_post(),
        smoother_peak(),
        alpha_adaptive_coarse_correction(typename BaseClass::SystemMatrixType::DataType(1))
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
       * \param[in] prol_operat
       * A reference to the prolongation operator.
       *
       * \param[in] rest_operat
       * A reference to the restriction operator.
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
        const ProlOperatorType& prol_operat,
        const RestOperatorType& rest_operat,
        std::shared_ptr<SolverType> smooth_pre,
        std::shared_ptr<SolverType> smooth_post,
        std::shared_ptr<SolverType> smooth_peak,
        std::shared_ptr<SolverType> crs_solver = nullptr,
        typename BaseClass::SystemMatrixType::DataType alpha = typename BaseClass::SystemMatrixType::DataType(1))
         :
        system_matrix(sys_matrix),
        system_filter(sys_filter),
        prol_operator(&prol_operat),
        rest_operator(&rest_operat),
        coarse_solver(crs_solver),
        smoother_pre(smooth_pre),
        smoother_post(smooth_post),
        smoother_peak(smooth_peak),
        alpha_adaptive_coarse_correction(alpha)
      {
      }

      /// virtal destructor
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

      /// \copydoc MultiGridLevelBase::get_prol_operator()
      virtual const ProlOperatorType* get_prol_operator() const override
      {
        return prol_operator;
      }

      /// \copydoc MultiGridLevelBase::get_rest_operator()
      virtual const RestOperatorType* get_rest_operator() const override
      {
        return rest_operator;
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

      /// \copydoc MultiGridLevelBase::get_alpha_adaptive_coarse_correction()
      virtual typename BaseClass::SystemMatrixType::DataType get_alpha_adaptive_coarse_correction() override
      {
        return alpha_adaptive_coarse_correction;
      }
    }; // class MultiGridLevelStd<...>

    // forward declaration
    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename ProlOperator_,
      typename RestOperator_>
    class MultiGrid;

    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename ProlOperator_,
      typename RestOperator_>
    class ScaRCMultiGrid;
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
      typename ProlOperator_,
      typename RestOperator_>
    class MultiGridHierarchy
    {
    public:
      /// MultiGrid is our friend
      friend class MultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>;
      friend class ScaRCMultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>;

      /// the level base class type
      typedef MultiGridLevelBase<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_> LevelType;
      /// the standard level class type
      typedef MultiGridLevelStd<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_> StdLevelType;

      /// the system matrix type
      typedef typename LevelType::SystemMatrixType SystemMatrixType;
      /// the system filter type
      typedef typename LevelType::SystemFilterType SystemFilterType;
      /// the system vector type
      typedef typename LevelType::SystemVectorType SystemVectorType;
      /// the prolongation operator type
      typedef typename LevelType::ProlOperatorType ProlOperatorType;
      /// the restriction operator tpye
      typedef typename LevelType::RestOperatorType RestOperatorType;
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
        SystemVectorType vec_rhs, vec_sol, vec_def, vec_cor;

        explicit LevelInfo(std::shared_ptr<LevelType> lvl) :
          level(lvl),
          unique_solvers()
        {
        }

        void init_solvers()
        {
          // build unique solver list
          unique_solvers.clear();
          _push_solver(level->get_coarse_solver().get());
          _push_solver(level->get_smoother_pre().get());
          _push_solver(level->get_smoother_post().get());
          _push_solver(level->get_smoother_peak().get());
        }

        void init_branch(const String parent)
        {
          // initialise all unique solvers in forward order
          for(auto it = unique_solvers.begin(); it != unique_solvers.end(); ++it)
            (*it)->init_branch(parent);
        }

        void init_symbolic()
        {
          // initialise all unique solvers in forward order
          for(auto it = unique_solvers.begin(); it != unique_solvers.end(); ++it)
            (*it)->init_symbolic();

          // get level system matrix
          const SystemMatrixType& matrix = level->get_system_matrix();

          // create our vectors
          vec_rhs = matrix.create_vector_r();
          vec_sol = matrix.create_vector_r();
          vec_def = matrix.create_vector_r();
          vec_cor = matrix.create_vector_r();
        }

        void done_symbolic()
        {
          // release our vectors
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
          // initialise all unique solvers in forward order
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
        XASSERTM(lvl < get_num_levels(), "invalid level index");
        return _levels.at(std::size_t(lvl));
      }

    protected:
      /// the deque of all level info objects
      std::deque<LevelInfo> _levels;
      /// specifies whether init_branch() has been called
      bool _have_init_branch;
      /// specifies whether init_symbolic() has been called
      bool _have_init_symbolic;
      /// specifies whether init_numeric() has been called
      bool _have_init_numeric;

    public:
      /**
       * \brief Constructor
       */
      MultiGridHierarchy() :
        _levels(),
        _have_init_branch(false),
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
        XASSERTM(!_have_init_branch,   "cannot push new level, init_branch() was already called");
        XASSERTM(!_have_init_symbolic, "cannot push new level, init_symbolic() was already called");
        XASSERTM(!_have_init_numeric,  "cannot push new level, init_numeric() was already called");
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
        XASSERTM(_levels.empty(), "cannot push coarse level as one already exists");
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
       * \param[in] prol_operator
       * A reference to the prolongation operator onto this level.
       *
       * \param[in] rest_operator
       * A reference to the restriction operator from this level.
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
        const ProlOperatorType& prol_operator,
        const RestOperatorType& rest_operator,
        std::shared_ptr<SolverType> smoother_pre,
        std::shared_ptr<SolverType> smoother_post,
        std::shared_ptr<SolverType> smoother_peak,
        std::shared_ptr<SolverType> coarse_solver = nullptr,
        typename SystemMatrixType::DataType alpha = typename SystemMatrixType::DataType(1)
        )
      {
        XASSERTM(!_levels.empty(), "cannot push level, no coarse level exists yet");
        push_level(std::make_shared<StdLevelType>(system_matrix, system_filter, prol_operator, rest_operator,
          smoother_pre, smoother_post, smoother_peak, coarse_solver, alpha));
      }

      /**
       * \brief Returns the number of levels in the hierarchy.
       */
      Index get_num_levels() const
      {
        return Index(_levels.size());
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
        XASSERTM(lvl < get_num_levels(), "invalid level index");
        return *(_levels.at(std::size_t(lvl)).level);
      }

      /** \copydoc get_level() */
      const LevelType& get_level(Index lvl) const
      {
        XASSERTM(lvl < get_num_levels(), "invalid level index");
        return *(_levels.at(std::size_t(lvl)).level);
      }

      /**
       * \brief Initialise solver branch description
       */
      void init_branch(String parent = "")
      {
        // loop over all levels in forward order
        for(auto it = _levels.begin(); it != _levels.end(); ++it)
        {
          it->init_solvers();
          it->init_branch(parent);
        }

        _have_init_branch = true;
      }

      /**
       * \brief Symbolic initialisation method
       *
       * This method is called to perform symbolic initialisation of the sub-solvers
       * in the multigrid hierarchy.
       */
      void init_symbolic()
      {
        XASSERTM(_have_init_branch, "init_branch() has not yet been called");
        XASSERTM(!_have_init_symbolic, "init_symbolic() was already called");
        XASSERTM(!_have_init_numeric, "init_numeric() was already called");

        // initialise all levels in forward order
        for(auto it = _levels.begin(); it != _levels.end(); ++it)
          it->init_symbolic();

        _have_init_symbolic = true;
      }

      /**
       * \brief Symbolic finalisation method
       *
       * This method is called to release any data allocated in the symbolic initialisation step
       * of the sub-solvers in the multigrid hierarchy.
       */
      void done_symbolic()
      {
        XASSERTM(_have_init_symbolic, "init_symbolic() has not yet been called");
        XASSERTM(!_have_init_numeric, "init_numeric() was already called");

        // release all levels in reverse order
        for(auto it = _levels.rbegin(); it != _levels.rend(); ++it)
          it->done_symbolic();

        _have_init_symbolic = false;
      }

      /**
       * \brief Numeric initialisation method
       *
       * This method is called to perform numeric initialisation of the sub-solvers
       * in the multigrid hierarchy.\n
       * Before this function can be called, the symbolic initialisation must be performed.
       */
      void init_numeric()
      {
        XASSERTM(_have_init_symbolic, "init_symbolic() has not yet been called");
        XASSERTM(!_have_init_numeric, "init_numeric() was already called");

        // initialise all levels in forward order
        for(auto it = _levels.begin(); it != _levels.end(); ++it)
          it->init_numeric();

        _have_init_numeric = true;
      }

      /**
       * \brief Numeric finalisation method
       *
       * This method is called to release any data allocated in the numeric initialisation step
       * of the sub-solvers in the multigrid hierarchy.
       */
      void done_numeric()
      {
        XASSERTM(_have_init_numeric, "init_numeric() has not yet been called");

        // release all levels in reverse order
        for(auto it = _levels.rbegin(); it != _levels.rend(); ++it)
          it->done_numeric();

        _have_init_numeric = false;
      }

      /**
       * \brief Initialisation method
       *
       * This function performs both the symbolic and numeric initialisation, i.e. it simply performs
         \verbatim
         this->init_branch();
         this->init_symbolic();
         this->init_numeric();
         \endverbatim
       */
      void init()
      {
        init_branch();
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
      void done()
      {
        done_numeric();
        done_symbolic();
      }

      /**
       * \brief Specifies whether the symbolic initialisation was performed.
       *
       * \returns
       * \c true between a call of init_symbolic() and done_symbolic(), otherwise \c false.
       */
      bool have_symbolic() const
      {
        return _have_init_symbolic;
      }

      /**
       * \brief Specifies whether the symbolic initialisation was performed.
       *
       * \returns
       * \c true between a call of init_symbolic() and done_symbolic(), otherwise \c false.
       */
      bool have_numeric() const
      {
        return _have_init_numeric;
      }
    }; // class MultiGridHierarchy<...>

    /**
     * \brief Multigrid preconditioner implementation
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
      typename ProlOperator_,
      typename RestOperator_>
    class MultiGrid :
      public SolverBase<typename SystemMatrix_::VectorTypeR>
    {
    public:
      /// our base-class
      typedef SolverBase<typename SystemMatrix_::VectorTypeR> BaseClass;

      /// our compatible multigrid hierarchy class
      typedef MultiGridHierarchy<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_> HierarchyType;

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
      /// the prolongation operator type
      typedef typename LevelType::ProlOperatorType ProlOperatorType;
      /// the restriction operator tpye
      typedef typename LevelType::RestOperatorType RestOperatorType;
      /// the sub-solver type
      typedef typename LevelType::SolverType SolverType;

      /// our data type
      typedef typename MatrixType::DataType DataType;

    protected:
      /// the multigrid hierarchy object
      std::shared_ptr<HierarchyType> _hierarchy;
      /// the multigrid cycle
      MultiGridCycle _cycle;
      /// the top-level of this multigrid
      Index _top_level;
      /// the coarse-level of this multigrid
      Index _crs_level;

      /// W-cycle level counter vector
      std::vector<int> _counters;

      /// array containing toe for each processed level
      std::vector<double> _toes;
      /// array containing toe of mpi execution for each processed level
      std::vector<double> _mpi_execs;
      /// array containing toe of mpi wait for each processed level
      std::vector<double> _mpi_waits;

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
       * Set to -1 to use the finest level in the multigrid hierarchy.
       *
       * \param[in] crs_level
       * The desired coarse-level for this multigrid solver.\n
       * Set to 0 to use the coarsest level in the multigrid hierarchy.
       *
       * \note
       * The cycle may be changed anytime by using the set_cycle() function.
       */
      explicit MultiGrid(std::shared_ptr<HierarchyType> hierarchy, MultiGridCycle cycle, int top_level = -1, int crs_level = 0) :
        BaseClass(),
        _hierarchy(hierarchy),
        _cycle(cycle),
        _top_level(top_level >= 0 ? Index(top_level) : hierarchy->get_num_levels()-1),
        _crs_level(Index(crs_level))
      {
        XASSERTM(crs_level >= 0, "invalid coarse level");
        XASSERTM(_top_level >= _crs_level, "invalid top-/coarse-level combination");
      }

      /// virtual destructor
      virtual ~MultiGrid()
      {
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
       * Set to -1 to use the finest level in the multigrid hierarchy.
       *
       * \param[in] crs_level
       * The desired coarse-level for this multigrid solver.\n
       * Set to 0 to use the coarsest level in the multigrid hierarchy.
       */
      void set_levels(int top_level, int crs_level)
      {
        XASSERTM(crs_level >= 0, "invalid coarse level");
        _top_level = (top_level >= 0 ? Index(top_level) : _hierarchy->get_num_levels()-1);
        _crs_level = Index(crs_level);
        XASSERTM(_top_level >= _crs_level, "invalid top-/coarse-level combination");
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
        return "MultiGrid";
      }

      /// \todo implement this correctly
      virtual void init_branch(String parent = "") override
      {
        BaseClass::init_branch(parent);
      }

      /**
       * \brief Symbolic initialisation function.
       *
       * \attention
       * This function does \b not perform the actual symbolic initialisation
       * of the multigrid hierarchy. You have to perform the symbolic initialisation
       * explicitly by calling the init_symbolic function of the corresponding
       * MultiGridHierarchy object \b before initialising the solver tree!
       */
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // ensure that the hierarchy was initialised
        /// \todo does this work in nested multigrid?
        XASSERTM(_hierarchy->have_symbolic(), "init_symbolic() of multigrid hierarchy was not called yet");

        const std::size_t num_lvls = _hierarchy->get_num_levels();

        // allocate counter vector
        _counters.resize(std::size_t(num_lvls), 0);

        // allocate statistics vectors
        _toes.resize(std::size_t(num_lvls), 0.0);
        _mpi_execs.resize(std::size_t(num_lvls), 0.0);
        _mpi_waits.resize(std::size_t(num_lvls), 0.0);
      }

      /**
       * \brief Symbolic finalisation function.
       *
       * \attention
       * This function does \b not perform the actual symbolic finalisation
       * of the multigrid hierarchy. You have to perform the symbolic finalisation
       * explicitly by calling the done_symbolic function of the corresponding
       * MultiGridHierarchy object \b after finalising the solver tree!
       */
      virtual void done_symbolic() override
      {
        _mpi_waits.clear();
        _mpi_execs.clear();
        _toes.clear();
        _counters.clear();

        BaseClass::done_symbolic();
      }

      /**
       * \brief Numeric initialisation function.
       *
       * \attention
       * This function does \b not perform the actual numeric initialisation
       * of the multigrid hierarchy. You have to perform the numeric initialisation
       * explicitly by calling the init_numeric function of the corresponding
       * MultiGridHierarchy object \b before initialising the solver tree!
       */
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // ensure that the hierarchy was initialised
        /// \todo does this work in nested multigrid?
        XASSERTM(_hierarchy->have_numeric(), "init_numeric() of multigrid hierarchy was not called yet");
      }

      /**
       * \brief Numeric finalisation function.
       *
       * \attention
       * This function does \b not perform the actual numeric finalisation
       * of the multigrid hierarchy. You have to perform the numeric finalisation
       * explicitly by calling the done_numeric function of the corresponding
       * MultiGridHierarchy object \b after finalising the solver tree!
       */
      virtual void done_numeric() override
      {
        BaseClass::done_numeric();
      }

      /// \todo be more verbosive on different pre/post/peak smoother infos
      virtual String get_formatted_solver_tree() const override
      {
        String result;
        result += this->name() + "-" + stringify(_cycle);
        result += " ( ";
        std::shared_ptr<SolverType> smoother = _hierarchy->_get_level_info(_top_level).level->get_smoother_pre();
        if(smoother)
        {
          result += "S: " + smoother->get_formatted_solver_tree() + " / ";
        }
        else
        {
          result += "S: None / ";
        }

        std::shared_ptr<SolverType> coarse_solver = _hierarchy->_get_level_info(_crs_level).level->get_coarse_solver();
        if(coarse_solver)
        {
          result += "C: " + coarse_solver->get_formatted_solver_tree();

        }
        else
        {
          result += "C: Identity";
        }
        result += " ) ";

        return result;
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
        // insert -1 to signal new starting cycle
        Statistics::add_solver_toe(this->_branch, double(-1));
        Statistics::add_solver_mpi_execute(this->_branch, double(-1));
        Statistics::add_solver_mpi_wait(this->_branch, double(-1));

        // reset statistics counters
        for(std::size_t i(0); i <= std::size_t(_top_level); ++i)
        {
          _toes[i] = 0.0;
          _mpi_execs[i] = 0.0;
          _mpi_waits[i] = 0.0;
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
          status = Status::aborted;
          break;
        }

        // copy solution to correction vector
        vec_cor.copy(lvl_top.vec_sol);

        // propagate solver statistics
        for(std::size_t i(0); i <= std::size_t(_top_level); ++i)
        {
          Statistics::add_solver_toe(this->_branch, _toes.at(i));
          Statistics::add_solver_mpi_execute(this->_branch, _mpi_execs.at(i));
          Statistics::add_solver_mpi_wait(this->_branch, _mpi_waits.at(i));
        }

        // okay
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

        // restrict (apply top-level pre-smoother if if it exists)
        this->_apply_rest(_top_level, true);

        // F-cycle intermediate peak-level loop
        for(Index peak_lvl(_crs_level+1); peak_lvl < _top_level; ++peak_lvl)
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

        // solve coarse level
        this->_apply_coarse();

        // prolongate (apply top-level post-smoother if if it exists)
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

        // restrict (apply top-level pre-smoother if if it exists)
        this->_apply_rest(_top_level, true);

        // solve coarse level
        this->_apply_coarse();

        // reset all level counters to 0
        for(std::size_t i(_crs_level); i <= std::size_t(_top_level); ++i)
        {
          _counters[i] = 0;
        }

        // compute the total number of coarse-grid solves to perform:
        // that's = 2^(top_level - crs_level)
        const int num_cgs = 1 << (_top_level - _crs_level);

        // W-cycle loop: count the number of coarse-grid solves
        for(int cgs(1); cgs < num_cgs; ++cgs)
        {
          // find the next peak level; that's the lowest level
          // above the coarse level whose counter is 0
          std::size_t peak_lvl = std::size_t(_crs_level+1);
          for(; peak_lvl <= std::size_t(_top_level); ++peak_lvl)
          {
            // if the counter is 0, we have our next peak level
            if(_counters[peak_lvl] == 0)
              break;
          }

          // >>>>> DEBUG >>>>>
          //std::cout << stringify(cgs).pad_front(3) << " ";
          //std::cout << stringify(peak_lvl).pad_front(2) << ":";
          //for(std::size_t i(0); i <= std::size_t(top_lvl); ++i)
          //  std::cout << " " << _counters[i];
          // <<<<< DEBUG <<<<<

          // sanity check: do not go beyond our top level
          XASSERTM(peak_lvl <= _top_level, "W-cycle sanity check failed: invalid peak level");

          // reset all counters below our peak level to 0
          for(std::size_t i = std::size_t(_crs_level+1); i < peak_lvl; ++i)
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
          this->_apply_prol(peak_lvl, false);

          // apply peak-smoother
          this->_apply_smooth_peak(peak_lvl);

          // restrict from current peak-level without pre-smoothing
          this->_apply_rest(peak_lvl, false);

          // solve coarse level
          this->_apply_coarse();
        }

        // sanity check: all level counters above the coarse level must be 1
        for(std::size_t i(_crs_level+1); i <= std::size_t(_top_level); ++i)
        {
          XASSERTM(this->_counters[i] == 1, "W-cycle sanity check failed: invalid counters");
        }

        // prolongate (apply top-level post-smoother if if it exists)
        this->_apply_prol(_top_level, true);

        // okay
        return Status::success;
      }

      /**
       * \brief Applies the coarse grid solver on the coarse level.
       */
      Status _apply_coarse()
      {
        // process the coarse grid level
        TimeStamp at;
        double mpi_exec_start(Statistics::get_time_mpi_execute());
        double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
        double mpi_wait_start_spmv(Statistics::get_time_mpi_wait_spmv());

        // get the coarse level info
        LevelInfo& lvl_crs = _hierarchy->_get_level_info(_crs_level);

        // get the system filter
        const FilterType& system_filter = lvl_crs.level->get_system_filter();

        // get the coarse grid solver
        std::shared_ptr<SolverType> coarse_solver = lvl_crs.level->get_coarse_solver();

        // if the have a coarse grid solver, apply it
        if(coarse_solver)
        {
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

        TimeStamp bt;
        _toes.at(std::size_t(_crs_level)) = bt.elapsed(at);
        double mpi_exec_stop(Statistics::get_time_mpi_execute());
        double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
        double mpi_wait_stop_spmv(Statistics::get_time_mpi_wait_spmv());
        _mpi_execs.at(std::size_t(_crs_level)) = mpi_exec_stop - mpi_exec_start;
        _mpi_waits.at(std::size_t(_crs_level)) = (mpi_wait_stop_reduction - mpi_wait_start_reduction) + (mpi_wait_stop_spmv - mpi_wait_start_spmv);

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
        if(!status_success(smoother.apply(lvl.vec_cor, lvl.vec_def)))
          return false;

        // apply correction filter
        system_filter.filter_cor(lvl.vec_cor);

        // update solution vector
        lvl.vec_sol.axpy(lvl.vec_cor, lvl.vec_sol);

        // re-compute defect
        system_matrix.apply(lvl.vec_def, lvl.vec_sol, lvl.vec_rhs, -DataType(1));

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
        double mpi_exec_start(Statistics::get_time_mpi_execute());
        double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
        double mpi_wait_start_spmv(Statistics::get_time_mpi_wait_spmv());

        // get the level info
        LevelInfo& lvl = _hierarchy->_get_level_info(cur_lvl);

        // get the system matrix and filter
        const MatrixType& system_matrix = lvl.level->get_system_matrix();
        const FilterType& system_filter = lvl.level->get_system_filter();

        // compute defect
        system_matrix.apply(lvl.vec_def, lvl.vec_sol, lvl.vec_rhs, -DataType(1));

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
        _toes.at(std::size_t(cur_lvl))= bt.elapsed(at);
        double mpi_exec_stop(Statistics::get_time_mpi_execute());
        double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
        double mpi_wait_stop_spmv(Statistics::get_time_mpi_wait_spmv());
        _mpi_execs.at(std::size_t(cur_lvl)) = mpi_exec_stop - mpi_exec_start;
        _mpi_waits.at(std::size_t(cur_lvl)) = (mpi_wait_stop_reduction - mpi_wait_start_reduction) + (mpi_wait_stop_spmv - mpi_wait_start_spmv);

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
        // restriction loop: from current level to corse level
        for(Index i(cur_lvl); i > _crs_level; --i)
        {
          TimeStamp at;
          double mpi_exec_start(Statistics::get_time_mpi_execute());
          double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_start_spmv(Statistics::get_time_mpi_wait_spmv());

          // get our fine and coarse levels
          LevelInfo& lvl_f = _hierarchy->_get_level_info(i);
          LevelInfo& lvl_c = _hierarchy->_get_level_info(i-1);

          // get system matrix and filters
          const MatrixType& system_matrix   = lvl_f.level->get_system_matrix();
          const FilterType& system_filter_f = lvl_f.level->get_system_filter();
          const FilterType& system_filter_c = lvl_c.level->get_system_filter();

          // get our restriction operator
          const RestOperatorType* rest_operator = lvl_f.level->get_rest_operator();
          XASSERTM(rest_operator != nullptr, "restriction operator is missing");

          // get our pre-smoother
          std::shared_ptr<SolverType> smoother = lvl_f.level->get_smoother_pre();

          // The following if-block ensures that we never apply the pre-smoother
          // when restricting from an intermediate peak-level in the F- or W-cycle,
          // as otherwise the following statements would discard the solution vector
          // approximation that has already been computed by prior restriction and
          // peak-smoothing operations.
          if(cur_smooth || (i < cur_lvl))
          {
            // apply pre-smoother if we have one
            if(smoother)
            {
              // apply pre-smoother
              if(!status_success(smoother->apply(lvl_f.vec_sol, lvl_f.vec_rhs)))
                return Status::aborted;

              // compute defect
              system_matrix.apply(lvl_f.vec_def, lvl_f.vec_sol, lvl_f.vec_rhs, -DataType(1));
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

          // restrict onto coarse level
          rest_operator->apply(lvl_c.vec_rhs, lvl_f.vec_def);

          // filter coarse fefect
          system_filter_c.filter_def(lvl_c.vec_rhs);

          TimeStamp bt;
          _toes.at((size_t)i)= bt.elapsed(at);
          double mpi_exec_stop(Statistics::get_time_mpi_execute());
          double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_stop_spmv(Statistics::get_time_mpi_wait_spmv());
          _mpi_execs.at((size_t)i) = mpi_exec_stop - mpi_exec_start;
          _mpi_waits.at((size_t)i) = (mpi_wait_stop_reduction - mpi_wait_start_reduction) + (mpi_wait_stop_spmv - mpi_wait_start_spmv);

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
        // prolongation loop: from coarse level to current level
        for(Index i(_crs_level+1); i <= cur_lvl; ++i)
        {
          TimeStamp at;
          double mpi_exec_start(Statistics::get_time_mpi_execute());
          double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_start_spmv(Statistics::get_time_mpi_wait_spmv());

          // get our level and the coarse level
          LevelInfo& lvl_f = _hierarchy->_get_level_info(i);
          LevelInfo& lvl_c = _hierarchy->_get_level_info(i-1);

          // get system matrix and filters
          const MatrixType& system_matrix   = lvl_f.level->get_system_matrix();
          const FilterType& system_filter_f = lvl_f.level->get_system_filter();

          // get our prolongation operator
          const ProlOperatorType* prol_operator = lvl_f.level->get_prol_operator();
          XASSERTM(prol_operator != nullptr, "prolongation operator is missing");

          // prolongate coarse grid solution
          prol_operator->apply(lvl_f.vec_cor, lvl_c.vec_sol);

          // apply correction filter
          system_filter_f.filter_cor(lvl_f.vec_cor);

          // update our solution vector
          lvl_f.vec_sol.axpy(lvl_f.vec_cor, lvl_f.vec_sol, lvl_f.level->get_alpha_adaptive_coarse_correction());

          // get our post-smoother
          std::shared_ptr<SolverType> smoother = lvl_f.level->get_smoother_post();

          // apply post-smoother if we have one
          if(smoother && (cur_smooth || (i < cur_lvl)))
          {
            // compute new defect
            system_matrix.apply(lvl_f.vec_def, lvl_f.vec_sol, lvl_f.vec_rhs, -DataType(1));

            // apply defect filter
            system_filter_f.filter_def(lvl_f.vec_def);

            // apply post-smoother
            if(!status_success(smoother->apply(lvl_f.vec_cor, lvl_f.vec_def)))
              return Status::aborted;

            // update solution vector
            lvl_f.vec_sol.axpy(lvl_f.vec_cor, lvl_f.vec_sol);
          }

          TimeStamp bt;
          _toes.at((size_t)i) += bt.elapsed(at);
          double mpi_exec_stop(Statistics::get_time_mpi_execute());
          double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_stop_spmv(Statistics::get_time_mpi_wait_spmv());
          _mpi_execs.at((size_t)i) += mpi_exec_stop - mpi_exec_start;
          _mpi_waits.at((size_t)i) += (mpi_wait_stop_reduction - mpi_wait_start_reduction) + (mpi_wait_stop_spmv - mpi_wait_start_spmv);

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
     * Set to -1 to use the finest level in the multigrid hierarchy.
     *
     * \param[in] crs_level
     * The desired coarse-level for this multigrid solver.\n
     * Set to 0 to use the coarsest level in the multigrid hierarchy.
     *
     * \returns
     * A shared pointer to a new MultiGrid object.
     */
    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename ProlOperator_,
      typename RestOperator_>
    std::shared_ptr<MultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>> new_multigrid(
      std::shared_ptr<MultiGridHierarchy<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>> hierarchy,
      MultiGridCycle cycle = MultiGridCycle::V,
      int top_level = -1,
      int crs_level = 0)
    {
      return std::make_shared<MultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>>
        (hierarchy, cycle, top_level, crs_level);
    }

    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename ProlOperator_,
      typename RestOperator_>
    class ScaRCMultiGrid :
      public MultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>
    {
      ///Use MultiGrid's CTORs
      using MultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>::MultiGrid;

      public:
      typedef MultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_> BaseClass;

      /**
       * \brief Prolongates from the coarse level onto the current level
       *
       * \param[in] cur_lvl
       * The level onto which to prolongate.
       *
       * \param[in] cur_smooth
       * Specifies whether to apply the post-smoother on the current level.
       */
      virtual Status _apply_prol(const Index cur_lvl, bool cur_smooth) override
      {
        // prolongation loop: from coarse level to current level
        for(Index i(this->_crs_level+1); i <= cur_lvl; ++i)
        {
          TimeStamp at;
          double mpi_exec_start(Statistics::get_time_mpi_execute());
          double mpi_wait_start_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_start_spmv(Statistics::get_time_mpi_wait_spmv());

          // get our level and the coarse level
          typename BaseClass::LevelInfo& lvl_f = this->_hierarchy->_get_level_info(i);
          typename BaseClass::LevelInfo& lvl_c = this->_hierarchy->_get_level_info(i-1);

          // get system matrix and filters
          const typename BaseClass::MatrixType& system_matrix   = lvl_f.level->get_system_matrix();
          const typename BaseClass::FilterType& system_filter_f = lvl_f.level->get_system_filter();

          // get our prolongation operator
          const typename BaseClass::ProlOperatorType* prol_operator = lvl_f.level->get_prol_operator();
          XASSERTM(prol_operator != nullptr, "prolongation operator is missing");

          // prolongate coarse grid solution
          prol_operator->apply(lvl_f.vec_cor, lvl_c.vec_sol);

          // apply correction filter
          system_filter_f.filter_cor(lvl_f.vec_cor);

          // update our solution vector
          auto alpha(lvl_f.level->get_alpha_adaptive_coarse_correction());
          decltype(lvl_f.vec_cor) tmp(lvl_f.vec_cor.size());
          lvl_f.level->get_system_matrix().apply(tmp, lvl_f.vec_cor);
          alpha = lvl_f.vec_cor.dot(lvl_f.vec_def) / lvl_f.vec_cor.dot(tmp);
          lvl_f.vec_sol.axpy(lvl_f.vec_cor, lvl_f.vec_sol, alpha);

          // get our post-smoother
          std::shared_ptr<typename BaseClass::SolverType> smoother = lvl_f.level->get_smoother_post();

          // apply post-smoother if we have one
          if(smoother && (cur_smooth || (i < cur_lvl)))
          {
            // compute new defect
            system_matrix.apply(lvl_f.vec_def, lvl_f.vec_sol, lvl_f.vec_rhs, -typename BaseClass::DataType(1));

            // apply defect filter
            system_filter_f.filter_def(lvl_f.vec_def);

            // apply post-smoother
            if(!status_success(smoother->apply(lvl_f.vec_cor, lvl_f.vec_def)))
              return Status::aborted;

            // update solution vector
            lvl_f.vec_sol.axpy(lvl_f.vec_cor, lvl_f.vec_sol);
          }

          TimeStamp bt;
          this->_toes.at((size_t)i) += bt.elapsed(at);
          double mpi_exec_stop(Statistics::get_time_mpi_execute());
          double mpi_wait_stop_reduction(Statistics::get_time_mpi_wait_reduction());
          double mpi_wait_stop_spmv(Statistics::get_time_mpi_wait_spmv());
          this->_mpi_execs.at((size_t)i) += mpi_exec_stop - mpi_exec_start;
          this->_mpi_waits.at((size_t)i) += (mpi_wait_stop_reduction - mpi_wait_start_reduction) + (mpi_wait_stop_spmv - mpi_wait_start_spmv);

          // ascend to next level
        }

        return Status::success;
      }
    };

    template<
      typename SystemMatrix_,
      typename SystemFilter_,
      typename ProlOperator_,
      typename RestOperator_>
    std::shared_ptr<MultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>> new_scarcmultigrid(
      std::shared_ptr<MultiGridHierarchy<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>> hierarchy,
      MultiGridCycle cycle = MultiGridCycle::V,
      int top_level = -1,
      int crs_level = 0)
    {
      return std::make_shared<ScaRCMultiGrid<SystemMatrix_, SystemFilter_, ProlOperator_, RestOperator_>>
        (hierarchy, cycle, top_level, crs_level);
    }

  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_MULTIGRID_HPP
