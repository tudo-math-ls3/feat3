// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/solver/adp_solver_base.hpp>

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    /**
     * \brief DirectSparseSolver backend wrapper namespace
     *
     * This namespace encapsulates the various backend core wrapper functions used by the
     * DirectSparseSolver class. There is no need for you to look any further into this namespace
     * unless you want to add another backend to the DirectSparseSolver class.
     */
    namespace DSS
    {
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // NVIDIA CUDSS
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef FEAT_HAVE_CUDSS
      // cuDSS is available as backend
      static constexpr bool have_cudss = true;

      /// native data type of cuDSS solver
      typedef double CUDSS_DT;

      /// native index type of cuDSS solver
      typedef std::int32_t CUDSS_IT;

      void* create_cudss_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes);

      void destroy_cudss_core(void* core);

      CUDSS_IT* get_cudss_row_ptr(void* core);
      CUDSS_IT* get_cudss_col_idx(void* core);
      CUDSS_DT* get_cudss_mat_val(void* core);
      CUDSS_DT* get_cudss_rhs_val(void* core);
      CUDSS_DT* get_cudss_sol_val(void* core);

      void init_cudss_symbolic(void* core);
      void init_cudss_numeric(void* core);

      void solve_cudss(void* core);

      std::int64_t get_peak_mem_cudss_host(void* core);
      std::int64_t get_peak_mem_cudss_device(void* core);
#else
      // cuDSS is not available as backend
      static constexpr bool have_cudss = false;
#endif // FEAT_HAVE_CUDSS

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // INTEL MKL-DSS
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef FEAT_HAVE_MKL
      // MKL-DSS is available as backend
      static constexpr bool have_mkldss = true;

      /// native data type for MKL CSS solver
      typedef double MKLDSS_DT;

      /// native index type for MKL DSS solver
      /// must be identical to 'MKL_INT'
#ifdef MKL_ILP64
      typedef long long MKLDSS_IT;
#else
      typedef int MKLDSS_IT;
#endif

      void* create_mkldss_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes);

      void destroy_mkldss_core(void* core);

      MKLDSS_IT* get_mkldss_row_ptr(void* core);
      MKLDSS_IT* get_mkldss_col_idx(void* core);
      MKLDSS_DT* get_mkldss_mat_val(void* core);
      MKLDSS_DT* get_mkldss_rhs_val(void* core);
      MKLDSS_DT* get_mkldss_sol_val(void* core);

      void init_mkldss_symbolic(void* core);
      void init_mkldss_numeric(void* core);

      void solve_mkldss(void* core);

      std::int64_t get_peak_mem_mkldss(void* core);
#else
      // MKL-DSS is not available as backend
      static constexpr bool have_mkldss = false;
#endif // FEAT_HAVE_MKL

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // MUMPS
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef FEAT_HAVE_MUMPS
      // SuperLU is available as backend
      static constexpr bool have_mumps = true;

      /// native data type of MUMPS solver
      typedef double MUMPS_DT;

      /// index type of MUMPS solver; this is independent of 'MUMPS_INT', because we have to convert the
      /// 0-based CSR arrays to 1-based COO arrays anyways
      typedef std::int64_t MUMPS_IT;

      void* create_mumps_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes);

      void destroy_mumps_core(void* core);

      MUMPS_IT* get_mumps_row_ptr(void* core);
      MUMPS_IT* get_mumps_col_idx(void* core);
      MUMPS_DT* get_mumps_mat_val(void* core);
      MUMPS_DT* get_mumps_vector(void* core);

      void init_mumps_symbolic(void* core);
      void init_mumps_numeric(void* core);
      void done_mumps_numeric(void* core);

      void solve_mumps(void* core);
#else
      // MUMPS is not available as backend
      static constexpr bool have_mumps = false;
#endif // FEAT_HAVE_MUMPS_DIST

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // SUPERLU_DIST
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef FEAT_HAVE_SUPERLU_DIST
      // SuperLU is available as backend
      static constexpr bool have_superlu = true;

      /// native data type of SuperLU solver
      typedef double SUPERLU_DT;

      /// native index type of SuperLU solver
      /// must be identical to 'int_t' in <superlu_defs.h>
#ifdef FEAT_TPL_SUPERLU_INT64
      typedef std::int64_t SUPERLU_IT;
#else
      typedef int SUPERLU_IT;
#endif

      void* create_superlu_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes);

      void destroy_superlu_core(void* core);

      SUPERLU_IT* get_superlu_row_ptr(void* core);
      SUPERLU_IT* get_superlu_col_idx(void* core);
      SUPERLU_DT* get_superlu_mat_val(void* core);
      SUPERLU_DT* get_superlu_vector(void* core);

      void init_superlu_symbolic(void* core);
      void init_superlu_numeric(void* core);
      void done_superlu_numeric(void* core);

      void solve_superlu(void* core);
#else
      // SuperLU is not available as backend
      static constexpr bool have_superlu = false;
#endif // FEAT_HAVE_SUPERLU_DIST

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // UMFPACK (SuiteSparse)
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef FEAT_HAVE_UMFPACK
      // UMFPACK is available as backend
      static constexpr bool have_umfpack = true;

      /// native data type of UMFPACK solver
      typedef double UMFPACK_DT;

      /// native index type of UMFPACK solver
      /// note: UMFPACK always provides 32 and 64 bit interfaces, so we can choose which one we want to use here
      typedef std::int32_t UMFPACK_IT;

      void* create_umfpack_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes);

      void destroy_umfpack_core(void* core);

      UMFPACK_IT* get_umfpack_row_ptr(void* core);
      UMFPACK_IT* get_umfpack_col_idx(void* core);
      UMFPACK_DT* get_umfpack_mat_val(void* core);
      UMFPACK_DT* get_umfpack_rhs_val(void* core);
      UMFPACK_DT* get_umfpack_sol_val(void* core);

      void init_umfpack_symbolic(void* core);
      void init_umfpack_numeric(void* core);
      void done_umfpack_numeric(void* core);

      void solve_umfpack(void* core);
#else
      // UMFPACK is not available as backend
      static constexpr bool have_umfpack = false;
#endif // FEAT_HAVE_UMFPACK

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      /// specifies whether at least one backend for local systems is available
      static constexpr bool have_backend_local = have_cudss || have_mkldss || have_superlu || have_umfpack;

      /// specifies whether at least one backend for global systems is available
      static constexpr bool have_backend_global = have_cudss || have_mkldss || have_superlu;

    } // namespace DSS
    /// \endcond

    /**
     * \brief Enumeration of allowed backends for DirectSparseSolver class
     */
    enum class DSSBackend : int
    {
      none     = 0x0000, ///< no backend allowed
      cudss    = 0x0001, ///< cuDSS backend
      mumps    = 0x0002, ///< MUMPS backend
      mkldss   = 0x0004, ///< MKL-DSS backend
      superlu  = 0x0008, ///< SuperLU backend
      umfpack  = 0x0010, ///< UMFPACK backend
      all      = 0x001F  ///< all backends allowed
    };

    /// bit-wise OR operator for DSSBackend enum values
    static constexpr inline DSSBackend operator|(DSSBackend a, DSSBackend b)
    {
      return static_cast<DSSBackend>(static_cast<int>(a) | static_cast<int>(b));
    }

    /// bit-wise AND operator for DSSBackend enum values
    static constexpr inline DSSBackend operator&(DSSBackend a, DSSBackend b)
    {
      return static_cast<DSSBackend>(static_cast<int>(a) & static_cast<int>(b));
    }

    /// checks whether at least one of the bits in DSSBackend is selected
    static constexpr inline bool operator*(DSSBackend a)
    {
      // just check the interesting bits
      return (static_cast<int>(a) & static_cast<int>(DSSBackend::all)) != 0;
    }

    /// output operator for DSSBackend
    static std::ostream& operator<<(std::ostream& os, DSSBackend a)
    {
      int n = 0;
      if(*(a & DSSBackend::cudss))
      {
        if(++n > 1)
          os << '|';
        os << "cudss";
      }
      if(*(a & DSSBackend::mkldss))
      {
        if(++n > 1)
          os << '|';
        os << "mkldss";
      }
      if(*(a & DSSBackend::mumps))
      {
        if(++n > 1)
          os << '|';
        os << "mumps";
      }
      if(*(a & DSSBackend::superlu))
      {
        if(++n > 1)
          os << '|';
        os << "superlu";
      }
      if(*(a & DSSBackend::umfpack))
      {
        if(++n > 1)
          os << '|';
        os << "umfpack";
      }
      if(n <= 0)
        os << "none";
      return os;
    }

    //enum class DSSSystemHint
    //{
    //  none             = 0x0000,
    //  sym_struct       = 0x0001, ///< matrix has symmetric structure
    //  sym_values       = 0x0002, ///< matrix is symmetric
    //  pos_definite     = 0x0004  ///< matrix is positive definite
    //};

    /**
     * \brief Exception base class for errors occurring in DirectSparseSolver
     */
    class DirectSparseSolverException :
      public SolverException
    {
    public:
      explicit DirectSparseSolverException(const String& msg) :
        SolverException("DirectSparseSolver: " + msg)
      {
      }

      explicit DirectSparseSolverException(const String& backend, const String& msg) :
        SolverException("DirectSparseSolver [" + backend + "]: " + msg)
      {
      }
    }; // class DirectSparseSolverException

    /**
     * \brief DirectSparseSolver backend not found exception
     *
     * This exception is thrown if a requested backend is not available or if no suitable backend
     * was found at all.
     */
    class DSSBackendNotFoundException :
      public DirectSparseSolverException
    {
    public:
      DSSBackendNotFoundException() :
        DirectSparseSolverException("DirectSparseSolver: no suitable backend found")
      {
      }

      explicit DSSBackendNotFoundException(DSSBackend b) :
        DirectSparseSolverException("DirectSparseSolver: Backend not available: " + stringify(b))
      {
      }
    }; // class DSSBackendNotFoundException

    /**
     * \brief Front-end wrapper class for (parallel) third-party direct sparse solvers
     *
     * This class acts as a unified front-end that enables the use of various direct sparse solvers
     * implemented in third-party libraries. This class derives from the ADPSolverBase class and
     * can therefore be applied to any type of system that is supported by the ADPSolverBase class.
     *
     * In particular, this class currently offers the following backends:
     * - cuDSS: A MPI-parallel direct solver for CUDA, which utilizes the GPU to solve the system.
     *   Obviously requires CUDA and the cuDSS library and can only utilize a single GPU per MPI
     *   process.
     * - MKL-DSS: A MPI-parallel direct solver, which is part of the Intel oneAPI MKL libraries,
     *   which is typically the best choice if the solver is applied on multiple processes.
     * - MUMPS: A MPI-parallel direct solver written in Fortran.
     * - SuperLU: A MPI-parallel direct solver implemented in the SuperLU_Dist library. This solver
     *   is usually slower that UMFPACK in the single-process case and slower that MKL-DSS in the
     *   multi-process case, so it serves as a fallback if the other ones are not available or not
     *   suitable.
     * - UMFPACK: A sequential direct solver, which is typically the best choice if the solver is
     *   applied only on a single MPI process or in sequential builds.
     *
     * \attention
     * If you intend to use the cuDSS library in combination with MPI in "multi-node-multi-gpu"
     * mode, then you need to ensure that the <c>CUDSS_COMM_LIB</c> environment variable is set
     * to the absolute path of cuDSS's <c>libcudss_commlayer_openmpi.so</c> library, because the
     * code will otherwise crash with obscure illegal memory access errors when trying to run the
     * solver with more than 1 MPI process!
     * See https://docs.nvidia.com/cuda/cudss/advanced_features.html#multi-gpu-multi-node-mgmn-mode
     * for more details.
     *
     * \note
     * Technically, the MKL-DSS backend includes two different backends, namely a backend for the
     * Intel oneAPI MKL "Direct Sparse Solver" (aka "DSS") solver for non-MPI builds, as well as
     * for the Intel oneAPI MKL "Parallel Direct Sparse Solver for Clusters" for MPI builds.
     * The reason why both backends are included in this wrapper class is that the latter parallel
     * solver backend cannot be used in non-MPI builds and trying to do so would result in crashes
     * due to invalid MPI calls from within the library, therefore we also had to include the "DSS"
     * backend for non-MPI builds. To avoid confusion, this wrapper always denotes the MKL backend
     * as "MKL-DSS" and does not distinguish between the two MKL backends to avoid user confusion.
     *
     * <u><b>Backend Selection Process:</b></u>\n
     * The actual backend selection process is performed during the execution of the init_symbolic()
     * method and the selection process depends on three main factors:
     * - which solver backend third-party libraries are available
     * - whether the user has specified (a deque of) allowed backends and, if so, what the user
     *   has actually specified
     * - the value of Backend::get_preferred_backend() at the time of the construction of this
     *   object, i.e. at the time of its constructor call
     *
     * Independently of all other factors, a backend can only be selected if its corresponding
     * third-party library is available for obvious reasons. Then, there are two main scenarios:
     *
     * <b>Scenario A:</b> If the user specified a list of allowed backend (combinations) by calling
     * one of the push_allowed_backend_list() or push_allowed_backends() functions, then the
     * algorithm will start with the first allowed backend (combination) from the list and it will
     * check if any of the allowed backends in the combination is available and, if so, it will
     * select that backend. If the backend combination contains more than a single allowed backend,
     * the algorithm will check them in the following order:
     * -# cuDSS
     * -# MKL-DSS
     * -# MUMPS
     * -# UMFPACK
     * -# SuperLU
     *
     * If none of the backends in a given backend combination is available, the algorithm continues
     * with the next combination in the list until it has either found a suitable backend or until
     * it reaches the end of the list, in which case the init_symbolic() function throws an
     * instance of the DSSBackendNotFoundException exception to indicate that it failed to find
     * a suitable backend that meets the users list of allowed backends.
     *
     * As a simple example, assume the user specified <c>$"mkldss|superlu cudss umfpack"</c> as a
     * list of allowed backend combinations, then the algorithm will first try to create the MKL-DSS
     * or the SuperLU backends and, if none of the two are available, it will then try to create
     * the cuDSS backend and, if that isn't available either, it will finally try to create the
     * UMFPACK backend.
     *
     * <b>Scenario B:</b> If the user did \b not specify any allowed backend combinations, then
     * the algorithm will check which PreferredBackend was selected at the time that the
     * DirectSparseSolver object was constructed and it will perform the following steps in the
     * following order
     * -# If the preferred backend was PreferredBackend::cuda, then the algorithm will try to create
     *    a cuDSS backend. If the cuDSS backend is not available -- which could happen if FEAT was
     *    compiled with CUDA enabled but without the additional cuDSS third party library -- a
     *    DSSBackendNotFoundException is thrown.
     * -# If the preferred backend was PreferredBackend::mkl, then the algorithm will try to create
     *    a MKL-DSS backend. If the MKL-DSS backend is not available -- which should never happen --
     *    a DSSBackendNotFoundException is thrown.
     * -# If the preferred backend is any other than the ones listed above, then the algorithm will
     *    try to select a suitable backend by trying out all available backends in the following
     *    order:
     *    -# cuDSS
     *    -# MKL-DSS
     *    -# MUMPS
     *    -# UMFPACK
     *    -# SuperLU
     * -# If none of the above backends are available, it throws a DSSBackendNotFoundException.
     *
     * Note: If you want to skip the check of the preferred backend and you just want the algorithm
     * to select \e any suitable backend, then you can simply push DSSBackend::all as the allowed
     * backend.
     *
     * \todo check for unsupported filters in is_available()
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class DirectSparseSolver :
      public ADPSolverBase<Matrix_, Filter_>
    {
    public:
      /// our base class
      typedef ADPSolverBase<Matrix_, Filter_> BaseClass;
      /// our (global or local) system matrix type
      typedef typename BaseClass::MatrixType MatrixType;
      /// our (global or local) system vector type
      typedef typename BaseClass::VectorType VectorType;
      /// our (global or local) system filter type
      typedef typename BaseClass::FilterType FilterType;

      /// specifies whether the cuDSS solver backend is available
      static constexpr bool have_backend_cudss = DSS::have_cudss;
      /// specifies whether the MKL-DSS solver backend is available
      static constexpr bool have_backend_mkldss = DSS::have_mkldss;
      /// specifies whether the MUMPS solver backend is available
      static constexpr bool have_backend_mumps = DSS::have_mumps;
      /// specifies whether the SuperLU solver backend is available
      static constexpr bool have_backend_superlu = DSS::have_superlu;
      /// specifies whether the UMFPACK solver backend is available
      static constexpr bool have_backend_umfpack = DSS::have_umfpack;
      /// specifies whether at least one backend for local systems is available
      static constexpr bool have_backend_local = DSS::have_backend_local;
      /// specifies whether at least one backend for global systems is available
      static constexpr bool have_backend_global = DSS::have_backend_global;

    private:
      /// backends that the user allowed us to use
      std::deque<DSSBackend> _allowed_backends;
      /// preferred backend of runtime at time of object creation
      PreferredBackend _preferred_backend;
      /// Nvidia cuDSS core
      void* _core_cudss;
      /// Intel MKLDSS core
      void* _core_mkldss;
      /// MUMPS core
      void* _core_mumps;
      /// SuperLU core
      void* _core_superlu;
      /// UMFPACK core
      void* _core_umfpack;

    public:
      /**
       * \brief Creates a new DirectSparseSolver object for the given system
       *
       * \param[in] matrix
       * A \resident reference to the system matrix.
       *
       * \param[in] filter
       * A \resident reference to the system filter.
       *
       * \note This constructor queries the result of Backend::get_preferred_backend(), which
       * is used during the symbolic initialization phase to determine which of the available
       * solver backends to choose if the user did not manually specify a set of allowed backends.
       */
      explicit DirectSparseSolver(const MatrixType& matrix, const FilterType& filter) :
        BaseClass(matrix, filter),
        _allowed_backends(),
        _preferred_backend(Backend::get_preferred_backend()),
        _core_cudss(nullptr),
        _core_mkldss(nullptr),
        _core_mumps(nullptr),
        _core_superlu(nullptr),
        _core_umfpack(nullptr)
      {
      }

      /// virtual destructor
      virtual ~DirectSparseSolver()
      {
      }

      /// Returns the name of the solver
      virtual String name() const override
      {
        String s = "DirectSparseSolver";
        if(_core_cudss)
          s += " [cuDSS]";
        if(_core_mkldss)
          s += " [MKL-DSS]";
        if(_core_mumps)
          s += " [MUMPS]";
        if(_core_superlu)
          s += " [SuperLU]";
        if(_core_umfpack)
          s += " [UMFPACK]";
        return s;
      }

      /**
       * \brief Checks whether a suitable backend for the given matrix and filter
       *
       * \param[in] matrix
       * A \transient reference to the system matrix
       *
       * \param[in] filter
       * A \transient reference to the system filter
       *
       * \returns \c true if at least one direct sparse solver backend is available for this
       * system, otherwise \c false
       */
      static bool is_available(const MatrixType& matrix, const FilterType& DOXY(filter))
      {
        // For purely local matrices, everything except for SuperLU_DIST is applicable
        if constexpr (MatrixType::is_local)
          return have_backend_local;

        if constexpr (MatrixType::is_global)
        {
          // potentially parallel case: check if we have a comm with more than 1 process
          const Dist::Comm* comm = matrix.get_comm();
          if((comm == nullptr) || (comm->size() <= 1))
            return have_backend_local;

          // MPI-parallel case with more than 1 process
          return have_backend_global;
        }

        // we should never come out here
        XABORTM("matrix type is neither local nor global!");
        return false;
      }

      /**
       * \brief Pushes a list of allowed backends into the backend queue
       *
       * \param[in] backends
       * A string containing a list of backends combinations separated by whitespaces
       */
      void push_allowed_backend_list(const String& backends)
      {
        std::deque<String> backs = backends.split_by_whitespaces();
        push_allowed_backend_list(backs);
      }

      /**
       * \brief Pushes a deque of allowed backends into the backend queue
       *
       * \param[in] backends
       * A \transient reference to a deque of backend strings to push
       */
      void push_allowed_backend_list(const std::deque<String>& backends)
      {
        for(const auto& s : backends)
          push_allowed_backend(s);
      }

      /**
       * \brief Pushes a combination of allowed backends into the backend queue
       *
       * \param[in] backend
       * A combination of backends separated by a vertical line character |
       *
       * The following strings can be used to select the corresponding backends:
       * - DSSBackend::none: "none"
       * - DSSBackend::all: "all"
       * - DSSBackend::cudss: "cudss"
       * - DSSBackend::mkldss: "mkldss" or "mkl-dss"
       * - DSSBackend::mumps: "mumps"
       * - DSSBackend::superlu: "superlu"
       * - DSSBackend::umfpack: "umfpack"
       */
      void push_allowed_backend(const String& backend)
      {
        DSSBackend b = DSSBackend::none;
        std::deque<String> sbacks = backend.split_by_charset("|");
        for(const String& s : sbacks)
        {
          if(s.compare_no_case("none") == 0)
            b = b | DSSBackend::none;
          else if(s.compare_no_case("all") == 0)
            b = b | DSSBackend::all;
          else if(s.compare_no_case("cudss") == 0)
            b = b | DSSBackend::cudss;
          else if(s.compare_no_case("mkldss") == 0)
            b = b | DSSBackend::mkldss;
          else if(s.compare_no_case("mkl-dss") == 0)
            b = b | DSSBackend::mkldss;
          else if(s.compare_no_case("mumps") == 0)
            b = b | DSSBackend::mumps;
          else if(s.compare_no_case("superlu") == 0)
            b = b | DSSBackend::superlu;
          else if(s.compare_no_case("umfpack") == 0)
            b = b | DSSBackend::umfpack;
          else
            throw DirectSparseSolverException("DirectSparseSolver: unknown DSS backend: '" + backend + "'");
        }
        push_allowed_backend(b);
      }

      /**
       * \brief Pushes a combination of allowed backends into the backend queue
       *
       * \param[in] backend
       * A combination of backends
       */
      void push_allowed_backend(DSSBackend backend)
      {
        _allowed_backends.push_back(backend);
      }

      /// Clears the deque of allowed backends
      void clear_allowed_backends()
      {
        _allowed_backends.clear();
      }

      /// Returns a reference to the deque of allowed backends
      const std::deque<DSSBackend> get_allowed_backends() const
      {
        return _allowed_backends;
      }

      /**
       * \brief Returns the selected backend
       *
       * Note that the backend is selected during the execution of init_symbolic() and therefore
       * this function only returns the actually selected backend after init_symbolic() was called.
       * This function always returns DSSBackend::none before the call to init_symbolic() and after
       * the call to done_symbolic().
       */
      DSSBackend get_selected_backend() const
      {
        if(this->_core_umfpack)
          return DSSBackend::umfpack;
        if(this->_core_superlu)
          return DSSBackend::superlu;
        if(this->_core_mkldss)
          return DSSBackend::mkldss;
        if(this->_core_mumps)
          return DSSBackend::mumps;
        if(this->_core_cudss)
          return DSSBackend::cudss;
        return DSSBackend::none;
      }

      /**
       * \brief Performs the symbolic initialization of the solver
       *
       * This function selects the actual backend before starting the symbolic initialization.
       * Please check the documentation of this class for details about how the backend selection
       * process is performed.
       *
       * This function may throw instances of DirectSparseSolver exceptions if anything goes wrong.
       */
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // First of all, let's see if we're in a single-process scenario here
        const Dist::Comm* comm = this->_get_comm();
        const bool single_process = ((comm == nullptr) || (comm->size() <= 1));

        // if the user didn't specify any allowed backends, then we'll try a default setting
        if(_allowed_backends.empty())
        {
          // First, let's see if the user wanted a cuda-based solver
          if(_preferred_backend == PreferredBackend::cuda)
          {
            // Ok, try to create a cuDSS core then
            if(!_create_core_cudss())
            {
              // Nope, apparently, cuDSS is not available
              throw DSSBackendNotFoundException(DSSBackend::cudss);
            }
          }

          // Next, let's see if the user wanted a MKL-based solver
          if(_preferred_backend == PreferredBackend::mkl)
          {
            // Ok, try to create a MKL-DSS core then
            if(!_create_core_mkldss())
            {
              // Nope, apparently, MKL-DSS is not available
              throw DSSBackendNotFoundException(DSSBackend::mkldss);
            }
          }

          // Now we're left with the generic case, in which we may choose any backend we see fit

          // Let's try cuDSS first
          if(_create_core_cudss())
            return;

          // Let's try MKL-DSS next
          if(_create_core_mkldss())
            return;

          // Let's try UMFPACK if we're in a single-process case next
          if(single_process && _create_core_umfpack())
            return;

          // Let's try MUMPS next
          if(_create_core_mumps())
            return;

          // Let's try SuperLU last
          if(_create_core_superlu())
            return;

          // If we come out here, then we have failed to find a suitable backend...
          throw DSSBackendNotFoundException();
        }

        // Let's loop over all allowed backends
        for(auto ab : _allowed_backends)
        {
          // Let's try cuDSS first (if allowed)
          if(*(ab & DSSBackend::cudss) && _create_core_cudss())
            return;

          // Let's try MKL-DSS next (if allowed)
          if(*(ab & DSSBackend::mkldss) && _create_core_mkldss())
            return;

          // Let's try UMFPACK if we're in a single-process case next (if allowed)
          if(*(ab & DSSBackend::umfpack) && single_process && _create_core_umfpack())
            return;

          // Let's try MUMPS next (if allowed)
          if(*(ab & DSSBackend::mkldss) && _create_core_mumps())
            return;

          // Let's try SuperLU last (if allowed)
          if(*(ab & DSSBackend::superlu) && _create_core_superlu())
            return;
        }

        // If we come out here, then we have failed to find a suitable backend...
        throw DSSBackendNotFoundException();
      }

      /**
       * \brief Releases the symbolic factorization data and frees the backend
       */
      virtual void done_symbolic() override
      {
#ifdef FEAT_HAVE_CUDSS
        if(this->_core_cudss)
        {
          DSS::destroy_cudss_core(this->_core_cudss);
          this->_core_cudss = nullptr;
        }
#endif // FEAT_HAVE_CUDSS
#ifdef FEAT_HAVE_MKL
        if(this->_core_mkldss)
        {
          DSS::destroy_mkldss_core(this->_core_mkldss);
          this->_core_mkldss = nullptr;
        }
#endif // FEAT_HAVE_MKL
#ifdef FEAT_HAVE_MUMPS
        if(this->_core_mumps)
        {
          DSS::destroy_mumps_core(this->_core_mumps);
          this->_core_mumps = nullptr;
        }
#endif // FEAT_HAVE_MUMPS
#ifdef FEAT_HAVE_SUPERLU_DIST
        if(this->_core_superlu)
        {
          DSS::destroy_superlu_core(this->_core_superlu);
          this->_core_superlu = nullptr;
        }
#endif // FEAT_HAVE_SUPERLU_DIST
#ifdef FEAT_HAVE_UMFPACK
        if(this->_core_umfpack)
        {
          DSS::destroy_umfpack_core(this->_core_umfpack);
          this->_core_umfpack = nullptr;
        }
#endif // FEAT_HAVE_UMFPACK
        BaseClass::done_symbolic();
      }

      /**
       * \brief Performs the numeric initialization of the solver
       *
       * This function may throw instances of DirectSparseSolver exceptions if anything goes wrong.
       */
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
#ifdef FEAT_HAVE_CUDSS
        if(this->_core_cudss)
        {
          this->_upload_numeric(
            DSS::get_cudss_mat_val(this->_core_cudss),
            DSS::get_cudss_row_ptr(this->_core_cudss),
            DSS::get_cudss_col_idx(this->_core_cudss));

          DSS::init_cudss_numeric(this->_core_cudss);
        }
#endif // FEAT_HAVE_CUDSS
#ifdef FEAT_HAVE_MKL
        if(this->_core_mkldss)
        {
          this->_upload_numeric(
            DSS::get_mkldss_mat_val(this->_core_mkldss),
            DSS::get_mkldss_row_ptr(this->_core_mkldss),
            DSS::get_mkldss_col_idx(this->_core_mkldss));

          DSS::init_mkldss_numeric(this->_core_mkldss);
        }
#endif // FEAT_HAVE_MKL
#ifdef FEAT_HAVE_MUMPS
        if(this->_core_mumps)
        {
          this->_upload_numeric(
            DSS::get_mumps_mat_val(this->_core_mumps),
            DSS::get_mumps_row_ptr(this->_core_mumps),
            DSS::get_mumps_col_idx(this->_core_mumps));

          DSS::init_mumps_numeric(this->_core_mumps);
        }
#endif // FEAT_HAVE_MUMPS
#ifdef FEAT_HAVE_SUPERLU_DIST
        if(this->_core_superlu)
        {
          this->_upload_numeric(
            DSS::get_superlu_mat_val(this->_core_superlu),
            DSS::get_superlu_row_ptr(this->_core_superlu),
            DSS::get_superlu_col_idx(this->_core_superlu));

          DSS::init_superlu_numeric(this->_core_superlu);
        }
#endif // FEAT_HAVE_SUPERLU_DIST
#ifdef FEAT_HAVE_UMFPACK
        if(this->_core_umfpack)
        {
          this->_upload_numeric(
            DSS::get_umfpack_mat_val(this->_core_umfpack),
            DSS::get_umfpack_row_ptr(this->_core_umfpack),
            DSS::get_umfpack_col_idx(this->_core_umfpack));

          DSS::init_umfpack_numeric(this->_core_umfpack);
        }
#endif // FEAT_HAVE_UMFPACK
      }

      /**
       * \brief Releases the numeric initialization data of the solver
       */
      virtual void done_numeric() override
      {
#ifdef FEAT_HAVE_SUPERLU_DIST
        if(this->_core_superlu)
        {
          DSS::done_superlu_numeric(this->_core_superlu);
        }
#endif // FEAT_HAVE_SUPERLU_DIST
#ifdef FEAT_HAVE_UMFPACK
        if(this->_core_umfpack)
        {
          DSS::done_umfpack_numeric(this->_core_umfpack);
        }
#endif // FEAT_HAVE_UMFPACK
        BaseClass::done_numeric();
      }

      /**
       * \brief Solves a linear system with the factorized system matrix.
       *
       * \param[in,out] vec_cor
       * A reference to the solution vector. The vector must be allocated to the correct length, but its
       * initial contents are ignored.
       *
       * \param[in] vec_def
       * A reference to the right-hand-side of the linear system.
       */
      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // silence compiler warnings about unused variables if no backend is available
        (void)vec_cor;
        (void)vec_def;

#ifdef FEAT_HAVE_CUDSS
        if(this->_core_cudss)
        {
          // upload defect vector
          this->_upload_vector(DSS::get_cudss_rhs_val(this->_core_cudss), vec_def);

          // solve system via MKL DSS
          try
          {
            DSS::solve_cudss(this->_core_cudss);
          }
          catch(const DirectSparseSolverException&)
          {
            return Status::aborted;
          }
          catch(...)
          {
            throw;
          }

          // download correction vector
          this->_download_vector(vec_cor, DSS::get_cudss_sol_val(this->_core_cudss));

          return Status::success;
        }
#endif // FEAT_HAVE_CUDSS
#ifdef FEAT_HAVE_MKL
        if(this->_core_mkldss)
        {
          // upload defect vector
          this->_upload_vector(DSS::get_mkldss_rhs_val(this->_core_mkldss), vec_def);

          // solve system via MKL DSS
          try
          {
            DSS::solve_mkldss(this->_core_mkldss);
          }
          catch(const DirectSparseSolverException&)
          {
            return Status::aborted;
          }
          catch(...)
          {
            throw;
          }

          // download correction vector
          this->_download_vector(vec_cor, DSS::get_mkldss_sol_val(this->_core_mkldss));

          return Status::success;
        }
#endif // FEAT_HAVE_MKL
#ifdef FEAT_HAVE_MUMPS
        if(this->_core_mumps)
        {
          // upload defect vector
          this->_upload_vector(DSS::get_mumps_vector(this->_core_mumps), vec_def);

          // solve system via MUMPS
          try
          {
            DSS::solve_mumps(this->_core_mumps);
          }
          catch(const DirectSparseSolverException&)
          {
            return Status::aborted;
          }
          catch(...)
          {
            throw;
          }

          // download correction vector
          this->_download_vector(vec_cor, DSS::get_mumps_vector(this->_core_mumps));

          return Status::success;
        }
#endif // FEAT_HAVE_MUMPS
#ifdef FEAT_HAVE_SUPERLU_DIST
        if(this->_core_superlu)
        {
          // upload defect vector
          this->_upload_vector(DSS::get_superlu_vector(this->_core_superlu), vec_def);

          // solve system via SuperLU
          try
          {
            DSS::solve_superlu(this->_core_superlu);
          }
          catch(const DirectSparseSolverException&)
          {
            return Status::aborted;
          }
          catch(...)
          {
            throw;
          }

          // download correction vector
          this->_download_vector(vec_cor, DSS::get_superlu_vector(this->_core_superlu));

          return Status::success;
        }
#endif // FEAT_HAVE_SUPERLU_DIST
#ifdef FEAT_HAVE_UMFPACK
        if(this->_core_umfpack)
        {
          // upload defect vector
          this->_upload_vector(DSS::get_umfpack_rhs_val(this->_core_umfpack), vec_def);

          // solve system via UMFPACK
          try
          {
            DSS::solve_umfpack(this->_core_umfpack);
          }
          catch(const DirectSparseSolverException&)
          {
            return Status::aborted;
          }
          catch(...)
          {
            throw;
          }

          // download correction vector
          this->_download_vector(vec_cor, DSS::get_umfpack_sol_val(this->_core_umfpack));

          return Status::success;
        }
#endif // FEAT_HAVE_UMFPACK

        return Status::aborted;
      }

    private:
      /// Tries to create a cuDSS backend core and returns true, if it succeeded
      bool _create_core_cudss()
      {
#ifdef FEAT_HAVE_CUDSS
        // create cuDSS core
        this->_core_cudss = DSS::create_cudss_core(
          this->_get_comm(),
          this->_get_num_global_dofs(),
          this->_get_global_dof_offset(),
          this->_get_num_owned_dofs(),
          this->_get_adp_matrix_num_nzes(),
          this->_get_num_global_nonzeros());

        // upload matrix structure to cuDSS core
        BaseClass::_upload_symbolic(
          DSS::get_cudss_row_ptr(this->_core_cudss),
          DSS::get_cudss_col_idx(this->_core_cudss));

        // perform symbolic factorization of cuDSS core
        DSS::init_cudss_symbolic(this->_core_cudss);

        return true;
#else
        return false;
#endif // FEAT_HAVE_CUDSS
      }

      /// Tries to create a MKL-DSS backend core and returns true, if it succeeded
      bool _create_core_mkldss()
      {
#ifdef FEAT_HAVE_MKL
        // create MKL DSS core
        this->_core_mkldss = DSS::create_mkldss_core(
          this->_get_comm(),
          this->_get_num_global_dofs(),
          this->_get_global_dof_offset(),
          this->_get_num_owned_dofs(),
          this->_get_adp_matrix_num_nzes(),
          this->_get_num_global_nonzeros());

        // upload matrix structure to MKL DSS core
        BaseClass::_upload_symbolic(
          DSS::get_mkldss_row_ptr(this->_core_mkldss),
          DSS::get_mkldss_col_idx(this->_core_mkldss));

        // perform symbolic factorization of MKL DSS core
        DSS::init_mkldss_symbolic(this->_core_mkldss);

        return true;
#else
        return false;
#endif // FEAT_HAVE_MKL
      }

      /// Tries to create a MUMPS backend core and returns true, if it succeeded
      bool _create_core_mumps()
      {
#ifdef FEAT_HAVE_MUMPS
        // create MUMPS core
        this->_core_mumps = DSS::create_mumps_core(
          this->_get_comm(),
          this->_get_num_global_dofs(),
          this->_get_global_dof_offset(),
          this->_get_num_owned_dofs(),
          this->_get_adp_matrix_num_nzes(),
          this->_get_num_global_nonzeros());

        // upload matrix structure to MUMPS core
        BaseClass::_upload_symbolic(
          DSS::get_mumps_row_ptr(this->_core_mumps),
          DSS::get_mumps_col_idx(this->_core_mumps));

        // perform symbolic factorization of MUMPS core
        DSS::init_mumps_symbolic(this->_core_mumps);

        return true;
#else
        return false;
#endif // FEAT_HAVE_MUMPS
      }

      /// Tries to create a SuperLU backend core and returns true, if it succeeded
      bool _create_core_superlu()
      {
#ifdef FEAT_HAVE_SUPERLU_DIST
        // create SuperLU core
        this->_core_superlu = DSS::create_superlu_core(
          this->_get_comm(),
          this->_get_num_global_dofs(),
          this->_get_global_dof_offset(),
          this->_get_num_owned_dofs(),
          this->_get_adp_matrix_num_nzes(),
          this->_get_num_global_nonzeros());

        // upload matrix structure to SuperLU core
        BaseClass::_upload_symbolic(
          DSS::get_superlu_row_ptr(this->_core_superlu),
          DSS::get_superlu_col_idx(this->_core_superlu));

        // perform symbolic factorization of SuperLU core
        DSS::init_superlu_symbolic(this->_core_superlu);

        return true;
#else
        return false;
#endif // FEAT_HAVE_SUPERLU_DIST
      }

      /// Tries to create an UMFPACK backend core and returns true, if it succeeded
      bool _create_core_umfpack()
      {
#ifdef FEAT_HAVE_UMFPACK
        // create UMFPACK core
        this->_core_umfpack = DSS::create_umfpack_core(
          this->_get_comm(),
          this->_get_num_global_dofs(),
          this->_get_global_dof_offset(),
          this->_get_num_owned_dofs(),
          this->_get_adp_matrix_num_nzes(),
          this->_get_num_global_nonzeros());

        // upload matrix structure to UMFPACK core
        BaseClass::_upload_symbolic(
          DSS::get_umfpack_row_ptr(this->_core_umfpack),
          DSS::get_umfpack_col_idx(this->_core_umfpack));

        // perform symbolic factorization of UMFPACK core
        DSS::init_umfpack_symbolic(this->_core_umfpack);

        return true;
#else
        return false;
#endif // FEAT_HAVE_UMFPACK
      }
    }; // class DirectSparseSolver<...>

    /**
     * \brief Creates a new DirectSparseSolver object
     *
     * \param[in] matrix
     * A \resident reference to the system matrix.
     *
     * \param[in] filter
     * A \resident reference to the system filter.
     *
     * \returns
     * A shared pointer to a new DirectSparseSolver object.
     */
    template<typename Matrix_, typename Filter_>
    std::shared_ptr<DirectSparseSolver<Matrix_, Filter_>> new_direct_sparse_solver(const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<DirectSparseSolver<Matrix_, Filter_>>(matrix, filter);
    }

    /**
     * \brief Checks whether at least one DirectSparseSolver backend is available
     *
     * \param[in] matrix
     * A \transient reference to the system matrix.
     *
     * \param[in] filter
     * A \transient reference to the system filter.
     *
     * \returns \c true, if at least one backend is available for the given system,
     * otherwise \c false.
     */
    template<typename Matrix_, typename Filter_>
    bool have_direct_sparse_solver(const Matrix_& matrix, const Filter_& filter)
    {
      return DirectSparseSolver<Matrix_, Filter_>::is_available(matrix, filter);
    }
  } // namespace Solver
} // namespace FEAT
