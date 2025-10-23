// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/solver/direct_sparse_solver.hpp>

#ifdef FEAT_HAVE_UMFPACK

FEAT_DISABLE_WARNINGS
#include <umfpack.h>
FEAT_RESTORE_WARNINGS

// use 32-bit indices for UMFPACK
#define FEAT_DSS_UMFPACK_I32

namespace FEAT
{
  namespace Solver
  {
    namespace DSS
    {
#ifdef FEAT_DSS_UMFPACK_I64
      static_assert(sizeof(UMFPACK_IT) == 8, "DirectSparseSolver: UMFPACK: index type size mismatch!");
#else
      static_assert(sizeof(UMFPACK_IT) == 4, "DirectSparseSolver: UMFPACK: index type size mismatch!");
#endif

      /**
       * \brief UMFPACK core wrapper class
       *
       * \author Peter Zajac
       */
      class UMFPACK_Core
      {
      public:
        // UMFPACK control array
        double control[UMFPACK_CONTROL];
        // UMFPACK info array
        double info[UMFPACK_INFO];
        // UMFPACK opaque symbolic object
        void* symbolic;
        // UMFPACK opaque numeric object
        void* numeric;

        UMFPACK_IT num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes;
        std::vector<UMFPACK_IT> row_ptr, col_idx;
        std::vector<UMFPACK_DT> mat_val, rhs_val, sol_val;

        explicit UMFPACK_Core(const Dist::Comm& comm, Index num_global_dofs_, Index dof_offset_,
          Index num_owned_dofs_, Index num_owned_nzes_, Index DOXY(num_global_nzes_)) :
          symbolic(nullptr),
          numeric(nullptr),
          num_global_dofs(UMFPACK_IT(num_global_dofs_)),
          dof_offset(UMFPACK_IT(dof_offset_)),
          num_owned_dofs(UMFPACK_IT(num_owned_dofs_)),
          num_owned_nzes(UMFPACK_IT(num_owned_nzes_))
        {
#ifdef FEAT_HAVE_MPI
          XASSERTM(comm.size() == 1, "DirectSparseSolver UMFPACK core can only be used on 1 MPI process!");
#else
          // prevent unused parameter warnings
          (void)comm;
#endif

          // allocate matrix and vector arrays
          row_ptr.resize(num_owned_dofs+1u, 0u);
          col_idx.resize(num_owned_nzes, 0u);
          mat_val.resize(num_owned_nzes, 0.0);
          rhs_val.resize(num_owned_dofs, 0.0);
          sol_val.resize(num_owned_dofs, 0.0);

#ifdef FEAT_DSS_UMFPACK_I64
          ::umfpack_dl_defaults(control);
#else
          ::umfpack_di_defaults(control);
#endif
          control[UMFPACK_PRL] = 1;
        }

        ~UMFPACK_Core()
        {
          if(numeric != nullptr)
          {
#ifdef FEAT_DSS_UMFPACK_I64
            ::umfpack_dl_free_numeric(&numeric);
#else
            ::umfpack_di_free_numeric(&numeric);
#endif
          }
          if(symbolic != nullptr)
          {
#ifdef FEAT_DSS_UMFPACK_I64
            ::umfpack_dl_free_symbolic(&symbolic);
#else
            ::umfpack_di_free_symbolic(&symbolic);
#endif
          }
        }

        void init_symbolic()
        {
          XASSERT(symbolic == nullptr);

          int status =
#ifdef FEAT_DSS_UMFPACK_I64
          ::umfpack_dl_symbolic
#else
          ::umfpack_di_symbolic
#endif
            (num_owned_dofs, num_owned_dofs, row_ptr.data(), col_idx.data(), mat_val.data(), &symbolic, control, info);

          // check status code
          switch(status)
          {
          case UMFPACK_OK:
            break;
          case UMFPACK_ERROR_out_of_memory:
            throw DirectSparseSolverException("UMFPACK", "out of memory");
          case UMFPACK_ERROR_invalid_matrix:
            throw DirectSparseSolverException("UMFPACK", "invalid matrix structure");
          case UMFPACK_ERROR_internal_error:
            throw DirectSparseSolverException("UMFPACK", "internal error");
          default:
            throw DirectSparseSolverException("UMFPACK", "unknown error");
          }
        }

        void init_numeric()
        {
          XASSERT(symbolic != nullptr);
          XASSERT(numeric == nullptr);

          int status =
#ifdef FEAT_DSS_UMFPACK_I64
            ::umfpack_dl_numeric
#else
            ::umfpack_di_numeric
#endif
            (row_ptr.data(), col_idx.data(), mat_val.data(), symbolic, &numeric, control, info);

          // check status code
          switch(status)
          {
          case UMFPACK_OK:
            break;
          case UMFPACK_ERROR_out_of_memory:
            throw DirectSparseSolverException("UMFPACK", "out of memory");
          case UMFPACK_ERROR_invalid_matrix:
            throw DirectSparseSolverException("UMFPACK", "invalid matrix structure");
          case UMFPACK_WARNING_singular_matrix:
            throw DirectSparseSolverException("UMFPACK", "singular matrix");
          case UMFPACK_ERROR_different_pattern:
            throw DirectSparseSolverException("UMFPACK", "different pattern");
          case UMFPACK_ERROR_internal_error:
            throw DirectSparseSolverException("UMFPACK", "internal error");
          default:
            throw DirectSparseSolverException("UMFPACK", "unknown error");
          }
        }

        void done_numeric()
        {
          XASSERT(numeric != nullptr);
#ifdef FEAT_DSS_UMFPACK_I64
          ::umfpack_dl_free_numeric(&numeric);
#else
          ::umfpack_di_free_numeric(&numeric);
#endif
          numeric = nullptr;
        }

        void solve()
        {
          int status =
#ifdef FEAT_DSS_UMFPACK_I64
          ::umfpack_dl_solve
#else
          ::umfpack_di_solve
#endif
            (UMFPACK_At, row_ptr.data(), col_idx.data(), mat_val.data(), sol_val.data(), rhs_val.data(), numeric, control, info);

          // check status code
          switch(status)
          {
          case UMFPACK_OK:
            break;
          case UMFPACK_ERROR_out_of_memory:
            throw DirectSparseSolverException("UMFPACK", "out of memory");
          case UMFPACK_WARNING_singular_matrix:
            throw DirectSparseSolverException("UMFPACK", "singular matrix");
          case UMFPACK_ERROR_internal_error:
            throw DirectSparseSolverException("UMFPACK", "internal error");
          default:
            throw DirectSparseSolverException("UMFPACK", "unknown error");
          }
        }
      };

      void* create_umfpack_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes)
      {
        return new UMFPACK_Core(*comm, num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes, num_global_nzes);
      }

      void destroy_umfpack_core(void* core)
      {
        XASSERT(core != nullptr);
        delete reinterpret_cast<UMFPACK_Core*>(core);
      }

      UMFPACK_IT* get_umfpack_row_ptr(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<UMFPACK_Core*>(core)->row_ptr.data();
      }

      UMFPACK_IT* get_umfpack_col_idx(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<UMFPACK_Core*>(core)->col_idx.data();
      }

      UMFPACK_DT* get_umfpack_mat_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<UMFPACK_Core*>(core)->mat_val.data();
      }

      UMFPACK_DT* get_umfpack_rhs_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<UMFPACK_Core*>(core)->rhs_val.data();
      }

      UMFPACK_DT* get_umfpack_sol_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<UMFPACK_Core*>(core)->sol_val.data();
      }

      void init_umfpack_symbolic(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<UMFPACK_Core*>(core)->init_symbolic();
      }

      void init_umfpack_numeric(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<UMFPACK_Core*>(core)->init_numeric();
      }

      void done_umfpack_numeric(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<UMFPACK_Core*>(core)->done_numeric();
      }

      void solve_umfpack(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<UMFPACK_Core*>(core)->solve();
      }
    } // namespace DSS
  } // namespace Solver
} // namespace FEAT

#else // no FEAT_HAVE_UMFPACK

void feat_direct_sparse_solver_umfpack_dummy()
{
}

#endif // FEAT_HAVE_UMFPACK
