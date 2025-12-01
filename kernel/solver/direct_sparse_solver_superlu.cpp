// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/solver/direct_sparse_solver.hpp>

#if defined(FEAT_HAVE_SUPERLU_DIST) && defined(FEAT_HAVE_MPI)

FEAT_DISABLE_WARNINGS
#define VTUNE 0
#include <superlu_ddefs.h>
FEAT_RESTORE_WARNINGS

#include <vector>

namespace FEAT
{
  namespace Solver
  {
    namespace DSS
    {
      static_assert(sizeof(SUPERLU_IT) == sizeof(int_t), "DirectSparseSolver: SuperLU: index type size mismatch!");

      /**
       * \brief SuperLU_DIST core wrapper class
       *
       * \author Peter Zajac
       */
      class SuperLU_Core
      {
      public:
        /// MPI communicator
        MPI_Comm _comm;
        /// SuperLU grid info structure
        gridinfo_t slu_grid;
        /// options for SuperLU
        superlu_dist_options_t slu_opts;
        /// LU statistics
        SuperLUStat_t slu_stats;
        /// SuperLU matrix structure
        SuperMatrix slu_matrix;
        /// actual matrix storage
        NRformat_loc slu_matrix_store;
        /// SuperLU scale permutation structure
        dScalePermstruct_t slu_scale_perm;
        /// SuperLU LU structure
        dLUstruct_t slu_lu_struct;
        /// SuperLU solve structure
        dSOLVEstruct_t slu_solve_struct;
        /// returned info of last solve call
        int slu_info;
        const int_t slu_dof_offset, slu_num_owned_dofs, slu_num_global_dofs, slu_num_nonzeros;

        /// row-pointer and column-index arrays; we need a backup of the column-index array
        /// because SuperLU permutes the first one without asking for our permission first :-/
        std::vector<int_t> slu_row_ptr, slu_col_idx, slu_col_idx2;
        /// matrix value array
        std::vector<double> slu_mat_val;
        /// temporary vector
        std::vector<double> slu_vector;
        /// temporary vector
        std::vector<double> slu_v_berr;

        /// remember whether we already performed a numeric factorization; if so, then we can
        /// reuse the symbolic factorization of the pattern in the next numeric factorization
        bool was_factorized;

        explicit SuperLU_Core(const Dist::Comm& comm, Index num_global_dofs, Index my_dof_offset,
          Index num_owned_dofs, Index num_owned_nzes, Index DOXY(num_global_nzes)) :
          _comm(comm.mpi_comm()),
          slu_info(0),
          slu_dof_offset(static_cast<int_t>(my_dof_offset)),
          slu_num_owned_dofs(static_cast<int_t>(num_owned_dofs)),
          slu_num_global_dofs(static_cast<int_t>(num_global_dofs)),
          slu_num_nonzeros(static_cast<int_t>(num_owned_nzes)),
          was_factorized(false)
        {
          // first of all, memclear all SuperLU structures just to be safe
          memset(&slu_grid, 0, sizeof(gridinfo_t));
          memset(&slu_opts, 0, sizeof(superlu_dist_options_t));
          memset(&slu_stats, 0, sizeof(SuperLUStat_t));
          memset(&slu_matrix, 0, sizeof(SuperMatrix));
          memset(&slu_matrix_store, 0, sizeof(NRformat_loc));
          memset(&slu_scale_perm, 0, sizeof(dScalePermstruct_t));
          memset(&slu_lu_struct, 0, sizeof(dLUstruct_t));
          memset(&slu_solve_struct, 0, sizeof(dSOLVEstruct_t));

          // set up grid
          int ranks(0);
          MPI_Comm_size(_comm, &ranks);
          superlu_gridinit(_comm, ranks, 1, &slu_grid);

          // allocate matrix and vector arrays
          slu_row_ptr.resize(num_owned_dofs+1u);
          slu_col_idx.resize(num_owned_nzes);
          slu_col_idx2.resize(num_owned_nzes);
          slu_mat_val.resize(num_owned_nzes);
          slu_vector.resize(num_owned_dofs);
          slu_v_berr.resize(num_owned_dofs);

          // setup SuperLU matrix
          slu_matrix.Stype = SLU_NR_loc; // CSR, local
          slu_matrix.Dtype = SLU_D;      // double
          slu_matrix.Mtype = SLU_GE;     // generic matrix
          slu_matrix.nrow = slu_num_global_dofs;
          slu_matrix.ncol = slu_num_global_dofs; // always square
          slu_matrix.Store = &slu_matrix_store;

          // setup SuperLU matrix storage
          slu_matrix_store.nnz_loc = slu_num_nonzeros;
          slu_matrix_store.m_loc   = slu_num_owned_dofs;
          slu_matrix_store.fst_row = slu_dof_offset;
          slu_matrix_store.nzval   = slu_mat_val.data();
          slu_matrix_store.rowptr  = slu_row_ptr.data();
          slu_matrix_store.colind  = slu_col_idx2.data(); // use copy here
        }

        ~SuperLU_Core()
        {
          PStatFree(&slu_stats);
          dLUstructFree(&slu_lu_struct);
          dScalePermstructFree(&slu_scale_perm);
          dSolveFinalize(&slu_opts, &slu_solve_struct);
          superlu_gridexit(&slu_grid);
        }

        void init_symbolic()
        {
          // copy column indices arrays, because SuperLU permutes them
          memcpy(slu_col_idx2.data(), slu_col_idx.data(), sizeof(int_t)*slu_col_idx.size());

          // set up default options
          set_default_options_dist(&slu_opts);

          // don't print statistics to cout
          slu_opts.PrintStat = NO;

          // replace tiny pivots
          slu_opts.ReplaceTinyPivot = YES;

          // initialize scale perm structure
          dScalePermstructInit(slu_num_global_dofs, slu_num_global_dofs, &slu_scale_perm);

          // initialize LU structure
          dLUstructInit(slu_num_global_dofs, &slu_lu_struct);

          // initialize the statistics
          PStatInit(&slu_stats);
        }

        void init_numeric()
        {
          // reset column-index array because it might have been
          // overwritten by the previous factorization call
          memcpy(slu_col_idx2.data(), slu_col_idx.data(), sizeof(int_t)*slu_col_idx.size());

          // have we already factorized before?
          slu_opts.Fact = (was_factorized ? SamePattern : DOFACT);

          // call the solver routine to factorize the matrix
          pdgssvx(
            &slu_opts,
            &slu_matrix,
            &slu_scale_perm,
            nullptr,
            int(slu_num_owned_dofs),
            0,
            &slu_grid,
            &slu_lu_struct,
            &slu_solve_struct,
            slu_v_berr.data(),
            &slu_stats,
            &slu_info);

          // factorization successful?
          if(slu_info == 0)
          {
            // okay, remember that we have factorized
            was_factorized = true;
          }
          else if(slu_info < 0)
          {
            throw InternalError(__func__, __FILE__, __LINE__, "invalid argument for SuperLU pdgssvx");
          }
          else if(slu_info < int(slu_num_global_dofs))
          {
            throw DirectSparseSolverException("SuperLU", "zero pivot");
          }
          else
          {
            throw DirectSparseSolverException("SuperLU", "out of memory");
          }
        }

        void done_numeric()
        {
          dDestroy_LU(slu_num_global_dofs, &slu_grid, &slu_lu_struct);
        }

        void solve()
        {
          // matrix was already factorized in set_matrix_values() call
          slu_opts.Fact = FACTORED;

          // call the solver routine
          pdgssvx(
            &slu_opts,
            &slu_matrix,
            &slu_scale_perm,
            slu_vector.data(),
            int(slu_num_owned_dofs),
            1,
            &slu_grid,
            &slu_lu_struct,
            &slu_solve_struct,
            slu_v_berr.data(),
            &slu_stats,
            &slu_info);

          // check return value
          if(slu_info == 0)
          {
            // ok, everthing's fine
          }
          else if(slu_info < 0)
          {
            throw InternalError(__func__, __FILE__, __LINE__, "invalid argument for SuperLU pdgssvx");
          }
          else if(slu_info < int(slu_num_global_dofs))
          {
            throw DirectSparseSolverException("SuperLU", "zero pivot");
          }
          else
          {
            throw DirectSparseSolverException("SuperLU", "out of memory");
          }
        }
      }; // class SuperLU_Core

      void* create_superlu_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes)
      {
        return new SuperLU_Core(*comm, num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes, num_global_nzes);
      }

      void destroy_superlu_core(void* core)
      {
        XASSERT(core != nullptr);
        delete reinterpret_cast<SuperLU_Core*>(core);
      }

      SUPERLU_IT* get_superlu_row_ptr(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<SuperLU_Core*>(core)->slu_row_ptr.data();
      }

      SUPERLU_IT* get_superlu_col_idx(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<SuperLU_Core*>(core)->slu_col_idx.data();
      }

      SUPERLU_DT* get_superlu_mat_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<SuperLU_Core*>(core)->slu_mat_val.data();
      }

      SUPERLU_DT* get_superlu_vector(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<SuperLU_Core*>(core)->slu_vector.data();
      }

      void init_superlu_symbolic(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<SuperLU_Core*>(core)->init_symbolic();
      }

      void init_superlu_numeric(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<SuperLU_Core*>(core)->init_numeric();
      }

      void done_superlu_numeric(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<SuperLU_Core*>(core)->done_numeric();
      }

      void solve_superlu(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<SuperLU_Core*>(core)->solve();
      }
    } // namespace DSS
  } // namespace Solver
} // namespace FEAT

#else // no FEAT_HAVE_SUPERLU_DIST

void feat_direct_sparse_solver_superlu_dummy()
{
}

#endif // FEAT_HAVE_SUPERLU_DIST
