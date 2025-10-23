// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>

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
    namespace SuperLU_Aux
    {
      /// This is the actual wrapper class around the SuperLU structures.
      class Core
      {
      public:
        typedef int_t IT;
        typedef double DT;

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
        const IT slu_dof_offset, slu_num_owned_dofs, slu_num_global_dofs, slu_num_nonzeros;

        /// row-pointer and column-index arrays; we need a backup of the column-index array
        /// because SuperLU permutes the first one without asking for our permission first :-/
        std::vector<IT> slu_row_ptr, slu_col_idx, slu_col_idx2;
        /// matrix value array
        std::vector<double> slu_mat_val;
        /// temporary vector
        std::vector<double> slu_vector;
        /// temporary vector
        std::vector<double> slu_v_berr;

        /// remember whether we already performed a numeric factorization; if so, then we can
        /// reuse the symbolic factorization of the pattern in the next numeric factorization
        bool was_factorized;

      public:
        explicit Core(const void* comm, Index num_global_dofs, Index my_dof_offset,
          Index num_owned_dofs, Index num_nonzeros) :
          _comm(*reinterpret_cast<const MPI_Comm*>(comm)),
          slu_info(0),
          slu_dof_offset(IT(my_dof_offset)),
          slu_num_owned_dofs(IT(num_owned_dofs)),
          slu_num_global_dofs(IT(num_global_dofs)),
          slu_num_nonzeros(IT(num_nonzeros)),
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
          slu_col_idx.resize(num_nonzeros);
          slu_col_idx2.resize(num_nonzeros);
          slu_mat_val.resize(num_nonzeros);
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

        ~Core()
        {
          PStatFree(&slu_stats);
          dDestroy_LU(slu_matrix.nrow, &slu_grid, &slu_lu_struct);
          dLUstructFree(&slu_lu_struct);
          dScalePermstructFree(&slu_scale_perm);
          dSolveFinalize(&slu_opts, &slu_solve_struct);
          superlu_gridexit(&slu_grid);
        }

        void init_symbolic()
        {
          // copy column indices arrays, because SuperLU permutes them
          memcpy(slu_col_idx2.data(), slu_col_idx.data(), sizeof(IT)*slu_col_idx.size());

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

        int init_numeric()
        {
          // reset column-index array because it might have been
          // overwritten by the previous factorization call
          memcpy(slu_col_idx2.data(), slu_col_idx.data(), sizeof(IT)*slu_col_idx.size());

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

          // okay
          return slu_info;
        }

        int solve()
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

          // return the result
          return slu_info;
        }
      }; // class Core

      void* create_core(const void* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_nonzeros)
      {
        return new Core(comm, num_global_dofs, dof_offset, num_owned_dofs, num_nonzeros);
      }

      void destroy_core(void* core)
      {
        delete reinterpret_cast<Core*>(core);
      }

      int_t* get_row_ptr(void* core)
      {
        return reinterpret_cast<Core*>(core)->slu_row_ptr.data();
      }

      int_t* get_col_idx(void* core)
      {
        return reinterpret_cast<Core*>(core)->slu_col_idx.data();
      }

      double* get_mat_val(void* core)
      {
        return reinterpret_cast<Core*>(core)->slu_mat_val.data();
      }

      double* get_vector(void* core)
      {
        return reinterpret_cast<Core*>(core)->slu_vector.data();
      }

      void init_symbolic(void* core)
      {
        reinterpret_cast<Core*>(core)->init_symbolic();
      }

      int init_numeric(void* core)
      {
        return reinterpret_cast<Core*>(core)->init_numeric();
      }

      int solve(void* core)
      {
        return reinterpret_cast<Core*>(core)->solve();
      }
    } // namespace SuperLU_Aux
  } // namespace Solver
} // namespace FEAT

#else
// insert dummy function to suppress linker warnings
void dummy_superlu_function() {}
#endif // defined(FEAT_HAVE_SUPERLU_DIST)  && defined(FEAT_HAVE_MPI)
