// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
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

        /// row-pointer and column-index arrays; we need a backup of the column-index array
        /// because SuperLU permutes the first one without asking for our permission first :-/
        std::vector<int_t> slu_row_ptr, slu_col_idx, slu_col_idx2;
        /// matrix value array
        std::vector<double> a_val;
        /// temporary vector
        std::vector<double> v_berr;

        /// remember whether we already performed a numeric factorization; if so, then we can
        /// reuse the symbolic factorization of the pattern in the next numeric factorization
        bool was_factorized;

      public:
        template<typename IT_>
        explicit Core(const void* comm, Index num_global_dofs, Index my_dof_offset, Index num_owned_dofs,
          const IT_* row_ptr, const IT_* col_idx) :
          _comm(*reinterpret_cast<const MPI_Comm*>(comm)),
          slu_info(0),
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

          // allocate vector for 'berr'
          v_berr.resize(num_owned_dofs);

          // copy pattern
          const Index num_nze(row_ptr[num_owned_dofs]);
          slu_row_ptr.resize(num_owned_dofs+1u);
          for(Index i(0); i <= num_owned_dofs; ++i)
            slu_row_ptr[i] = int_t(row_ptr[i]);
          slu_col_idx.resize(num_nze);
          for(Index i(0); i < num_nze; ++i)
            slu_col_idx[i] = int_t(col_idx[i]);
          slu_col_idx2 = slu_col_idx;

          // allocate matrix value array
          a_val.resize(num_nze);

          // setup SuperLU matrix
          slu_matrix.Stype = SLU_NR_loc; // CSR, local
          slu_matrix.Dtype = SLU_D;      // double
          slu_matrix.Mtype = SLU_GE;     // generic matrix
          slu_matrix.nrow = int_t(num_global_dofs);
          slu_matrix.ncol = int_t(num_global_dofs); // always square
          slu_matrix.Store = &slu_matrix_store;

          // setup SuperLU matrix storage
          slu_matrix_store.nnz_loc = int_t(num_nze);
          slu_matrix_store.m_loc   = int_t(num_owned_dofs);
          slu_matrix_store.fst_row = int_t(my_dof_offset);
          slu_matrix_store.nzval   = a_val.data();
          slu_matrix_store.rowptr  = slu_row_ptr.data();
          slu_matrix_store.colind  = slu_col_idx.data();

          // set up default options
          set_default_options_dist(&slu_opts);

          // don't print statistics to cout
          slu_opts.PrintStat = NO;

          // initialize scale perm structure
          dScalePermstructInit(int_t(num_global_dofs), int_t(num_global_dofs), &slu_scale_perm);

          // initialize LU structure
          dLUstructInit(int_t(num_global_dofs), &slu_lu_struct);

          // initialize the statistics
          PStatInit(&slu_stats);
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

        int init_numeric(const double* vals)
        {
          // copy matrix values
          for(std::size_t i(0); i < a_val.size(); ++i)
            a_val[i] = vals[i];

          // reset column-index array because it might have been
          // overwritten by the previous factorization call
          slu_col_idx = slu_col_idx2;

          // have we already factorized before?
          slu_opts.Fact = (was_factorized ? SamePattern : DOFACT);

          // call the solver routine to factorize the matrix
          pdgssvx(
            &slu_opts,
            &slu_matrix,
            &slu_scale_perm,
            nullptr,
            int(slu_matrix_store.m_loc),
            0,
            &slu_grid,
            &slu_lu_struct,
            &slu_solve_struct,
            v_berr.data(),
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

        int solve(double* x)
        {
          // matrix was already factorized in set_matrix_values() call
          slu_opts.Fact = FACTORED;

          // call the solver routine
          pdgssvx(
            &slu_opts,
            &slu_matrix,
            &slu_scale_perm,
            x,
            int(slu_matrix_store.m_loc),
            1,
            &slu_grid,
            &slu_lu_struct,
            &slu_solve_struct,
            v_berr.data(),
            &slu_stats,
            &slu_info);

          // return the result
          return slu_info;
        }
      }; // class Core

      void* create_core(const void* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, const unsigned int* row_ptr, const unsigned int* col_idx)
      {
        return new Core(comm, num_global_dofs, dof_offset, num_owned_dofs, row_ptr, col_idx);
      }

      void* create_core(const void* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, const unsigned long* row_ptr, const unsigned long* col_idx)
      {
        return new Core(comm, num_global_dofs, dof_offset, num_owned_dofs, row_ptr, col_idx);
      }

      void* create_core(const void* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, const unsigned long long* row_ptr, const unsigned long long* col_idx)
      {
        return new Core(comm, num_global_dofs, dof_offset, num_owned_dofs, row_ptr, col_idx);
      }

      void destroy_core(void* core)
      {
        delete reinterpret_cast<Core*>(core);
      }

      int init_numeric(void* core, const double* vals)
      {
        return reinterpret_cast<Core*>(core)->init_numeric(vals);
      }

      int solve(void* core, double* x)
      {
        return reinterpret_cast<Core*>(core)->solve(x);
      }
    } // namespace SuperLU_Aux
  } // namespace Solver
} // namespace FEAT

#else
// insert dummy function to suppress linker warnings
void dummy_superlu_function() {}
#endif // defined(FEAT_HAVE_SUPERLU_DIST)  && defined(FEAT_HAVE_MPI)
