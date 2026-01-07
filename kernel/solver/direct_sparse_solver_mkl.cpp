// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/solver/direct_sparse_solver.hpp>

#ifdef FEAT_HAVE_MKL

FEAT_DISABLE_WARNINGS
#include <mkl_dss.h>
#include <mkl_cluster_sparse_solver.h>
FEAT_RESTORE_WARNINGS

namespace FEAT
{
  namespace Solver
  {
    namespace DSS
    {
      static_assert(sizeof(MKLDSS_IT) == sizeof(MKL_INT), "DirectSparseSolver: MKL-DSS: index type size mismatch!");

      /**
       * \brief Intel MKL Cluster Sparse Solver core wrapper class
       *
       * \author Peter Zajac
       */
      class MKLDSS_Core
      {
      public:
#ifdef FEAT_HAVE_MPI
        /// MKL-DSS handle array for MPI builds
        int mpi_comm;
        /// MKL-DSS wants an array of 64 void pointers for whatever reason o_O
        void* css_handle[64];
        /// MKL-DSS parameter array
        MKL_INT css_iparm[64];
#else // no MPI
        /// MKL DSS handle for non-MPI builds
        void* dss_handle;
#endif // FEAT_HAVE_MPI
        /// matrix type
        MKL_INT matrix_type;
        /// number of right-hand-sides
        MKL_INT num_rhs;
        /// system data
        MKL_INT num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes;
        /// peak memory during symbolic and numeric factorization
        std::int64_t peak_mem_sym, peak_mem_num;
        /// matrix structure arrays
        std::vector<MKL_INT> row_ptr, col_idx;
        /// matrix and vector data arrays
        std::vector<double> mat_val, rhs_val, sol_val;

        explicit MKLDSS_Core(const Dist::Comm& comm, Index num_global_dofs_, Index dof_offset_,
          Index num_owned_dofs_, Index num_owned_nzes_, Index DOXY(num_global_nzes_)) :
#ifdef FEAT_HAVE_MPI
          mpi_comm(MPI_Comm_c2f(comm.mpi_comm())),
#else
          dss_handle(nullptr),
#endif // FEAT_HAVE_MPI
          matrix_type(11), // real unsymmmetric
          num_rhs(1),
          num_global_dofs(static_cast<MKL_INT>(num_global_dofs_)),
          dof_offset(static_cast<MKL_INT>(dof_offset_)),
          num_owned_dofs(static_cast<MKL_INT>(num_owned_dofs_)),
          num_owned_nzes(static_cast<MKL_INT>(num_owned_nzes_)),
          peak_mem_sym(0),
          peak_mem_num(0)
        {
#ifdef FEAT_HAVE_MPI
          memset(css_handle, 0, sizeof(void*)*64);
          memset(css_iparm, 0, sizeof(MKL_INT)*64);

          // set parameters; see
          // https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/cluster-sparse-solver-iparm-parameter.html
          css_iparm[ 0] =  1; // override default parameters
          css_iparm[ 1] = 10; // use MPI-parallel factorization
          css_iparm[ 4] =  0; // let MKL handle the permutation
          css_iparm[ 5] =  0; // write solution into sol vector
          css_iparm[ 7] =  0; // default number of iterative refinement steps
          css_iparm[ 9] = 13; // perturb small pivots by 10^{-13} (default value)
          css_iparm[10] =  0; // disable diagonal scaling (would require values in symbolic asm)
          css_iparm[11] =  0; // solve A*X=B
          css_iparm[12] =  0; // disable weighted matching (would require values in symbolic asm)
          css_iparm[17] = -1; // enable reporting of non-zero elements in pivots
          css_iparm[26] =  1; // check matrix for errors
          css_iparm[27] =  0; // data type is double precision
          css_iparm[30] =  0; // no partial solve
          css_iparm[34] =  1; // zero-based indexing (yes, '1' REALLY stands for zero-based and '0' for one-based indexing)
          css_iparm[35] =  0; // something about Schur-complements...
          css_iparm[36] =  0; // matrix storage is CSR
          css_iparm[39] =  2; // distributed CSR and distributed vectors
          css_iparm[40] = dof_offset;
          css_iparm[41] = dof_offset + num_owned_dofs - 1;
          css_iparm[59] =  0; // in-core mode
#else // no MPI
          (void)comm;
          // create DSS handle
          MKL_INT opt = MKL_DSS_TERM_LVL_ERROR + MKL_DSS_ZERO_BASED_INDEXING;
#ifdef DEBUG
          opt += MKL_DSS_MSG_LVL_WARNING;
#endif
          MKL_INT ret = ::dss_create(dss_handle, opt);

          switch(ret)
          {
          case MKL_DSS_SUCCESS:
            break;
          case MKL_DSS_OUT_OF_MEMORY:
            throw DirectSparseSolverException("MKL-DSS", "out of memory");
          default:
            throw DirectSparseSolverException("MKL-DSS", "unknown error");
          }
#endif // FEAT_HAVE_MPI

          // allocate matrix and vector arrays
          row_ptr.resize(std::size_t(num_owned_dofs+1), 0u);
          col_idx.resize(std::size_t(num_owned_nzes), 0u);
          mat_val.resize(std::size_t(num_owned_nzes), 0.0);
          rhs_val.resize(std::size_t(num_owned_dofs), 0.0);
          sol_val.resize(std::size_t(num_owned_dofs), 0.0);
        }

        ~MKLDSS_Core()
        {
#ifdef FEAT_HAVE_MPI
          MKL_INT maxfct =  1; // must be set to 1 (is ignored)
          MKL_INT mnum   =  1; // must be set to 1 (is ignored)
          MKL_INT phase  = -1; // cleanup phase
          MKL_INT msglvl =  0; // print statistics
          MKL_INT error  =  0; // error code output variable

          // call solver routine to perform cleanup
          cluster_sparse_solver(
            css_handle,
            &maxfct,
            &mnum,
            &matrix_type,
            &phase,
            &num_global_dofs,
            mat_val.data(),
            row_ptr.data(),
            col_idx.data(),
            nullptr,
            &num_rhs,
            css_iparm,
            &msglvl,
            rhs_val.data(),
            sol_val.data(),
            &mpi_comm,
            &error);
#else // no MPI
          if(dss_handle)
          {
            MKL_INT opt = MKL_DSS_TERM_LVL_ERROR;
#ifdef DEBUG
            opt += MKL_DSS_MSG_LVL_WARNING;
#endif
            ::dss_delete(dss_handle, opt);
            dss_handle = nullptr;
          }
#endif // FEAT_HAVE_MPI
        }

        void init_symbolic()
        {
#ifdef FEAT_HAVE_MPI
          MKL_INT maxfct =  1; // must be set to 1 (is ignored)
          MKL_INT mnum   =  1; // must be set to 1 (is ignored)
          MKL_INT phase  = 11; // symbolical factorization phase
          MKL_INT msglvl =  0; // print statistics
          MKL_INT error  =  0; // error code output variable

          // call factorization routine
          cluster_sparse_solver(
            css_handle,
            &maxfct,
            &mnum,
            &matrix_type,
            &phase,
            &num_global_dofs,
            mat_val.data(),
            row_ptr.data(),
            col_idx.data(),
            nullptr,
            &num_rhs,
            css_iparm,
            &msglvl,
            rhs_val.data(),
            sol_val.data(),
            &mpi_comm,
            &error);

          switch(error)
          {
          case 0: // no error
            return;

          case -2: // out of memory
            throw DirectSparseSolverException("MKL-DSS", "out of memory");

          default:
            throw DirectSparseSolverException("MKL-DSS", "unknown symbolic factorization error");
          }

          // collect statistics
          peak_mem_sym = css_iparm[14] * 1024ll; // size is given in kB
          peak_mem_num = css_iparm[16] * 1024ll; // ditto
#else // no MPI
          // set matrix structure
          MKL_INT opt = MKL_DSS_NON_SYMMETRIC;
          MKL_INT ret = ::dss_define_structure(dss_handle, opt, row_ptr.data(), num_owned_dofs,
            num_owned_dofs, col_idx.data(), num_owned_nzes);

          switch(ret)
          {
          case MKL_DSS_SUCCESS:
            break;
          case MKL_DSS_OUT_OF_MEMORY:
            throw DirectSparseSolverException("MKL-DSS", "out of memory");
          case MKL_DSS_STRUCTURE_ERR:
            throw DirectSparseSolverException("MKL-DSS", "invalid matrix structure");
          default:
            throw DirectSparseSolverException("MKL-DSS", "unknown error");
          }

          // reorder matrix
          opt = MKL_DSS_AUTO_ORDER;
          ret = ::dss_reorder(dss_handle, opt, nullptr);

          switch(ret)
          {
          case MKL_DSS_SUCCESS:
            break;
          case MKL_DSS_OUT_OF_MEMORY:
            throw DirectSparseSolverException("MKL-DSS", "out of memory");
          default:
            throw DirectSparseSolverException("MKL-DSS", "unknown error");
          }

          // collect statistics
          opt = 0;
          double stats[3] = {0.0, 0.0, 0.0};
          dss_statistics(dss_handle, opt, "Peakmem,Factormem,Solvemem", stats);
          peak_mem_sym = std::int64_t(stats[0] * 1024.0);
          peak_mem_num = std::int64_t((stats[1]+stats[2]) * 1024.0);
          //peak_mem_sym
#endif // FEAT_HAVE_MPI
        }

        void init_numeric()
        {
#ifdef FEAT_HAVE_MPI
          MKL_INT maxfct =  1; // must be set to 1 (is ignored)
          MKL_INT mnum   =  1; // must be set to 1 (is ignored)
          MKL_INT phase  = 22; // symbolical factorization phase
          MKL_INT msglvl =  0; // print statistics
          MKL_INT error  =  0; // error code output variable

          // call factorization routine
          cluster_sparse_solver(
            css_handle,
            &maxfct,
            &mnum,
            &matrix_type,
            &phase,
            &num_global_dofs,
            mat_val.data(),
            row_ptr.data(),
            col_idx.data(),
            nullptr,
            &num_rhs,
            css_iparm,
            &msglvl,
            rhs_val.data(),
            sol_val.data(),
            &mpi_comm,
            &error);

          switch(error)
          {
          case 0: // no error
            return;

          case -2: // out of memory
            throw DirectSparseSolverException("MKL-DSS", "out of memory");

          case -4: // zero pivot
            throw DirectSparseSolverException("MKL-DSS", "zero pivot");

          default:
            throw DirectSparseSolverException("MKL-DSS", "unknown numeric factorization error");
          }
#else // no MPI
          MKL_INT opt = MKL_DSS_INDEFINITE;
          MKL_INT ret = ::dss_factor_real(dss_handle, opt, mat_val.data());

          switch(ret)
          {
          case MKL_DSS_SUCCESS:
            break;
          case MKL_DSS_OUT_OF_MEMORY:
            throw DirectSparseSolverException("MKL-DSS", "out of memory");
          case MKL_DSS_ZERO_PIVOT:
            throw DirectSparseSolverException("MKL-DSS", "zero pivot");
          default:
            throw DirectSparseSolverException("MKL-DSS", "unknown error");
          }
#endif // FEAT_HAVE_MPI
        }

        void solve()
        {
#ifdef FEAT_HAVE_MPI
          MKL_INT maxfct =  1; // must be set to 1 (is ignored)
          MKL_INT mnum   =  1; // must be set to 1 (is ignored)
          MKL_INT phase  = 33; // solution phase
          MKL_INT msglvl =  0; // print statistics
          MKL_INT error  =  0; // error code output variable

          // call factorization routine
          cluster_sparse_solver(
            css_handle,
            &maxfct,
            &mnum,
            &matrix_type,
            &phase,
            &num_global_dofs,
            mat_val.data(),
            row_ptr.data(),
            col_idx.data(),
            nullptr,
            &num_rhs,
            css_iparm,
            &msglvl,
            rhs_val.data(),
            sol_val.data(),
            &mpi_comm,
            &error);

          switch(error)
          {
          case 0: // no error
            return;

          case -2: // out of memory
            throw DirectSparseSolverException("MKL-DSS", "out of memory");

          case -4: // zero pivot
            throw DirectSparseSolverException("MKL-DSS", "zero pivot");

          default:
            throw DirectSparseSolverException("MKL-DSS", "unknown solve error");
          }
#else
          MKL_INT opt = 0;
          MKL_INT ret = ::dss_solve_real(dss_handle, opt, rhs_val.data(), num_rhs, sol_val.data());
          switch(ret)
          {
          case MKL_DSS_SUCCESS:
            break;
          default:
            throw DirectSparseSolverException("MKL-DSS", "unknown solve error");
          }
#endif // FEAT_HAVE_MPI
        }

        std::int64_t get_peak_mem() const
        {
          return Math::max(peak_mem_num, peak_mem_sym);
        }
      }; // class MKLDSS_Core

      void* create_mkldss_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes)
      {
        return new MKLDSS_Core(*comm, num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes, num_global_nzes);
      }

      void destroy_mkldss_core(void* core)
      {
        XASSERT(core != nullptr);
        delete reinterpret_cast<MKLDSS_Core*>(core);
      }

      MKLDSS_IT* get_mkldss_row_ptr(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MKLDSS_Core*>(core)->row_ptr.data();
      }

      MKLDSS_IT* get_mkldss_col_idx(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MKLDSS_Core*>(core)->col_idx.data();
      }

      MKLDSS_DT* get_mkldss_mat_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MKLDSS_Core*>(core)->mat_val.data();
      }

      MKLDSS_DT* get_mkldss_rhs_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MKLDSS_Core*>(core)->rhs_val.data();
      }

      MKLDSS_DT* get_mkldss_sol_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MKLDSS_Core*>(core)->sol_val.data();
      }

      void init_mkldss_symbolic(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<MKLDSS_Core*>(core)->init_symbolic();
      }

      void init_mkldss_numeric(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<MKLDSS_Core*>(core)->init_numeric();
      }

      void solve_mkldss(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<MKLDSS_Core*>(core)->solve();
      }

      std::int64_t get_peak_mem_mkldss(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MKLDSS_Core*>(core)->get_peak_mem();
      }
    } // namespace DSS
  } // namespace Solver
} // namespace FEAT

#else // no FEAT_HAVE_MKL

void feat_direct_sparse_solver_mkldss_dummy()
{
}

#endif // FEAT_HAVE_MKL
