// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/solver/direct_sparse_solver.hpp>
#include <kernel/util/omp_util.hpp>

#ifdef FEAT_HAVE_MUMPS

#ifdef FEAT_HAVE_MPI
#include <mpi.h>
#endif // FEAT_HAVE_MPI

#include <dmumps_c.h>

#include <vector>

namespace FEAT
{
  namespace Solver
  {
    namespace DSS
    {
      /**
       * \brief MUMPS core wrapper class
       *
       * \author Peter Zajac
       */
      class  MUMPS_Core
      {
      public:
        const Dist::Comm& comm;
#ifdef FEAT_HAVE_MPI
        /// MPI communicator
        MPI_Comm _comm;
#endif
        /// main MUMPS data structure
        DMUMPS_STRUC_C id;

        const MUMPS_INT num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes, num_global_nzes;
        MUMPS_INT max_owned_dofs;

        /// internal CSR arrays
        std::vector<std::int64_t> row_ptr, col_idx, vec_ptr;
        /// matrix COO arrays for MUMPS in 1-based indices
        std::vector<MUMPS_INT> mumps_row_idx, mumps_col_idx;
        /// matrix value array
        std::vector<double> mat_val;
        /// temporary vectors
        std::vector<double> vector;

        explicit MUMPS_Core(const Dist::Comm& comm_, Index num_global_dofs_, Index dof_offset_,
          Index num_owned_dofs_, Index num_owned_nzes_, Index num_global_nzes_) :
          comm(comm_),
          num_global_dofs(num_global_dofs_),
          dof_offset(dof_offset_),
          num_owned_dofs(num_owned_dofs_),
          num_owned_nzes(num_owned_nzes_),
          num_global_nzes(num_global_nzes_),
          max_owned_dofs(num_owned_dofs_),
          row_ptr(num_owned_dofs_+1),
          col_idx(num_owned_nzes_),
          mumps_row_idx(num_owned_nzes_),
          mumps_col_idx(num_owned_nzes_),
          mat_val(num_owned_nzes_),
          vector()
        {
          // initialize the MUMPS data structure
          memset(&id, 0, sizeof(DMUMPS_STRUC_C));
#ifdef FEAT_HAVE_MPI
          id.comm_fortran = (MUMPS_INT)MPI_Comm_c2f(comm.mpi_comm());
#endif
          id.sym = 0; // unsymmetric matrix
          id.par = 1; // host working

          // initialize MUMPS
          id.job = -1;
          dmumps_c(&id);

          // print errors in debug mode and always keep quiet in opt mode
#ifdef DEBUG
          id.icntl[4-1] = 1; // only error messages
#else
          id.icntl[4-1] = 0; // no messages
#endif

#ifdef FEAT_HAVE_MPI
          std::size_t num_procs = std::size_t(comm.size());
          // allocate vector pointer on rank 0
          if(comm.rank() == 0)
            vec_ptr.resize(num_procs + 1u, 0);

          // gather local vector sizes on rank 0
          std::int64_t loc_size(num_owned_dofs);
          comm.gather(&loc_size, std::size_t(1), vec_ptr.data(), std::size_t(1), 0);

          if(comm.rank() == 0)
          {
            // perform exclusive scan to obtain vector pointer array
            feat_omp_ex_scan(num_procs + 1u, vec_ptr.data(), vec_ptr.data());

            // ensure that the number of global dofs adds up
            XASSERT(vec_ptr.back() == std::int64_t(num_global_dofs));

            // allocate global vector only on rank 0
            vector.resize(num_global_dofs);
          }
          else
          {
            vector.resize(num_owned_dofs);
          }
#else // no FEAT_HAVE_MPI
          vector.resize(num_global_dofs);
#endif // FEAT_HAVE_MPI
        }

        ~MUMPS_Core()
        {
          // release MUMPS
          id.job = -2;
          dmumps_c(&id);
        }

        void init_symbolic()
        {
          // convert 0-based CSR index arrays to 1-based COO index arrays
          FEAT_PRAGMA_OMP(parallel for)
          for(std::int64_t i = 0; i < num_owned_dofs; ++i)
          {
            for(std::int64_t j = row_ptr[i]; j < row_ptr[i+1]; ++j)
            {
              mumps_row_idx[j] = static_cast<MUMPS_INT>(dof_offset + i + 1);
              mumps_col_idx[j] = static_cast<MUMPS_INT>(col_idx[j] + 1);
            }
          }

          // analyze matrix
          id.job = 1;
          // note: the -1 is to account for the 1-based indexing in the MUMPS documentation
          id.icntl[ 5 - 1] = 0; // assembled matrix format
          id.icntl[ 6 - 1] = 0; // column permutation not available for distributed matrices
          id.icntl[18 - 1] = 3; // matrix distributed by user
          id.icntl[20 - 1] = 0; // centratized RHS
          id.icntl[21 - 1] = 0; // centralized solution
          id.n = (MUMPS_INT)num_global_dofs;
          id.nnz_loc = (MUMPS_INT)num_owned_nzes;
          id.irn_loc = mumps_row_idx.data();
          id.jcn_loc = mumps_col_idx.data();
          id.a_loc = mat_val.data();
          id.nrhs = 1;
          id.lrhs = (MUMPS_INT)num_global_dofs;
          id.rhs = vector.data();
          dmumps_c(&id);

          if(id.infog[0] == 0)
            return;

          throw DirectSparseSolverException("mumps",
            "MUMPS Symbolic factorization error with INFO(1) = " + stringify(id.info[0]) + " and INFO(2) = " + stringify(id.info[1]));
        }

        void init_numeric()
        {
          // factorize matrix
          id.job = 2;
          dmumps_c(&id);

          if(id.infog[0] == 0)
            return;

          throw DirectSparseSolverException("mumps",
            "MUMPS Numeric factorization error with INFO(1) = " + stringify(id.info[0]) + " and INFO(2) = " + stringify(id.info[1]));
        }

        void done_numeric()
        {
          // release MUMPS factorization
          id.job = -4;
          dmumps_c(&id);
        }

        void solve()
        {
#ifdef FEAT_HAVE_MPI
          // gather vector on rank 0
          if(comm.size() > 1)
          {
            if(comm.rank() == 0)
            {
              // post receives
              Dist::RequestVector reqs(comm.size());
              for(int i = 1; i < comm.size(); ++i)
                reqs[i] = comm.irecv(&vector[vec_ptr[i]], std::size_t(vec_ptr[i+1] - vec_ptr[i]), i);

              // process all receives
              reqs.wait_all();
            }
            else
            {
              // send our local vector to rank 0
              comm.send(vector.data(), vector.size(), 0);
            }
          }
#endif // FEAT_HAVE_MPI

          // solve system
          id.job = 3;
          dmumps_c(&id);

#ifdef FEAT_HAVE_MPI
          // scatter vector from rank 0
          if(comm.size() > 1)
          {
            if(comm.rank() == 0)
            {
              // post sends
              Dist::RequestVector reqs(comm.size());
              for(int i = 1; i < comm.size(); ++i)
                reqs[i] = comm.isend(&vector[vec_ptr[i]], std::size_t(vec_ptr[i+1] - vec_ptr[i]), i);

              // process all sens
              reqs.wait_all();
            }
            else
            {
              // receive our local vector from rank 0
              comm.recv(vector.data(), vector.size(), 0);
            }
          }
#endif // FEAT_HAVE_MPI
        }
      }; // class MUMPS_Core

      void* create_mumps_core(const Dist::Comm* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, Index num_owned_nzes, Index num_global_nzes)
      {
        return new MUMPS_Core(*comm, num_global_dofs, dof_offset, num_owned_dofs, num_owned_nzes, num_global_nzes);
      }

      void destroy_mumps_core(void* core)
      {
        XASSERT(core != nullptr);
        delete reinterpret_cast<MUMPS_Core*>(core);
      }

      MUMPS_IT* get_mumps_row_ptr(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MUMPS_Core*>(core)->row_ptr.data();
      }

      MUMPS_IT* get_mumps_col_idx(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MUMPS_Core*>(core)->col_idx.data();
      }

      MUMPS_DT* get_mumps_mat_val(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MUMPS_Core*>(core)->mat_val.data();
      }

      MUMPS_DT* get_mumps_vector(void* core)
      {
        XASSERT(core != nullptr);
        return reinterpret_cast<MUMPS_Core*>(core)->vector.data();
      }

      void init_mumps_symbolic(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<MUMPS_Core*>(core)->init_symbolic();
      }

      void init_mumps_numeric(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<MUMPS_Core*>(core)->init_numeric();
      }

      void done_mumps_numeric(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<MUMPS_Core*>(core)->done_numeric();
      }

      void solve_mumps(void* core)
      {
        XASSERT(core != nullptr);
        reinterpret_cast<MUMPS_Core*>(core)->solve();
      }
    } // namespace DSS
  } // namespace Solver
} // namespace FEAT

#else // no FEAT_HAVE_MUMPS

void feat_direct_sparse_solver_mumps_dummy()
{
}

#endif // FEAT_HAVE_MUMPS
