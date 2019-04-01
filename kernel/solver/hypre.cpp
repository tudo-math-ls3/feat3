// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>

#ifdef FEAT_HAVE_HYPRE

// preprocessor macro sanity checks
#ifdef FEAT_HAVE_MPI
#  ifdef HYPRE_SEQUENTIAL
#    error FEAT_HAVE_MPI and HYPRE_SEQUENTIAL are both defined
#  endif
#  ifndef HYPRE_HAVE_MPI
#    define HYPRE_HAVE_MPI 1
#  endif
#else
#  ifdef HYPRE_HAVE_MPI
#    error FEAT_HAVE_MPI is not defined, but HYPRE_HAVE_MPI is defined
#  endif
#  ifndef HYPRE_SEQUENTIAL
#    define HYPRE_SEQUENTIAL 1
#  endif
#endif // FEAT_HAVE_MPI

// includes, HYPRE third-party lib
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Solver
  {
    namespace Hypre
    {
      /**
       * \brief HYPRE wrapper core data class
       */
      class Core
      {
      public:
        // Note: if compiled without MPI, HYPRE defines its own MPI_Comm dummy type
        MPI_Comm _comm;
        /// HYPRE global DOF indices array
        std::vector<HYPRE_Int> _hypre_dof_idx;
        /// HYPRE number of non-zero entries per row array
        std::vector<HYPRE_Int> _hypre_num_nze;
        /// HYPRE column indices array
        std::vector<HYPRE_Int> _hypre_col_idx;

        // note: the following members are actually pointers
        /// handle to the HYPRE IJ matrix object
        HYPRE_IJMatrix _hypre_matrix;
        /// handles to the HYPRE IJ vector objects
        HYPRE_IJVector _hypre_vec_def, _hypre_vec_cor;
        /// HYPRE ParCSRMatrix handle to the system matrix
        HYPRE_ParCSRMatrix _hypre_parcsr;
        /// HYPRE ParVector handles to the defect/correction vectors
        HYPRE_ParVector _hypre_par_def, _hypre_par_cor;


        template<typename IT_>
        explicit Core(const void* comm, Index my_dof_offset, Index num_owned_dofs, const IT_* row_ptr, const IT_* col_idx) :
#ifdef FEAT_HAVE_MPI
          _comm(*reinterpret_cast<const MPI_Comm*>(comm)),
#endif
          _hypre_dof_idx(),
          _hypre_num_nze(),
          _hypre_col_idx(),
          _hypre_matrix(nullptr),
          _hypre_vec_def(nullptr),
          _hypre_vec_cor(nullptr),
          _hypre_parcsr(nullptr),
          _hypre_par_def(nullptr),
          _hypre_par_cor(nullptr)
        {
#ifndef FEAT_HAVE_MPI
          (void)comm; // suppress unused parameter warnings
#endif
          // get dof offset and count
          const HYPRE_Int dof_offset = HYPRE_Int(my_dof_offset);
          const HYPRE_Int num_owned  = HYPRE_Int(num_owned_dofs);

          // set up HYPRE dof indices vector
          _hypre_dof_idx.resize((std::size_t)num_owned);
          for(HYPRE_Int i(0); i < num_owned; ++i)
            _hypre_dof_idx[std::size_t(i)] = dof_offset + i;

          // get lower and upper DOF bounds
          const HYPRE_Int ilower = _hypre_dof_idx.front();
          const HYPRE_Int iupper = _hypre_dof_idx.back();

          // create HYPRE matrix
          HYPRE_IJMatrixCreate(_comm, ilower, iupper, ilower, iupper, &_hypre_matrix);
          HYPRE_IJMatrixSetObjectType(_hypre_matrix, HYPRE_PARCSR);
          HYPRE_IJMatrixInitialize(_hypre_matrix);

          // create HYPRE vectors
          HYPRE_IJVectorCreate(_comm, ilower, iupper, &_hypre_vec_def);
          HYPRE_IJVectorCreate(_comm, ilower, iupper, &_hypre_vec_cor);
          HYPRE_IJVectorSetObjectType(_hypre_vec_def, HYPRE_PARCSR);
          HYPRE_IJVectorSetObjectType(_hypre_vec_cor, HYPRE_PARCSR);
          HYPRE_IJVectorInitialize(_hypre_vec_def);
          HYPRE_IJVectorInitialize(_hypre_vec_cor);

          // get matrix dimensions
          const HYPRE_Int num_nze = HYPRE_Int(row_ptr[num_owned_dofs]);

          // create HYPRE index vectors
          _hypre_num_nze.resize((std::size_t)num_owned);
          _hypre_col_idx.resize((std::size_t)num_nze);

          // loop over all owned matrix rows
          for(HYPRE_Int i(0); i < num_owned; ++i)
          {
            _hypre_num_nze[std::size_t(i)] = HYPRE_Int(row_ptr[i+1] - row_ptr[i]);
          }
          for(HYPRE_Int i(0); i < num_nze; ++i)
          {
            _hypre_col_idx[std::size_t(i)] = HYPRE_Int(col_idx[i]);
          }

          // vector for initial vector values
          std::vector<HYPRE_Real> va((std::size_t)num_nze, 0.0);
          std::vector<HYPRE_Real> vx((std::size_t)num_owned, 0.0);

          // set matrix structure + dummy values
          HYPRE_IJMatrixSetValues(_hypre_matrix, num_owned, _hypre_num_nze.data(), _hypre_dof_idx.data(),
            _hypre_col_idx.data(), va.data());

          // set vector structures + dummy values
          HYPRE_IJVectorSetValues(_hypre_vec_def, num_owned, _hypre_dof_idx.data(), vx.data());
          HYPRE_IJVectorSetValues(_hypre_vec_cor, num_owned, _hypre_dof_idx.data(), vx.data());

          // assemble matrix and get ParCSR object pointer
          HYPRE_IJMatrixAssemble(_hypre_matrix);
          HYPRE_IJMatrixGetObject(_hypre_matrix, (void**) &_hypre_parcsr);

          // assemble vectors and get Par object pointers
          HYPRE_IJVectorAssemble(_hypre_vec_def);
          HYPRE_IJVectorAssemble(_hypre_vec_cor);
          HYPRE_IJVectorGetObject(_hypre_vec_def, (void **) &_hypre_par_def);
          HYPRE_IJVectorGetObject(_hypre_vec_cor, (void **) &_hypre_par_cor);
        }

        ~Core()
        {
          HYPRE_IJVectorDestroy(_hypre_vec_cor);
          HYPRE_IJVectorDestroy(_hypre_vec_def);
          HYPRE_IJMatrixDestroy(_hypre_matrix);
        }

        void set_matrix_values(const double* vals)
        {
          HYPRE_IJMatrixSetValues(_hypre_matrix, HYPRE_Int(_hypre_dof_idx.size()), _hypre_num_nze.data(),
            _hypre_dof_idx.data(), _hypre_col_idx.data(), vals);
        }

        void set_vec_cor_values(const double* vals)
        {
          if(vals != nullptr)
            HYPRE_IJVectorSetValues(_hypre_vec_cor, HYPRE_Int(_hypre_dof_idx.size()), _hypre_dof_idx.data(), vals);
          else
            HYPRE_ParVectorSetConstantValues(_hypre_par_cor, HYPRE_Real(0));
        }

        void set_vec_def_values(const double* vals)
        {
          if(vals != nullptr)
            HYPRE_IJVectorSetValues(_hypre_vec_def, HYPRE_Int(_hypre_dof_idx.size()), _hypre_dof_idx.data(), vals);
          else
            HYPRE_ParVectorSetConstantValues(_hypre_par_def, HYPRE_Real(0));
        }

        void get_vec_cor_values(double* vals) const
        {
          HYPRE_IJVectorGetValues(_hypre_vec_cor, HYPRE_Int(_hypre_dof_idx.size()), _hypre_dof_idx.data(), vals);
        }

        void get_vec_def_values(double* vals) const
        {
          HYPRE_IJVectorGetValues(_hypre_vec_cor, HYPRE_Int(_hypre_dof_idx.size()), _hypre_dof_idx.data(), vals);
        }
      }; // class Core

      void* create_core(const void* comm, Index dof_offset, Index num_owned_dofs,
        const unsigned int* row_ptr, const unsigned int* col_idx)
      {
        return new Core(comm, dof_offset, num_owned_dofs, row_ptr, col_idx);
      }

      void* create_core(const void* comm, Index dof_offset, Index num_owned_dofs,
        const unsigned long* row_ptr, const unsigned long* col_idx)
      {
        return new Core(comm, dof_offset, num_owned_dofs, row_ptr, col_idx);
      }

      void* create_core(const void* comm, Index dof_offset, Index num_owned_dofs,
        const unsigned long long* row_ptr, const unsigned long long* col_idx)
      {
        return new Core(comm, dof_offset, num_owned_dofs, row_ptr, col_idx);
      }

      void destroy_core(void* core)
      {
        delete reinterpret_cast<Core*>(core);
      }

      void set_matrix_values(void* core, const double* vals)
      {
        reinterpret_cast<Core*>(core)->set_matrix_values(vals);
      }

      void set_vec_cor_values(void* core, const double* vals)
      {
        reinterpret_cast<Core*>(core)->set_vec_cor_values(vals);
      }

      void set_vec_def_values(void* core, const double* vals)
      {
        reinterpret_cast<Core*>(core)->set_vec_def_values(vals);
      }

      void get_vec_cor_values(const void* core, double* vals)
      {
        reinterpret_cast<const Core*>(core)->get_vec_cor_values(vals);
      }

      void get_vec_def_values(const void* core, double* vals)
      {
        reinterpret_cast<const Core*>(core)->get_vec_def_values(vals);
      }

      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */

      void* create_parasails(void* cr, int* iparam, double* dparam)
      {
        Core* core = reinterpret_cast<Core*>(cr);
        HYPRE_Solver* solver = new HYPRE_Solver();

        // create ParaSails preconditioner
        HYPRE_ParaSailsCreate(core->_comm, solver);

        // set ParaSails parameters
        HYPRE_ParaSailsSetParams(*solver, dparam[0], iparam[0]);
        HYPRE_ParaSailsSetFilter(*solver, dparam[1]);
        HYPRE_ParaSailsSetSym(*solver, iparam[1]);

        // setup ParaSails
        // according to the documentation, the two vectors are ignored by this function
        HYPRE_ParaSailsSetup(*solver, core->_hypre_parcsr, core->_hypre_par_def, core->_hypre_par_cor);

        return solver;
      }

      void destroy_parasails(void* slv)
      {
        HYPRE_Solver* solver = reinterpret_cast<HYPRE_Solver*>(slv);
        HYPRE_ParaSailsDestroy(*solver);
        delete solver;
      }

      void solve_parasails(void* cr, void* slv)
      {
        Core* core = reinterpret_cast<Core*>(cr);
        HYPRE_Solver* solver = reinterpret_cast<HYPRE_Solver*>(slv);

        HYPRE_ParaSailsSolve(*solver, core->_hypre_parcsr, core->_hypre_par_def, core->_hypre_par_cor);
      }

      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */

      void* create_euclid(void* cr, int* iparam, double* dparam)
      {
        Core* core = reinterpret_cast<Core*>(cr);
        HYPRE_Solver* solver = new HYPRE_Solver();

        // create Euclid preconditioner
        HYPRE_EuclidCreate(core->_comm, solver);

        // set Euclid parameters
        HYPRE_EuclidSetLevel(*solver, iparam[0]);
        HYPRE_EuclidSetSparseA(*solver, dparam[0]);

        // setup Euclid
        // according to the documentation, the two vectors are ignored by this function
        HYPRE_EuclidSetup(*solver, core->_hypre_parcsr, core->_hypre_par_def, core->_hypre_par_cor);

        return solver;
      }

      void destroy_euclid(void* slv)
      {
        HYPRE_Solver* solver = reinterpret_cast<HYPRE_Solver*>(slv);
        HYPRE_EuclidDestroy(*solver);
        delete solver;
      }

      void solve_euclid(void* cr, void* slv)
      {
        Core* core = reinterpret_cast<Core*>(cr);
        HYPRE_Solver* solver = reinterpret_cast<HYPRE_Solver*>(slv);

        HYPRE_EuclidSolve(*solver, core->_hypre_parcsr, core->_hypre_par_def, core->_hypre_par_cor);
      }

      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */

      void* create_boomeramg(void* cr, int* iparam, double* dparam)
      {
        Core* core = reinterpret_cast<Core*>(cr);
        HYPRE_Solver* solver = new HYPRE_Solver();

        // create BoomerAMG preconditioner
        HYPRE_BoomerAMGCreate(solver);

        // set BoomerAMG parameters
        HYPRE_BoomerAMGSetMaxIter(*solver, iparam[0]);
        HYPRE_BoomerAMGSetTol(*solver, dparam[0]);

        // setup BoomerAMG
        // according to the documentation, the two vectors are ignored by this function
        HYPRE_BoomerAMGSetup(*solver, core->_hypre_parcsr, core->_hypre_par_def, core->_hypre_par_cor);

        return solver;
      }

      void destroy_boomeramg(void* slv)
      {
        HYPRE_Solver* solver = reinterpret_cast<HYPRE_Solver*>(slv);
        HYPRE_BoomerAMGDestroy(*solver);
        delete solver;
      }

      void solve_boomeramg(void* cr, void* slv)
      {
        Core* core = reinterpret_cast<Core*>(cr);
        HYPRE_Solver* solver = reinterpret_cast<HYPRE_Solver*>(slv);

        HYPRE_BoomerAMGSolve(*solver, core->_hypre_parcsr, core->_hypre_par_def, core->_hypre_par_cor);
      }
    } // namespace Hypre
  } // namespace Solver
} // namespace FEAT
#else
// insert dummy function to suppress linker warnings
void dummy_hypre_function() {}
#endif // FEAT_HAVE_HYPRE
