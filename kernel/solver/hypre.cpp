// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
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
        const HYPRE_Int _dof_offset;
        const HYPRE_Int _num_owned_dofs;
        const HYPRE_Int _num_nonzeros;
        /// HYPRE global DOF indices array
        std::vector<HYPRE_Int> _hypre_dof_idx;
        /// HYPRE number of non-zero entries per row array
        std::vector<HYPRE_Int> _hypre_num_nze;
        /// HYPRE row pointer array
        std::vector<HYPRE_Int> _hypre_row_ptr;
        /// HYPRE column indices array
        std::vector<HYPRE_Int> _hypre_col_idx;
        /// HYPRE matrix values
        std::vector<HYPRE_Real> _hypre_mat_val;
        /// HYPRE defect vector values
        std::vector<HYPRE_Real> _hypre_def_val;
        /// HYPRE correction vector values
        std::vector<HYPRE_Real> _hypre_cor_val;

        // note: the following members are actually pointers
        /// handle to the HYPRE IJ matrix object
        HYPRE_IJMatrix _hypre_matrix;
        /// handles to the HYPRE IJ vector objects
        HYPRE_IJVector _hypre_vec_def, _hypre_vec_cor;
        /// HYPRE ParCSRMatrix handle to the system matrix
        HYPRE_ParCSRMatrix _hypre_parcsr;
        /// HYPRE ParVector handles to the defect/correction vectors
        HYPRE_ParVector _hypre_par_def, _hypre_par_cor;

        explicit Core(const void* comm, Index my_dof_offset, Index num_owned_dofs, Index num_nonzeros) :
#ifdef FEAT_HAVE_MPI
          _comm(*reinterpret_cast<const MPI_Comm*>(comm)),
#endif
          _dof_offset(HYPRE_Int(my_dof_offset)),
          _num_owned_dofs(HYPRE_Int(num_owned_dofs)),
          _num_nonzeros(HYPRE_Int(num_nonzeros)),
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
          // set up HYPRE dof indices vector
          _hypre_dof_idx.resize((std::size_t)_num_owned_dofs);
          for(HYPRE_Int i(0); i < _num_owned_dofs; ++i)
            _hypre_dof_idx[std::size_t(i)] = _dof_offset + i;

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

          // create HYPRE index vectors
          _hypre_num_nze.resize((std::size_t)_num_owned_dofs);
          _hypre_row_ptr.resize((std::size_t)(_num_owned_dofs+1u));
          _hypre_col_idx.resize((std::size_t)_num_nonzeros);
          _hypre_mat_val.resize((std::size_t)_num_nonzeros);
          _hypre_def_val.resize((std::size_t)_num_owned_dofs);
          _hypre_cor_val.resize((std::size_t)_num_owned_dofs);
        }

        ~Core()
        {
          HYPRE_IJVectorDestroy(_hypre_vec_cor);
          HYPRE_IJVectorDestroy(_hypre_vec_def);
          HYPRE_IJMatrixDestroy(_hypre_matrix);
        }

        void upload_symbolic()
        {
          // compute row non-zeros from row pointer
          std::size_t n = _hypre_num_nze.size();
          for(std::size_t i(0); i < n; ++i)
            _hypre_num_nze[i] = _hypre_row_ptr[i+1] - _hypre_row_ptr[i];

          // set matrix structure + dummy values
          HYPRE_IJMatrixSetValues(_hypre_matrix, _num_owned_dofs, _hypre_num_nze.data(), _hypre_dof_idx.data(),
            _hypre_col_idx.data(), _hypre_mat_val.data());

          // set vector structures + dummy values
          HYPRE_IJVectorSetValues(_hypre_vec_def, _num_owned_dofs, _hypre_dof_idx.data(), _hypre_def_val.data());
          HYPRE_IJVectorSetValues(_hypre_vec_cor, _num_owned_dofs, _hypre_dof_idx.data(), _hypre_cor_val.data());

          // assemble matrix and get ParCSR object pointer
          HYPRE_IJMatrixAssemble(_hypre_matrix);
          HYPRE_IJMatrixGetObject(_hypre_matrix, (void**) &_hypre_parcsr);

          // assemble vectors and get Par object pointers
          HYPRE_IJVectorAssemble(_hypre_vec_def);
          HYPRE_IJVectorAssemble(_hypre_vec_cor);
          HYPRE_IJVectorGetObject(_hypre_vec_def, (void **) &_hypre_par_def);
          HYPRE_IJVectorGetObject(_hypre_vec_cor, (void **) &_hypre_par_cor);
        }

        void upload_mat_val()
        {
          HYPRE_IJMatrixSetValues(_hypre_matrix, _num_owned_dofs, _hypre_num_nze.data(), _hypre_dof_idx.data(),
            _hypre_col_idx.data(), _hypre_mat_val.data());
        }

        void upload_vec_def()
        {
          HYPRE_IJVectorSetValues(_hypre_vec_def, _num_owned_dofs, _hypre_dof_idx.data(), _hypre_def_val.data());
        }

        void download_vec_cor()
        {
          HYPRE_IJVectorGetValues(_hypre_vec_cor, _num_owned_dofs, _hypre_dof_idx.data(), _hypre_cor_val.data());
        }

        void format_vec_cor()
        {
          HYPRE_ParVectorSetConstantValues(_hypre_par_cor, HYPRE_Real(0));
        }

        HYPRE_Int* get_row_ptr()
        {
          return _hypre_row_ptr.data();
        }

        HYPRE_Int* get_col_idx()
        {
          return _hypre_col_idx.data();
        }

        HYPRE_Real* get_mat_val()
        {
          return _hypre_mat_val.data();
        }

        HYPRE_Real* get_vec_def()
        {
          return _hypre_def_val.data();
        }

        HYPRE_Real* get_vec_cor()
        {
          return _hypre_cor_val.data();
        }
      }; // class Core

      void* create_core(const void* comm, Index dof_offset, Index num_owned_dofs, Index num_nonzeros)
      {
        return new Core(comm, dof_offset, num_owned_dofs, num_nonzeros);
      }

      void destroy_core(void* core)
      {
        delete reinterpret_cast<Core*>(core);
      }

      HYPRE_Int* get_row_ptr(void* core)
      {
        return reinterpret_cast<Core*>(core)->get_row_ptr();
      }

      HYPRE_Int* get_col_idx(void* core)
      {
        return reinterpret_cast<Core*>(core)->get_col_idx();
      }

      HYPRE_Real* get_mat_val(void* core)
      {
        return reinterpret_cast<Core*>(core)->get_mat_val();
      }

      HYPRE_Real* get_vec_def(void* core)
      {
        return reinterpret_cast<Core*>(core)->get_vec_def();
      }

      HYPRE_Real* get_vec_cor(void* core)
      {
        return reinterpret_cast<Core*>(core)->get_vec_cor();
      }

      void upload_symbolic(void* core)
      {
        reinterpret_cast<Core*>(core)->upload_symbolic();
      }

      void upload_mat_val(void* core)
      {
        reinterpret_cast<Core*>(core)->upload_mat_val();
      }

      void upload_vec_def(void* core)
      {
        reinterpret_cast<Core*>(core)->upload_vec_def();
      }

      void download_vec_cor(void* core)
      {
        reinterpret_cast<Core*>(core)->download_vec_cor();
      }

      void format_vec_cor(void* core)
      {
        reinterpret_cast<Core*>(core)->format_vec_cor();
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
