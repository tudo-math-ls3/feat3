// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_SUPERLU_HPP
#define KERNEL_SOLVER_SUPERLU_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/adp_solver_base.hpp>

#if defined(FEAT_HAVE_SUPERLU_DIST) || defined(DOXYGEN)

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    /**
     * \brief SuperLU wrapper function namespace
     *
     * This namespace encapsulates our opaque wrapper functions for SuperLU.
     */
    namespace SuperLU_Aux
    {
      /**
       * \brief Creates an opaque core wrapper object for SuperLU.
       *
       * The "core wrapper object", which is returned by this function, is basically a collection
       * of all SuperLU data structures, which are required to represent a linear system A*x=b, i.e.
       * a partitioned matrix as well as two vectors.
       *
       * \param[in] comm
       * A pointer to an MPI_Comm object that represents the communicator.
       *
       * \param[in] num_global_dofs
       * The total number of global DOFs.
       *
       * \param[in] dof_offset
       * The global DOF offset for this process.
       *
       * \param[in] num_owned_dofs
       * The number of global DOFs owned by this process.
       *
       * \param[in] row_ptr
       * The row-pointer array of the partitioned CSR matrix.
       *
       * \param[in] col_idx
       * The column-index array of the partitioned CSR matrix.
       *
       * \returns
       * A pointer to a newly allocated core wrapper object.
       */
      void* create_core(const void* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, const unsigned int* row_ptr, const unsigned int* col_idx);
      void* create_core(const void* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, const unsigned long* row_ptr, const unsigned long* col_idx);
      void* create_core(const void* comm, Index num_global_dofs, Index dof_offset,
        Index num_owned_dofs, const unsigned long long* row_ptr, const unsigned long long* col_idx);

      /// destroys a wrapper core
      void destroy_core(void* core);

      /// sets the matrix values and performs numeric factorization
      int init_numeric(void* core, const double* vals);

      /**
       * \brief Solves a linear system
       *
       * \param[inout] x
       * On entry, the right-hand side vector.
       * On exit, the solution vector.
       */
      int solve(void* core, double* x);
    } // namespace SuperLU_Aux
    /// \endcond

    /**
     * \brief (distributed) SuperLU direct sparse solver
     *
     * This class implements a wrapper around the distributed SuperLU direct sparse solver,
     * which is provided by the SuperLU third-party library.
     *
     * \attention
     * Because this class only wraps around the SuperLU_DIST library, it is only available when
     * FEAT is configured and linked with MPI support, i.e. there exists no non-MPI version of
     * this SuperLU class in serial builds. If you need a serial direct solver, look for UMFPACK.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class SuperLU :
      public ADPSolverBase<Matrix_, Filter_>
    {
    public:
      /// our base class
      typedef ADPSolverBase<Matrix_, Filter_> BaseClass;
      /// our vector type
      typedef typename Matrix_::VectorTypeL VectorType;

    protected:
      /// a pointer to our opaque SuperLU core object
      void* _core;

    public:
      /**
       * \brief Constructor
       */
      explicit SuperLU(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _core(nullptr)
      {
      }

      virtual String name() const override
      {
        return "SuperLU";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        XASSERT(this->_core == nullptr);

        // create our core wrapper object
        this->_core = SuperLU_Aux::create_core(
#ifdef FEAT_HAVE_MPI
          &this->_get_comm()->mpi_comm(),
#else
          nullptr, // no communicator for non-MPI builds
#endif
          this->_get_num_global_dofs(),
          this->_get_global_dof_offset(),
          this->_get_num_owned_dofs(),
          this->_get_mat_row_ptr(),
          this->_get_mat_col_idx());
      }

      virtual void done_symbolic() override
      {
        XASSERT(this->_core != nullptr);

        SuperLU_Aux::destroy_core(this->_core);
        this->_core = nullptr;

        BaseClass::done_symbolic();
      }

      virtual void init_numeric() override
      {
        XASSERT(this->_core != nullptr);
        BaseClass::init_numeric();

        // set the new matrix values
        const int info = SuperLU_Aux::init_numeric(this->_core, this->_get_mat_vals());

        // info = 0 means success
        // info < 0 means illegal argument; this should never happen
        XASSERTM(info >= 0, "SuperLU: illegal argument in init_numeric()");
        if(0 < info)
        {
          // 0 < info < N: zero pivot = singular matrix
          if(info <= int(this->_get_num_global_dofs()))
            throw SingularMatrixException("SuperLU: singular matrix");
          else // info > N: out of memory
            throw SolverException("SuperLU: out of memory");
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector into internal correction vector
        this->_upload_vec_cor(vec_def);

        // solve
        int info = SuperLU_Aux::solve(this->_core, this->_get_vec_cor_vals(vec_cor));

        // info != 0 means something went wrong, but this should never happen here...
        XASSERTM(info == 0, "SuperLU: solve returned non-zero info value");

        // download solution into correction vector
        this->_download_vec_cor(vec_cor);

        // okay
        return Status::success;
      }
    }; // class SuperLU

    /**
     * \brief Creates a new SuperLU solver object
     *
     * \param[in] matrix
     * The system matrix
     *
     * \param[in] filter
     * The system filter
     *
     * \returns
     * A shared pointer to a new SuperLU object.
     */
    template<typename Matrix_, typename Filter_>
    std::shared_ptr<SuperLU<Matrix_, Filter_>> new_superlu(const Matrix_& matrix, const Filter_& filter)
    {
#ifdef FEAT_HAVE_MPI
      return std::make_shared<SuperLU<Matrix_, Filter_>>(matrix, filter);
#else
      XABORT("new_superlu: SuperLU solver is only available in MPI builds")
#endif
    }
  } // namespace Solver
} // namespace FEAT

#endif // defined(FEAT_HAVE_SUPERLU_DIST) || defined(DOXYGEN)
#endif // KERNEL_SOLVER_SUPERLU_HPP
