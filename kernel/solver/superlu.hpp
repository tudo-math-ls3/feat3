// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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
      /// the data type used by SuperLU
      typedef double SLUDataType;

      /// the index type used by SuperLU
      /// must be identical to 'int_t' in <superlu_defs.h>
#ifdef FEAT_TPL_SUPERLU_INT64
      typedef std::int64_t SLUIndexType;
#else
      typedef int SLUIndexType;
#endif

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
       * \param[in] num_nonzeros
       * The number of non-zero entries for this process
       *
       * \returns
       * A pointer to a newly allocated core wrapper object.
       */
      void* create_core(const void* comm, Index num_global_dofs, Index dof_offset, Index num_owned_dofs, Index num_nonzeros);

      /// destroys a wrapper core
      void destroy_core(void* core);

      SLUIndexType* get_row_ptr(void* core);
      SLUIndexType* get_col_idx(void* core);
      SLUDataType* get_mat_val(void* core);
      SLUDataType* get_vector(void* core);

      void init_symbolic(void* core);
      int init_numeric(void* core);

      /**
       * \brief Solves a linear system
       */
      int solve(void* core);
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
          this->_get_adp_matrix_num_nzes());

        // upload matrix structure
        BaseClass::_upload_symbolic(SuperLU_Aux::get_row_ptr(this->_core), SuperLU_Aux::get_col_idx(this->_core));

        // initialize
        SuperLU_Aux::init_symbolic(this->_core);
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

        // upload matrix values
        this->_upload_numeric(SuperLU_Aux::get_mat_val(this->_core),
          SuperLU_Aux::get_row_ptr(this->_core), SuperLU_Aux::get_col_idx(this->_core));

        // set the new matrix values
        const int info = SuperLU_Aux::init_numeric(this->_core);

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
        this->_upload_vector(SuperLU_Aux::get_vector(this->_core), vec_def.local());

        // solve
        int info = SuperLU_Aux::solve(this->_core);

        // info != 0 means something went wrong, but this should never happen here...
        XASSERTM(info == 0, "SuperLU: solve returned non-zero info value");

        // download solution into correction vector
        this->_download_vector(vec_cor.local(), SuperLU_Aux::get_vector(this->_core));

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
