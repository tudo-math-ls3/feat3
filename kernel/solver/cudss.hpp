// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace FEAT
{
  namespace Solver
  {
#if defined(FEAT_HAVE_CUDSS) || defined(DOXYGEN)
    /// native data type of cuDSS solver
    typedef double cudss_native_data_type;
    /// native index type of cuDSS solver
    typedef std::uint32_t cudss_native_index_type;

    /**
     * \brief Nvidia CUDA Direct Sparse Solver wrapper
     *
     * This class provides an implementation of the SolverBase interface using the
     * Nvidia CUDA direct sparse solver "cuDSS" for doing the actual dirty work.
     *
     * \note
     * This solver can only be applied onto SparseMatrixCSR<double,Index> matrices.
     *
     * \attention
     * This class is only declared if FEAT was configured to build and link against the \c Nvidia cuDSS library.
     *
     * \author Peter Zajac
     */
    class CUDSS :
      public SolverBase<LAFEM::DenseVector<double, Index>>
    {
    public:
      /// compatible matrix type
      typedef LAFEM::SparseMatrixCSR<double, Index> MatrixType;
      /// compatible vector type
      typedef LAFEM::DenseVector<double, Index> VectorType;

      /// our base class
      typedef SolverBase<VectorType> BaseClass;

    private:
      /// system matrix
      const MatrixType& _system_matrix;
      void* _cudss_core;

    public:
      /**
      * \brief Constructor
      *
      * \param[in] system_matrix
      * A reference to the system matrix to be factorized.
      */
      explicit CUDSS(const MatrixType& system_matrix);

      /// virtual destructor
      virtual ~CUDSS();

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "CUDSS";
      }

      virtual void init_symbolic() override;
      virtual void done_symbolic() override;
      virtual void init_numeric() override;
      virtual void done_numeric() override;

      /**
      * \brief Solves a linear system with the factorized system matrix.
      *
      * \param[in,out] vec_sol
      * A reference to the solution vector. The vector must be allocated to the correct length, but its
      * initial contents are ignored.
      *
      * \param[in] vec_rhs
      * A reference to the right-hand-side of the linear system.
      */
      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override;
    }; // class CUDSS

    /**
    * \brief Creates a new CUDSS solver object
    *
    * \param[in] matrix
    * The system matrix.
    *
    * \returns
    * A shared pointer to a new CUDSS object.
    */
    inline std::shared_ptr<CUDSS> new_cudss(const LAFEM::SparseMatrixCSR<double, Index>& matrix)
    {
      return std::make_shared<CUDSS>(matrix);
    }
#endif // defined(FEAT_HAVE_CUDSS) || defined(DOXYGEN)
  } // namespace Solver
} // namespace FEAT
