// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Symmetric lumped Schur complement matrix
     *
     * \tparam LumpedMatrixA_
     * Type for the vector representing the diagonal middle matrix
     *
     * \tparam MatrixB_
     * Type for the right matrix
     *
     * \tparam MatrixD_
     * Type for the left matrix
     *
     * \tparam FilterA_
     * Filter for the diagonal middle matrix
     *
     * For given matrices \f$ B, A = diag(a), D = B^T \f$, this emulates the matrix
     * \f[ S = B^T A^{-1} B. \f]
     * This is called Lumped because the diagonal structure of A usually comes from lumping.
     *
     * \warning For efficiency reasons, \f$ B^T \f$ needs to be explicitly computed and passed as \f$ D \f$, but this
     * property is only checked very superficially (and only in debug mode).
     *
     * \note In the current implementation, all container template parameters must be Global:: as the synchronization
     * is done manually by using local() objects and calling sync_0()
     *
     * \author Jordi Paul
     */
    template
    <
      typename LumpedMatrixA_,
      typename MatrixB_,
      typename MatrixD_,
      typename FilterA_
    >
    class SymmetricLumpedSchurMatrix
    {
    public:
      /// The type of A = diag(a)
      typedef LumpedMatrixA_ LumpedMatrixTypeA;
      /// The type of B
      typedef MatrixB_ MatrixTypeB;
      /// The type of D = B^T
      typedef MatrixD_ MatrixTypeD;
      /// The filter for A
      typedef FilterA_ FilterTypeA;

      /// The floating point precision
      typedef typename LumpedMatrixA_::DataType DataType;
      /// The vector for left-multiplication with S
      typedef typename MatrixD_::VectorTypeL VectorTypeL;
      /// The vector for right-multiplication with S
      typedef typename MatrixB_::VectorTypeR VectorTypeR;

      /// The left-vector type for A
      typedef LumpedMatrixA_ VectorTypeML;
      /// The right-vector type for A
      typedef LumpedMatrixA_ VectorTypeMR;

      /// The row-gate type (used by SFINAE)
      typedef typename MatrixD_::GateRowType GateRowType;
      /// The column-gate type (used by SFINAE)
      typedef typename MatrixB_::GateColType GateColType;

      static constexpr bool is_global = true;
      static constexpr bool is_local = false;

      /// A = diag(a)
      LumpedMatrixA_ inv_lumped_matrix_a;
      /// B
      const MatrixB_& matrix_b;
      /// D = B^T
      const MatrixD_& matrix_d;
      /// The filter for a
      const FilterA_& filter_a;

    private:
      // These two need to be modified even when apply() (which is const) is called
      /// Buffer vector for receiving v <- x^T A
      mutable VectorTypeML _vec_ml;
      /// Buffer vector for receiving v <- A x
      mutable VectorTypeMR _vec_mr;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] lumped_matrix_a_
       * A = diag(a)
       *
       * \param[in] matrix_b_
       * B
       *
       * \param[in] matrix_d_
       * D = B^T
       *
       * \param[in] filter_a_
       * Filter for A
       */
      SymmetricLumpedSchurMatrix(
        const LumpedMatrixA_& lumped_matrix_a_,
        const MatrixB_& matrix_b_,
        const MatrixD_& matrix_d_,
        const FilterA_& filter_a_) :
        inv_lumped_matrix_a(lumped_matrix_a_.clone(LAFEM::CloneMode::Layout)),
        matrix_b(matrix_b_),
        matrix_d(matrix_d_),
        filter_a(filter_a_),
        _vec_ml(lumped_matrix_a_.clone(LAFEM::CloneMode::Layout)),
        _vec_mr(lumped_matrix_a_.clone(LAFEM::CloneMode::Layout))
      {
        ASSERT(matrix_d.num_cols() == matrix_b.num_rows());
        ASSERT(matrix_d.num_nzes() == matrix_b.num_nzes());

        inv_lumped_matrix_a.component_invert(lumped_matrix_a_);
      }

      /**
       * \brief Virtual destructor
       */
      virtual ~SymmetricLumpedSchurMatrix()
      {
      }

      const Dist::Comm* get_comm() const
      {
        return inv_lumped_matrix_a.get_comm();
      }

      void update_lumped_a(const LumpedMatrixA_& lumped_matrix_a_)
      {
        // If these were initialized empty (as it frequently happens with Global containers), adjust the sizes
        if(_vec_ml.local().size() == Index(0))
        {
          _vec_ml.local().clone(lumped_matrix_a_.local(), LAFEM::CloneMode::Layout);
        }
        if(_vec_mr.local().size() == Index(0))
        {
          _vec_mr.local().clone(lumped_matrix_a_.local(), LAFEM::CloneMode::Layout);
        }
        if(inv_lumped_matrix_a.local().size() == Index(0))
        {
          inv_lumped_matrix_a.local().clone(lumped_matrix_a_.local(), LAFEM::CloneMode::Layout);
        }

        inv_lumped_matrix_a.component_invert(lumped_matrix_a_);

      }

      /**
       * \brief Clone operation
       *
       * \param[in] mode
       * The mode for cloning
       *
       * \returns A clone of \c this
       */
      SymmetricLumpedSchurMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return SymmetricSchurMatrix(inv_lumped_matrix_a.clone(mode), matrix_b.clone(mode), matrix_d.clone(mode),
        filter_a.clone(mode));
      }

      /**
       * \brief Returns a left-vector
       *
       * \returns A new left-vector
       */
      VectorTypeL create_vector_l() const
      {
        return matrix_d.create_vector_l();
      }

      /**
       * \brief Returns a right-vector
       *
       * \returns A new right-vector
       */
      VectorTypeR create_vector_r() const
      {
        return matrix_b.create_vector_r();
      }

      /**
       * \brief Gets the total number of columns in this matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of columns
       */
      Index num_cols() const
      {
        return matrix_b.num_cols();
      }

      /**
       * \brief Gets the total number of rows in this matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of columns
       */
      Index num_rows() const
      {
        return matrix_d.num_rows();
      }

      /**
       * \brief Returns the total number of non-zeros in this matrix.
       *
       * \note This always returns the raw (or POD - Plain Old Data) count, as everything else is ambiguous.
       *
       * \returns The total number of nonzeros in this matrix
       */
      Index num_nzes() const
      {
        return inv_lumped_matrix_a.num_nzes() + matrix_b.num_nzes() + matrix_d.num_nzes();
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _vec_ml.bytes() + _vec_mr.bytes();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this~ x \f$
       *
       * \param[out] r
       * The vector that receives the result.
       *
       * \param[in] x
       * The vector to be multiplied by this matrix.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        matrix_b.apply(_vec_mr, x);

        filter_a.filter_def(_vec_mr);
        _vec_ml.component_product(_vec_mr, inv_lumped_matrix_a);
        filter_a.filter_cor(_vec_ml);

        matrix_d.apply(r, _vec_ml);
      }

      /**
       * \brief Calculate \f$ r \leftarrow this~ x \f$
       *
       * \param[out] r
       * The vector that receives the result.
       *
       * \param[in] x
       * The vector to be multiplied by this matrix.
       */
      void preproc_rhs(VectorTypeML& r, const VectorTypeMR& x) const
      {
        _vec_mr.copy(x);
        filter_a.filter_def(_vec_mr);
        _vec_ml.component_product(_vec_mr, inv_lumped_matrix_a);
        filter_a.filter_cor(_vec_ml);

        matrix_d.apply(r, _vec_ml);
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this~ x \f$
       *
       * \param[out] r
       * The vector that receives the result.
       *
       * \param[in] x
       * The vector to be multiplied by this matrix.
       *
       * \param[in] y
       * The summand vector.
       *
       * \param[in] alpha
       * A scalar to scale the product with.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const
      {
        // r <- y + alpha*(D A^(-1) B)*x
        matrix_b.apply(_vec_mr, x);

        filter_a.filter_def(_vec_mr);
        _vec_ml.component_product(_vec_mr, inv_lumped_matrix_a);
        filter_a.filter_cor(_vec_ml);

        matrix_d.apply(r, _vec_ml, y, alpha);
      }

      //LocalMatrix convert_to_1() const
      //{
      //  LocalMatrix locmat = _matrix.clone(LAFEM::CloneMode::Weak);
      //  if((_row_gate != nullptr) && (_col_gate != nullptr))
      //    synch_matrix(locmat, *_row_gate->_comm, _row_gate->_ranks, _row_gate->_mirrors, _col_gate->_mirrors);
      //  return locmat;
      //}

    }; // class SymmetricLumpedSchurMatrix<...>
  } // namespace Global
} // namespace FEAT
