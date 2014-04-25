#pragma once
#ifndef KERNEL_LAFEM_TRANSPOSITION_HPP
#define KERNEL_LAFEM_TRANSPOSITION_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct Transposition
    {
    };

    /**
     * \brief Summation calculations.
     *
     * This class transposes a given matrix ex-situ.
     *
     * \author Christoph Lohmann
     */
    template <>
    struct Transposition<Algo::Generic>
    {
      /**
       * \brief Calculate \f$A^T\f$
       *
       * \param[in] a The input matrix.
       * \returns The transpose of a.
       */
      template <typename DT_, typename IT_>
      static SparseMatrixCOO<Mem::Main, DT_, IT_> value(const SparseMatrixCOO<Mem::Main, DT_, IT_> & a)
      {
        typedef Mem::Main Mem_;

        const Index arows(a.rows());
        const Index acolumns(a.columns());
        const Index used_elements(a.used_elements());

        const IT_ * parow_ind(a.row_indices());
        const IT_ * pacol_ind(a.column_indices());
        const DT_ * paval(a.val());

        DenseVector<Mem_, IT_, IT_> tcol_ind(used_elements);
        DenseVector<Mem_, DT_, IT_> tval(used_elements);
        DenseVector<Mem_, IT_, IT_> trow_ind(used_elements);

        IT_ * ptcol_ind(tcol_ind.elements());
        DT_ * ptval(tval.elements());
        IT_ * ptrow_ind(trow_ind.elements());

        DenseVector<Mem_, IT_, IT_> trow_ptr(acolumns + 1, Index(0));
        IT_ * ptrow_ptr(trow_ptr.elements());

        for (Index i(0); i < used_elements; ++i)
        {
          ++ptrow_ptr[pacol_ind[i] + 1];
        }

        ptrow_ptr[0] = 0;

        for (Index i(1); i < acolumns - 1; ++i)
        {
          ptrow_ptr[i + 1] += ptrow_ptr[i];
        }

        for (Index i(0); i < used_elements; ++i)
        {
          const Index l(pacol_ind[i]);
          const Index j(ptrow_ptr[l]);
          ptval[j] = paval[i];
          ptcol_ind[j] = parow_ind[i];
          ptrow_ind[j] = pacol_ind[i];
          ++ptrow_ptr[l];
        }

        SparseMatrixCOO<Mem_, DT_, IT_> t(acolumns, arows, trow_ind, tcol_ind, tval);
        return t;
      }

      /**
       * \brief Calculate \f$A^T\f$
       *
       * \param[in] a The input matrix.
       * \returns The transpose of a.
       */
      template <typename DT_, typename IT_>
      static SparseMatrixCSR<Mem::Main, DT_, IT_> value(const SparseMatrixCSR<Mem::Main, DT_, IT_> & a)
      {
        typedef Mem::Main Mem_;

        const Index arows(a.rows());
        const Index acolumns(a.columns());
        const Index used_elements(a.used_elements());

        const IT_ * pacol_ind(a.col_ind());
        const IT_ * parow_ptr(a.row_ptr());
        const DT_ * paval(a.val());

        DenseVector<Mem_, IT_, IT_> tcol_ind(used_elements);
        DenseVector<Mem_, DT_, IT_> tval(used_elements);
        DenseVector<Mem_, IT_, IT_> trow_ptr(acolumns + 1, Index(0));

        IT_ * ptcol_ind(tcol_ind.elements());
        DT_ * ptval(tval.elements());
        IT_ * ptrow_ptr(trow_ptr.elements());

        ptrow_ptr[0] = 0;

        for (Index i(0); i < used_elements; ++i)
        {
          ++ptrow_ptr[pacol_ind[i] + 1];
        }

        for (Index i(1); i < acolumns - 1; ++i)
        {
          ptrow_ptr[i + 1] += ptrow_ptr[i];
        }

        for (Index i(0); i < arows; ++i)
        {
          for (Index k(parow_ptr[i]); k < parow_ptr[i+1]; ++k)
          {
            const Index l(pacol_ind[k]);
            const Index j(ptrow_ptr[l]);
            ptval[j] = paval[k];
            ptcol_ind[j] = i;
            ++ptrow_ptr[l];
          }
        }

        for (Index i(acolumns); i > 0; --i)
        {
          ptrow_ptr[i] = ptrow_ptr[i - 1];
        }
        ptrow_ptr[0] = 0;

        SparseMatrixCSR<Mem_, DT_, IT_> t(acolumns, arows, tcol_ind, tval, trow_ptr);
        return t;
      }

      /**
       * \brief Calculate \f$A^T\f$
       *
       * \param[in] a The input matrix.
       * \returns The transpose of a.
       */
      template <typename DT_, typename IT_>
      static SparseMatrixELL<Mem::Main, DT_, IT_> value(const SparseMatrixELL<Mem::Main, DT_, IT_> & a)
      {
        typedef Mem::Main Mem_;

        const Index arows(a.rows());
        const Index acolumns(a.columns());
        const Index used_elements(a.used_elements());
        const Index astride(a.stride());

        const DT_ * pax(a.Ax());
        const IT_ * paj(a.Aj());
        const IT_ * parl(a.Arl());

        const Index alignment(32);
        const Index tstride(alignment * ((acolumns + alignment - 1)/ alignment));

        DenseVector<Mem_, IT_, IT_> trl(acolumns, Index(0));
        IT_ * ptrl(trl.elements());

        for (Index i(0); i < arows; ++i)
        {
          for (Index j(i); j < i + parl[i] * astride; j += astride)
          {
            ++ptrl[paj[j]];
          }
        }

        Index num_cols_per_row(0);
        for (Index i(0); i < acolumns; ++i)
        {
          if (num_cols_per_row < ptrl[i])
          {
            num_cols_per_row = ptrl[i];
          }
          ptrl[i] = 0;
        }

        DenseVector<Mem_, IT_, IT_> tj(tstride * num_cols_per_row);
        DenseVector<Mem_, DT_, IT_> tx(tstride * num_cols_per_row);

        IT_ * ptj(tj.elements());
        DT_ * ptx(tx.elements());

        for (Index i(0); i < arows; ++i)
        {
          for (Index j(i); j < i + parl[i] * astride; j += astride)
          {
            const Index k(paj[j]);
            ptj[k + ptrl[k] * tstride] = i;
            ptx[k + ptrl[k] * tstride] = pax[j];
            ++ptrl[k];
          }
        }

        SparseMatrixELL<Mem_, DT_, IT_> t(acolumns, arows, tstride, num_cols_per_row, used_elements, tx, tj, trl);
        return t;
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SUM_HPP
