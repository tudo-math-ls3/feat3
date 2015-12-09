#pragma once
#ifndef KERNEL_SOLVER_ILU_PRECOND_HPP
#define KERNEL_SOLVER_ILU_PRECOND_HPP 1

// includes, FEAST
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>

namespace FEAST
{
  namespace Solver
  {
    namespace Intern
    {
      int cuda_ilu_apply(double * y, const double * x, double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo);
      void * cuda_ilu_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd);
      void cuda_ilu_init_numeric(double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo);
      void cuda_ilu_done(void * vinfo);

      int cuda_ilub_apply(double * y, const double * x, double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo);
      void * cuda_ilub_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd, const int blocksize);
      void cuda_ilub_init_numeric(double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo);
      void cuda_ilub_done(void * vinfo);
    }

    template<typename Matrix_, typename Filter_>
    class ILUPrecond;

    /**
     * \brief ILU(0) and ILU(p) preconditioner implementation
     *
     * This class implements a simple ILU(0) and ILU(p)preconditioner.
     *
     * This implementation works for the following matrix types:
     * - LAFEM::SparseMatrixCSR
     * - LAFEM::SparseMatrixELL
     * - LAFEM::SparseMatrixCOO
     *
     * Moreover, this implementation supports only Mem::Main.
     *
     * \author Dirk Ribbrock
     */
    template<typename DT_, typename IT_, typename Filter_>
    class ILUPrecond<LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>, Filter_> :
      public SolverBase<typename LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>::VectorTypeL>
    {
    public:
      typedef DT_ DataType;
      typedef IT_ IndexType;
      typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;

      /// Our memory architecture type
      typedef Mem::Main MemType;

    protected:
      const MatrixType& _A;
      MatrixType _LU;
      const FilterType& _filter;
      const Index _p;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] p level of fillin
       *           if p = 0, the layout of matrix is used for the ILU-decomposition
       */
      explicit ILUPrecond(const MatrixType& matrix, const FilterType& filter, const Index p = 0) :
        _A(matrix),
        _filter(filter),
        _p(p)
      {
        if (_A.columns() != _A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "ILU";
      }

      virtual void init_symbolic() override
      {
        if (_p == 0)
        {
          _LU = LAFEM::SparseMatrixCSR<MemType, DT_, IT_>(_A.layout());
        }
        else
        {
          _symbolic_lu_factorisation((int) _p);
        }
      }

      virtual void done_symbolic() override
      {
      }

      virtual void init_numeric() override
      {
        if (_p == 0)
        {
          _copy_entries(false);
        }
        else
        {
          _copy_entries();
        }
        _create_lu();
      }

      virtual void done_numeric() override
      {
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual Status apply(VectorType& out, const VectorType& in) override
      {
        // copy in-vector to out-vector
        out.copy(in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(_LU.val());
        const IT_ * pcol_ind(_LU.col_ind());
        const IT_ * prow_ptr(_LU.row_ptr());
        const Index n(_LU.rows());

        IT_ col;

        // __forward-insertion__
        // iteration over all rows
        for (Index i(0); i < n; ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          col = prow_ptr[i];

          while (pcol_ind[col] < i)
          {
            pout[i] -= pval[col] * pout[pcol_ind[col]];
            ++col;
          }
        }

        // __backward-insertion__
        // iteration over all rows
        for (Index i(n); i > 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          col = prow_ptr[i+1]-1;

          while (pcol_ind[col] > i)
          {
            pout[i] -= pval[col] * pout[pcol_ind[col]];
            --col;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        this->_filter.filter_cor(out);

        return Status::success;
      } // function apply

    private:
      void _create_lu()
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/2014
         *    section 5.5.3; Algo 5.13.; page 132
         */

        DT_ * plu(_LU.val());
        const IT_ * pcol(_LU.col_ind());
        const IT_ * prow_ptr(_LU.row_ptr());
        const Index n(_LU.rows());

        // integer work array of length n
        //   for saving the position of the diagonal-entries
        IT_ * pw = new IT_[n];

        IT_ row_start;
        IT_ row_end;

        // iteration over all columns
        for (Index i(0); i < n; ++i)
        {
          row_start = prow_ptr[i];
          row_end = prow_ptr[i + 1];

          // iteration over all elements on the left side of the main-diagonal
          // k -> \tilde a_{ik}
          IT_ k = row_start;
          while (pcol[k] < i)
          {
            plu[k] /= plu[pw[pcol[k]]];
            IT_ m(pw[pcol[k]] + 1);
            // m -> \tilde a_{kj}
            for (IT_ j(k+1); j < row_end; ++j)
            {
              // j -> \tilde a_{ij}
              while (m < prow_ptr[pcol[k] + 1])
              {
                if (pcol[m] == pcol[j])
                {
                  plu[j] -= plu[k] * plu[m];
                  ++m;
                  break;
                }
                else if (pcol[m] > pcol[j])
                {
                  break;
                }
                ++m;
              }
            }
            ++k;
          }
          // save the position of the diagonal-entry
          pw[i] = k;
        }

        delete[] pw;
      } // function _create_lu


      void _symbolic_lu_factorisation(int p)
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/2014
         *    section 5.5.6; Algo 5.19.; page 142
         */

        const Index n(_A.rows());
        const IT_ * pacol(_A.col_ind());
        const IT_ * parow(_A.row_ptr());

        // type of list-entries
        typedef std::pair<int, IT_> PAIR_;
        // lists for saving the non-zero entries of L per row
        std::list<PAIR_> * ll = new std::list<PAIR_>[n];

        // vector for saving the iterators to the diag-entries
        typename std::list<PAIR_>::iterator * pldiag = new typename std::list<PAIR_>::iterator[n];
        // auxillary-iterators
        typename std::list<PAIR_>::iterator it, it1, it2;

        IT_ col_begin, col_end, col, col2;
        int l, l2, neues_level;

        // fill list with non-zero entries of A
        // and save iterators to the diag- and last entry of each row
        for (Index row(0); row < n; ++row)
        {
          std::list<PAIR_> & ll_row (ll[row]);
          col_begin = parow[row];
          col_end = parow[row + 1] - 1;

          for (IT_ k(col_begin); k <= col_end; ++k)
          {
            col = pacol[k];

            ll_row.emplace_back((int) n, col);

            if (col == row)
            {
              pldiag[row] = std::prev(ll_row.end());
            }
          }
        }

        // calculate "new" entries of LU
        for (Index row(1); row < n; ++row)
        {
          // iterate from the beginning of the line to the diag-entry
          it = ll[row].begin();
          while (it != pldiag[row])
          {
            col = it->second;
            l = it->first;

            // search non-zero entries in the col-th row
            it1 = it;
            it2 = std::next(pldiag[col]);
            while (it2 != ll[col].end())
            {
              col2 = it2->second;
              l2 = it2->first;

              neues_level = 2* (int) n - l - l2 + 1;

              // if new entries must be created, find the correct position in the list
              if (neues_level <= p)
              {
                while (it1 != ll[row].end() && it1->second < col2)
                {
                  if (it1 != ll[row].end() && it1 == ll[row].end())
                  {
                    ll[row].emplace_back(neues_level, col2);
                    break;
                  }
                  ++it1;
                }
                if (it1 == ll[row].end() || it1->second != col2)
                {
                  it1 = ll[row].emplace(it1, neues_level, col2);
                }
              }
              ++it2;
            }
            ++it;
          }
        }

        // Create LU-matrix
        // calculate number of non-zero-elements
        Index nnz(0);
        for (Index i(0); i < n; ++i)
        {
          nnz += Index(ll[i].size());
        }

        LAFEM::DenseVector<MemType, DT_, IT_> val(nnz);
        LAFEM::DenseVector<MemType, IT_, IT_> col_ind(nnz);
        LAFEM::DenseVector<MemType, IT_, IT_> row_ptr(n+1);
        IT_ * pcol_ind(col_ind.elements());
        IT_ * prow_ptr(row_ptr.elements());

        IT_ k1(0);
        prow_ptr[0] = 0;

        for (Index i(0); i < n; ++i)
        {
          for (it = ll[i].begin(); it != ll[i].end(); ++it)
          {
            pcol_ind[k1] = it->second;
            ++k1;
          }
          prow_ptr[i+1] = k1;
        }

        _LU = LAFEM::SparseMatrixCSR<MemType, DT_, IT_>(n, n, col_ind, val, row_ptr);

        delete[] ll;
        delete[] pldiag;
      } // _symbolic_lu_factorisation

      void _copy_entries(bool check = true)
      {
        if (check == false)
        {
          DT_ * plu(_LU.val());
          const DT_ * pa(_A.val());
          const Index used_elements(_LU.used_elements());

          // initialize LU array to A
          for (Index i(0); i < used_elements; ++i)
          {
            plu[i] = pa[i];
          }
        }
        else
        {
          DT_ * plu(_LU.val());
          const IT_ * plucol(_LU.col_ind());
          const IT_ * plurow(_LU.row_ptr());

          const DT_ * pa(_A.val());
          const IT_ * pacol(_A.col_ind());
          const IT_ * parow(_A.row_ptr());

          const Index n(_LU.rows());
          IT_ k;
          // initialize LU array to A
          for (Index i(0); i < n; ++i)
          {
            k = parow[i];
            for (IT_ j(plurow[i]); j < plurow[i + 1]; ++j)
            {
              plu[j] = DT_(0.0);
              while (k < parow[i + 1] && plucol[j] >= pacol[k])
              {
                if (plucol[j] == pacol[k])
                {
                  plu[j] = pa[k];
                  ++k;
                  break;
                }
                ++k;
              }
            }
          }
        }
      } // function _copy_entries
    }; // class ILUPrecond<SparseMatrixCSR<Mem::Main>>

    /// \todo fix copydoc
    /// \copydoc FEAST::Solver::ILUPrecond< LAFEM::SparseMatrixCSR< Mem::Main, DT_, IT_ >, Filter_ >
    template<typename DT_, typename IT_, typename Filter_>
    class ILUPrecond<LAFEM::SparseMatrixELL<Mem::Main, DT_, IT_>, Filter_> :
      public SolverBase<typename LAFEM::SparseMatrixELL<Mem::Main, DT_, IT_>::VectorTypeL>
    {
    public:
      typedef DT_ DataType;
      typedef IT_ IndexType;
      typedef LAFEM::SparseMatrixELL<Mem::Main, DataType, IndexType> MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;

      /// Our memory architecture type
      typedef Mem::Main MemType;

    protected:
      const MatrixType& _A;
      MatrixType _LU;
      const FilterType& _filter;
      const Index _p;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] p level of fillin
       *           if p = 0, the layout of matrix is used for the ILU-decomposition
       */
      explicit ILUPrecond(const MatrixType& matrix, const FilterType& filter, const Index p = 0) :
        _A(matrix),
        _filter(filter),
        _p(p)
      {
        if (_A.columns() != _A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "ILU";
      }

      virtual void init_symbolic() override
      {
        if (_p == 0)
        {
          _LU = LAFEM::SparseMatrixELL<MemType, DT_, IT_>(_A.layout());
        }
        else
        {
          _symbolic_lu_factorisation((int) _p);
        }
      }

      virtual void done_symbolic() override
      {
      }

      virtual void init_numeric() override
      {
        if (_p == 0)
        {
          _copy_entries(false);
        }
        else
        {
          _copy_entries();
        }
        _create_lu();
      }

      virtual void done_numeric() override
      {
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual Status apply(VectorType& out, const VectorType& in) override
      {
        // copy in-vector to out-vector
        out.copy(in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(_LU.val());
        const IT_ * pcol_ind(_LU.col_ind());
        const IT_ * pcs(_LU.cs());
        const IT_ * prl(_LU.rl());
        const IT_ C((IT_(_LU.C())));
        const IT_ n((IT_(_LU.rows())));

        IT_ col;

        // __forward-insertion__
        // iteration over all rows
        for (IT_ i(0); i < n; ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          col = pcs[i/C] + i%C;

          while (pcol_ind[col] < i)
          {
            pout[i] -= pval[col] * pout[pcol_ind[col]];
            col += C;
          }
        }

        // __backward-insertion__
        // iteration over all rows
        for (IT_ i(n); i > 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          col = pcs[i/C] + i%C + C * (prl[i] - 1);

          while (pcol_ind[col] > i)
          {
            pout[i] -= pval[col] * pout[pcol_ind[col]];
            col -= C;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        this->_filter.filter_cor(out);

        return Status::success;
      } // function apply

    private:
      void _create_lu()
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/2014
         *    section 5.5.3; Algo 5.13.; page 132
         */

        const Index n(_A.rows());
        DT_ * pval(_LU.val());
        const IT_ * pcol_ind(_LU.col_ind());
        const IT_ * pcs(_LU.cs());
        const IT_ * prl(_LU.rl());
        IT_ C(IT_(_LU.C()));

        // integer work array of length n
        //   for saving the position of the diagonal-entries
        IT_ * pw = new IT_[n];

        // iteration over all columns
        for (IT_ i(0); i < IT_(n); ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          // k -> \tilde a_{ik}
          IT_ k(pcs[i/C] + i%C);
          while (pcol_ind[k] < i)
          {
            pval[k] /= pval[pw[pcol_ind[k]]];
            IT_ m(pw[pcol_ind[k]] + C);
            // m -> \tilde a_{kj}
            for (IT_ j(k + C); j < pcs[i/C] + i%C + C*prl[i]; j += C)
            {
              // j -> \tilde a_{ij}
              while (m < pcs[pcol_ind[k]/C] + pcol_ind[k]%C + prl[pcol_ind[k]] * C)
              {
                if (pcol_ind[m] == pcol_ind[j])
                {
                  pval[j] -= pval[k] * pval[m];
                  m += C;
                  break;
                }
                else if (pcol_ind[m] > pcol_ind[j])
                {
                  break;
                }
                m += C;
              }
            }
            k += C;
          }
          // save the position of the diagonal-entry
          pw[i] = k;
        }

        delete[] pw;
      } // function _create_lu

      void _symbolic_lu_factorisation(int p)
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/2014
         *    section 5.5.6; Algo 5.19.; page 142
         */

        const Index n(_A.rows());
        const IT_ * pcol_ind(_A.col_ind());
        const IT_ * pcs(_A.cs());
        const IT_ * prl(_A.rl());
        const Index C(_A.C());

        // type of list-entries
        typedef std::pair<int, IT_> PAIR_;
        // list for saving the non-zero entries of L
        std::list<PAIR_> * ll = new std::list<PAIR_>[n];

        // vector for saving the iterators to the diag-entries
        typename std::list<PAIR_>::iterator * pldiag = new typename std::list<PAIR_>::iterator[n];
        // auxillary-iterators
        typename std::list<PAIR_>::iterator it, it1, it2;

        IT_ col, col2;
        int l, l2, neues_level;

        // fill list with non-zero entries of A
        // and save iterators to the diag- and last entry of each row
        for (Index row(0); row < n; ++row)
        {
          std::list<PAIR_> & ll_row (ll[row]);
          for (Index k(0); k < prl[row]; ++k)
          {
            col = pcol_ind[pcs[row/C] + row%C + C * k];

            ll_row.emplace_back((int) n, col);

            if (col == row)
            {
              pldiag[row] = std::prev(ll_row.end());
            }
          }
        }

        // calculate "new" entries of LU
        for (Index row(1); row < n; ++row)
        {
          // iterate from the beginning of the line to the diag-entry
          it = ll[row].begin();
          while (it != pldiag[row])
          {
            col = it->second;
            l = it->first;

            // search non-zero entries in the col-th row
            it1 = it;
            it2 = std::next(pldiag[col]);
            while (it2 != ll[col].end())
            {
              col2 = it2->second;
              l2 = it2->first;

              neues_level = 2*(int) n - l - l2 + 1;

              // if new entries must be created, find the correct position in the list
              if (neues_level <= p)
              {
                while (it1 != ll[row].end() && it1->second < col2)
                {
                  if (it1 != ll[row].end() && it1 == ll[row].end())
                  {
                    ll[row].emplace_back(neues_level, col2);
                    break;
                  }
                  ++it1;
                }
                if (it1 == ll[row].end() || it1->second != col2)
                {
                  it1 = ll[row].emplace(it1, neues_level, col2);
                }
              }
              ++it2;
            }
            ++it;
          }
        }

        // Create LU-matrix
        // calculate cl-array and fill rl-array
        Index num_of_chunks(Index(ceil(n / float(C))));
        LAFEM::DenseVector<MemType, IT_, IT_> lucl(num_of_chunks, IT_(0));
        LAFEM::DenseVector<MemType, IT_, IT_> lucs(num_of_chunks + 1);
        LAFEM::DenseVector<MemType, IT_, IT_> lurl(n);
        IT_ * plucl(lucl.elements());
        IT_ * plucs(lucs.elements());
        IT_ * plurl(lurl.elements());

        Index nnz(0);

        for (Index i(0); i < n; ++i)
        {
          plurl[i] = IT_(ll[i].size());
          plucl[i/C] = Math::max(plucl[i/C], plurl[i]);
          nnz += Index(ll[i].size());
        }

        // calculuate cs-array
        plucs[0] = IT_(0);
        for (Index i(0); i < num_of_chunks; ++i)
        {
          plucs[i+1] = plucs[i] + IT_(C) * plucl[i];
        }

        Index val_size = Index(plucs[num_of_chunks]);

        LAFEM::DenseVector<MemType, DT_, IT_> luval(val_size);
        LAFEM::DenseVector<MemType, IT_, IT_> lucol_ind(val_size);
        DT_ * pluval    (luval.elements());
        IT_ * plucol_ind(lucol_ind.elements());

        IT_ k1;

        for (Index i(0); i < n; ++i)
        {
          k1 = IT_(0);
          for (it = ll[i].begin(); it != ll[i].end(); ++it, ++k1)
          {
            plucol_ind[plucs[i/C] + i%C + k1 * C] = it->second;
          }
          for (k1 = plurl[i]; k1 < plurl[i]; ++k1)
          {
            plucol_ind[plucs[i/C] + i%C + k1 * C] = IT_(0);
            pluval    [plucs[i/C] + i%C + k1 * C] = DT_(0);
          }
        }

        _LU = LAFEM::SparseMatrixELL<MemType, DT_, IT_>(n, n, nnz, luval, lucol_ind, lucs, lucl, lurl, C);

        delete[] ll;
        delete[] pldiag;
      } // _symbolic_lu_factorisation

      void _copy_entries(bool check = true)
      {
        if (check == false)
        {
          DT_ * pluval(_LU.val());
          const DT_ * paval(_A.val());
          const Index data_length(_LU.val_size());

          // initialize LU array to A
          for (Index i(0); i < data_length; ++i)
          {
            pluval[i] = paval[i];
          }
        }
        else
        {
          DT_ * pluval(_LU.val());
          const IT_ * plucol_ind(_LU.col_ind());
          const IT_ * plucs(_LU.cs());
          const IT_ * plurl(_LU.rl());

          const DT_ * paval(_A.val());
          const IT_ * pacol_ind(_A.col_ind());
          const IT_ * pacs(_A.cs());
          const IT_ * parl(_A.rl());
          const IT_ C((IT_(_A.C())));

          const IT_ n((IT_(_A.rows())));

          IT_ k;
          Index ctr;

          // iteration over all rows
          for (IT_ row(0); row < n; ++row)
          {
            k  = pacs[row/C] + row%C;
            ctr = 0;
            for (IT_ j(0); j < plurl[row]; ++j)
            {
              pluval[plucs[row/C] + row%C + j * C] = DT_(0.0);
              while (ctr < parl[row] && pacol_ind[k] <= plucol_ind[plucs[row/C] + row%C + j * C])
              {
                if (plucol_ind[plucs[row/C] + row%C + j * C] == pacol_ind[k])
                {
                  pluval[plucs[row/C] + row%C + j * C] = paval[k];
                  ++ctr;
                  k = k + C;
                  break;
                }
                ++ctr;
                k = k + C;
              }
            }
          }
        }
      } // function _copy_entries
    }; // class ILUPrecond<SparseMatrixELL<Mem::Main>>

    /**
     * \brief ILU(0) preconditioner implementation
     *
     * This class implements a simple ILU(0) preconditioner,
     * e.g. zero fill-in and no pivoting.
     *
     * This implementation works for the following matrix types and combinations thereof:
     * - LAFEM::SparseMatrixCSR
     *
     * Moreover, this implementation supports only CUDA and uint containers.
     *
     * \author Dirk Ribbrock
     */
    template<typename Filter_>
    class ILUPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned int>, Filter_> :
      public SolverBase<LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned int>::VectorTypeL>
    {
    public:
      typedef LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned int> MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef typename MatrixType::DataType DataType;

    protected:
      const MatrixType& _matrix;
      MatrixType _lu_matrix;
      const FilterType& _filter;
      void * cuda_info;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] filter
       * The filter to be used for the correction vector.
       *
       * \note The parameter p is not active in CUDA, i.e. only ILU(0) is supported.
       */
      explicit ILUPrecond(const MatrixType& matrix, const FilterType& filter, const Index = 0) :
        _matrix(matrix),
        _filter(filter)
      {
        if (_matrix.columns() != _matrix.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "ILU";
      }

      virtual void init_symbolic() override
      {
        _lu_matrix.clone(_matrix, LAFEM::CloneMode::Layout);

        cuda_info = Intern::cuda_ilu_init_symbolic((int)_lu_matrix.rows(), (int)_lu_matrix.used_elements(), _lu_matrix.val(),
          (int*)_lu_matrix.row_ptr(), (int*)_lu_matrix.col_ind());
      }

      virtual void done_symbolic() override
      {
        Intern::cuda_ilu_done(cuda_info);
      }

      virtual void init_numeric() override
      {
        _lu_matrix.copy(_matrix);

        Intern::cuda_ilu_init_numeric(_lu_matrix.val(), (int*)_lu_matrix.row_ptr(), (int*)_lu_matrix.col_ind(), cuda_info);
      }

      virtual void done_numeric() override
      {
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        ASSERT(_matrix.rows() == vec_cor.size(), "Error: matrix / vector size missmatch!");
        ASSERT(_matrix.rows() == vec_def.size(), "Error: matrix / vector size missmatch!");

        int status = Intern::cuda_ilu_apply(vec_cor.elements(), vec_def.elements(), _lu_matrix.val(), (int*)_lu_matrix.row_ptr(), (int*)_lu_matrix.col_ind(), cuda_info);

        this->_filter.filter_cor(vec_cor);
        return (status == 0) ? Status::success :  Status::aborted;
      }
    }; // class ILUPrecond<SparseMatrixCSR<Mem::CUDA>>

    /**
     * \brief ILU(0) preconditioner implementation
     *
     * This class implements a simple ILU(0) preconditioner,
     * e.g. zero fill-in and no pivoting.
     *
     * This implementation works for the following matrix types and combinations thereof:
     * - LAFEM::SparseMatrixBCSR
     *
     * Moreover, this implementation supports only CUDA and uint containers.
     *
     * \author Dirk Ribbrock
     */
    template<typename Filter_, int blocksize_>
    class ILUPrecond<LAFEM::SparseMatrixBCSR<Mem::CUDA, double, unsigned int, blocksize_, blocksize_>, Filter_> :
      public SolverBase<typename LAFEM::SparseMatrixBCSR<Mem::CUDA, double, unsigned int, blocksize_, blocksize_>::VectorTypeL>
    {
    public:
      typedef LAFEM::SparseMatrixBCSR<Mem::CUDA, double, unsigned int, blocksize_, blocksize_> MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef typename MatrixType::DataType DataType;

    protected:
      const MatrixType& _matrix;
      MatrixType _lu_matrix;
      const FilterType& _filter;
      void * cuda_info;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] filter
       * The filter to be used for the correction vector.
       *
       * \note The parameter p is not active in CUDA, i.e. only ILU(0) is supported.
       */
      explicit ILUPrecond(const MatrixType& matrix, const FilterType& filter, const Index = 0) :
        _matrix(matrix),
        _filter(filter)
      {
        if (_matrix.columns() != _matrix.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "ILU";
      }

      virtual void init_symbolic() override
      {
        _lu_matrix.clone(_matrix, LAFEM::CloneMode::Layout);

        cuda_info = Intern::cuda_ilub_init_symbolic((int)_lu_matrix.rows(), (int)_lu_matrix.used_elements(), _lu_matrix.template val<LAFEM::Perspective::pod>(),
          (int*)_lu_matrix.row_ptr(), (int*)_lu_matrix.col_ind(), blocksize_);
      }

      virtual void done_symbolic() override
      {
        Intern::cuda_ilub_done(cuda_info);
      }

      virtual void init_numeric() override
      {
        _lu_matrix.copy(_matrix);

        Intern::cuda_ilub_init_numeric(_lu_matrix.template val<LAFEM::Perspective::pod>(), (int*)_lu_matrix.row_ptr(), (int*)_lu_matrix.col_ind(), cuda_info);
      }

      virtual void done_numeric() override
      {
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        ASSERT(_matrix.rows() == vec_cor.size(), "Error: matrix / vector size missmatch!");
        ASSERT(_matrix.rows() == vec_def.size(), "Error: matrix / vector size missmatch!");

        int status = Intern::cuda_ilub_apply(vec_cor.template elements<LAFEM::Perspective::pod>(), vec_def.template elements<LAFEM::Perspective::pod>(),
          _lu_matrix.template val<LAFEM::Perspective::pod>(), (int*)_lu_matrix.row_ptr(), (int*)_lu_matrix.col_ind(), cuda_info);

        this->_filter.filter_cor(vec_cor);
        return (status == 0) ? Status::success :  Status::aborted;
      }
    }; // class ILUPrecond<SparseMatrixBCSR<Mem::CUDA>>

    /**
     * \brief Creates a new ILUPrecond solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] p level of fillin
     *           if p = 0, the layout of matrix is used for the ILU-decomposition
     *
     * \returns
     * A shared pointer to a new ILUPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<ILUPrecond<Matrix_, Filter_>> new_ilu_precond(
      const Matrix_& matrix, const Filter_& filter, const Index p = 0)
    {
      return std::make_shared<ILUPrecond<Matrix_, Filter_>>(matrix, filter, p);
    }

  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_ILU_PRECOND_HPP
