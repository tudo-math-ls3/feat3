// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    namespace Intern
    {
      int cuda_ilu_apply(double * y, const double * x, double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo);
      void * cuda_ilu_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd);
      void cuda_ilu_init_numeric(double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo);
      void cuda_ilu_done_symbolic(void * vinfo);

      int cuda_ilub_apply(double * y, const double * x, double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo);
      void * cuda_ilub_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd, const int blocksize);
      void cuda_ilub_init_numeric(double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo);
      void cuda_ilub_done_symbolic(void * vinfo);

      /**
       * \brief ILU(k) symbolic core
       *
       * This class is responsible for the computation and management
       * of the symbolic structure of an ILU factorization.
       *
       * \tparam IT_
       * The index type to be used.
       *
       * \author Peter Zajac
       */
      template<typename IT_>
      class ILUCoreSymbolic
      {
      protected:
        /// The size of the system matrix
        IT_ _n;
        /// CSR structure of L
        std::vector<IT_> _row_ptr_l, _col_idx_l;
        /// CSR structure of U
        std::vector<IT_> _row_ptr_u, _col_idx_u;
        //std::vector<IT_> _lvl_l, _lvl_u;

      public:
        /// Clears all symbolic data arrays
        void clear()
        {
          _n = IT_(0);
          _row_ptr_l.clear();
          _col_idx_l.clear();
          _row_ptr_u.clear();
          _col_idx_u.clear();
          //_lvl_l.clear();
          //_lvl_u.clear();
        }

        /// Returns the number of non-zeros in L.
        IT_ get_nnze_l() const
        {
          return _row_ptr_l.empty() ? IT_(0) : _row_ptr_l.back();
        }

        /// Returns the number of non-zeros in U.
        IT_ get_nnze_u() const
        {
          return _row_ptr_u.empty() ? IT_(0) : _row_ptr_u.back();
        }

        /// Returns the number of non-zeros in L, U and D.
        IT_ get_nnze() const
        {
          return _n + get_nnze_l() + get_nnze_u();
        }

        /// Returns the size of the symbolic factorization in bytes.
        std::size_t bytes_symbolic() const
        {
          return sizeof(IT_) * (_row_ptr_l.size() + _row_ptr_u.size() + _col_idx_l.size() + _col_idx_u.size());
        }

        /**
         * \brief Initializes the ILU(0) structure from a CSR input matrix
         *
         * \param[in] n
         * The size of the matrix.
         *
         * \param[in] row_ptr_a
         * The row-pointer array of the CSR input matrix.
         *
         * \param[in] col_idx_a
         * The column-index array of the CSR input matrix.
         */
        void set_struct_csr(const IT_ n, const IT_* row_ptr_a, const IT_* col_idx_a)
        {
          // First of all, let's guess the number of non-zeros for L and U:
          // In most of our cases, the sparsity pattern of A is symmetric,
          // so we have a guess of (nz(A)-n)/2 for both L and U:
          const IT_ nzul = (Math::max(row_ptr_a[n], n) - n) / IT_(2);

          // store size and clear containers
          _n = n;
          _row_ptr_l.clear();
          _col_idx_l.clear();
          _row_ptr_u.clear();
          _col_idx_u.clear();

          // reserve our non-zero guess
          _row_ptr_l.reserve(n+1);
          _row_ptr_u.reserve(n+1);
          _col_idx_l.reserve(nzul);
          _col_idx_u.reserve(nzul);

          // initialize row-pointers
          _row_ptr_l.push_back(0);
          _row_ptr_u.push_back(0);

          // loop over all matrix rows
          for(IT_ i(0); i < n; ++i)
          {
            IT_ j = row_ptr_a[i];
            const IT_ jend = row_ptr_a[i+1];
            for(; (j < jend) && (col_idx_a[j] < i); ++j)
            {
              _col_idx_l.push_back(col_idx_a[j]);
            }
            if(col_idx_a[j] != i)
              throw InvalidMatrixStructureException();
            for(++j; j < jend; ++j)
            {
              _col_idx_u.push_back(col_idx_a[j]);
            }
            _row_ptr_l.push_back(IT_(_col_idx_l.size()));
            _row_ptr_u.push_back(IT_(_col_idx_u.size()));
          }

          //_lvl_l.resize(_col_idx_l.size(), IT_(0));
          //_lvl_u.resize(_col_idx_u.size(), IT_(0));
        }

        /**
         * \brief Initializes the ILU(0) structure from a CSR input matrix
         *
         * \param[in] matrix
         * The CSR input matrix.
         */
        template<typename DT_>
        void set_struct(const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix)
        {
          this->set_struct_csr(IT_(matrix.rows()), matrix.row_ptr(), matrix.col_ind());
        }

        /**
         * \brief Initializes the ILU(0) structure from a BCSR input matrix
         *
         * \param[in] matrix
         * The BCSR input matrix.
         */
        template<typename DT_, int m_, int n_>
        void set_struct(const LAFEM::SparseMatrixBCSR<DT_, IT_, m_, n_>& matrix)
        {
          this->set_struct_csr(IT_(matrix.rows()), matrix.row_ptr(), matrix.col_ind());
        }

        /**
         * \brief Performs symbolic ILU(p) factorization
         *
         * This function performs the symbolic factorization of the ILU(p) sparsity patterns.
         * This function treats the existing pattern (usually initialized by a call to one of
         * the #set_struct() functions) as the level-0 pattern and computes the level-p pattern
         * from it.
         *
         * \param[in] p
         * The maximum allowed level of fill-in.
         *
         * \note
         * This function does nothing if p < 1.
         */
        void factorize_symbolic(const int p)
        {
          if(p < 1)
            return;

          // our new level-p structures for L and U
          std::vector<IT_> new_ptr_l, new_idx_l, new_ptr_u, new_idx_u;

          // our temporary level arrays
          std::vector<int> new_lvl_l, new_lvl_u;

          // initialize new row pointers
          new_ptr_l.push_back(IT_(0));
          new_ptr_u.push_back(IT_(0));

          // loop over all rows
          for(IT_ i(0); i < _n; ++i)
          {
            // first off, add all level 0 entries to L and U
            for(IT_ j(_row_ptr_l[i]); j < _row_ptr_l[i+1]; ++j)
            {
              new_idx_l.push_back(_col_idx_l[j]);
              new_lvl_l.push_back(0);
            }
            for(IT_ j(_row_ptr_u[i]); j < _row_ptr_u[i+1]; ++j)
            {
              new_idx_u.push_back(_col_idx_u[j]);
              new_lvl_u.push_back(0);
            }

            // loop over row L_i
            for(IT_ j(new_ptr_l[i]); j < IT_(new_idx_l.size()); ++j)
            {
              // get column index and level of L_ij
              const IT_ cj = new_idx_l[j];
              const int lj = new_lvl_l[j];
              IT_ olj = j;
              IT_ ouj = new_ptr_u[i];

              // loop over row U_j
              for(IT_ k(new_ptr_u[cj]); k < new_ptr_u[cj+1]; ++k)
              {
                // get column index and level of U_jk
                const IT_ ck = new_idx_u[k];
                const int lk = new_lvl_u[k];

                // compute level of entry L/U_ik
                const int ll = lj + lk + 1;
                if(ll > p)
                  continue;

                // insert into L or U
                if(ck < i)
                  olj = _insert(new_idx_l, new_lvl_l, olj, ck, ll);
                else if(ck > i)
                  ouj = _insert(new_idx_u, new_lvl_u, ouj, ck, ll);
              }
            }

            // update row-pointers
            new_ptr_l.push_back(IT_(new_idx_l.size()));
            new_ptr_u.push_back(IT_(new_idx_u.size()));
          }

          // replace old arrays
          _row_ptr_l = std::move(new_ptr_l);
          _col_idx_l = std::move(new_idx_l);
          _row_ptr_u = std::move(new_ptr_u);
          _col_idx_u = std::move(new_idx_u);

          //_lvl_l = new_lvl_l;
          //_lvl_u = new_lvl_u;
        }

      protected:
        /**
         * \brief Linear insertion function
         *
         * This function inserts a new column-index into the column-index vector
         * of the to-be-computed sparsity pattern of L or U.
         *
         * \param[inout] idx
         * The new column-index vector of L or U.
         *
         * \param[inout] lvl
         * The level vector of L or U.
         *
         * \param[in] i
         * The starting position for the search in idx/lvl.
         *
         * \param[in] j
         * The column index of the entry to be inserted
         *
         * \param[in] l
         * The level of the entry to be inserted
         *
         * \returns
         * The starting position for the next column entry.
         */
        static IT_ _insert(std::vector<IT_>& idx, std::vector<int>& lvl, IT_ i, const IT_ j, const int l)
        {
          // get the current length of the column-index vector
          IT_ n = IT_(idx.size());

          // loop over the row and try to find our entry
          for(; i < n; ++i)
          {
            if(idx[i] == j)
            {
              // the entry that we want to add already exists in the sparsity pattern,
              // so check whether we need to adapt its level
              if(l < lvl[i])
                lvl[i] = l;

              // return next entry starting position
              return ++i;
            }
            else if(j < idx[i])
            {
              // we have reached a column index beyond the entry we want to add,
              // so stop searching
              break;
            }
          }

          // the entry we want to add does not yet exist in the sparsity pattern,
          // so we add it at position i by linear insertion:
          idx.push_back(IT_(0));
          lvl.push_back(0);

          // shift all entries backward by one position
          for(IT_ k(n); k > i; --k)
          {
            idx[k] = idx[k-1];
            lvl[k] = lvl[k-1];
          }

          // insert at desired position
          idx[i] = j;
          lvl[i] = l;

          // return next entry starting position
          return ++i;
        }
      }; // class ILUCoreSymbolic

      /**
       * \brief ILU Core for Scalar Matrices
       *
       * This class is responsible for the factorization and management of an ILU(p) factorization.
       *
       * \tparam DT_
       * The data-type to be used.
       *
       * \tparam IT_
       * The index-type to be used.
       *
       * \author Peter Zajac
       */
      template<typename DT_, typename IT_>
      class ILUCoreScalar :
        public ILUCoreSymbolic<IT_>
      {
      protected:
        typedef ILUCoreSymbolic<IT_> BaseClass;

        /// The data arrays of L, U and D.
        std::vector<DT_> _data_l, _data_u, _data_d;

      public:
        /// Clears all data arrays.
        void clear()
        {
          BaseClass::clear();
          _data_l.clear();
          _data_u.clear();
          _data_d.clear();
        }

        /// Allocates the data arrays for the numeric factorization
        void alloc_data()
        {
          // resize data arrays
          _data_d.resize(this->_n);
          _data_l.resize(this->_col_idx_l.size());
          _data_u.resize(this->_col_idx_u.size());
        }

        /// Returns the size of the numeric factorization in bytes
        std::size_t bytes_numeric() const
        {
          return sizeof(DT_) * (_data_l.size() + _data_u.size() + _data_d.size());
        }

        /// Returns the total size of the factorization in bytes
        std::size_t bytes() const
        {
          return this->bytes_numeric() + this->bytes_symbolic();
        }

        /**
         * \brief Copies the data arrays from an CSR input matrix
         *
         * \param[in] row_ptr_a
         * The row-pointer array of the CSR input matrix.
         *
         * \param[in] col_idx_a
         * The column-index array of the CSR input matrix.
         *
         * \param[in] data_a
         * The data arrays of the CSR input matrix.
         */
        void copy_data_csr(const IT_* row_ptr_a, const IT_* col_idx_a, const DT_* data_a)
        {
          // loop over all rows
          for(IT_ i(0); i < this->_n; ++i)
          {
            // get row pointer of a
            IT_ ra = row_ptr_a[i];
            IT_ xa = row_ptr_a[i+1];

            // fetch row i of L
            for(IT_ j(this->_row_ptr_l[i]); j < this->_row_ptr_l[i+1]; ++j)
            {
              if(this->_col_idx_l[j] == col_idx_a[ra])
                _data_l[j] = data_a[ra++];
              else //if(this->_col_idx_l[j] < col_idx_a[ra])
                _data_l[j] = DT_(0);
            }

            // fetch diagonal of a
            _data_d[i] = data_a[ra++];

            // fetch row i of U
            for(IT_ j(this->_row_ptr_u[i]); j < this->_row_ptr_u[i+1]; ++j)
            {
              if((ra < xa) && (this->_col_idx_u[j] == col_idx_a[ra]))
                _data_u[j] = data_a[ra++];
              else// if(this->_col_idx_u[j] < col_idx_a[ra])
                _data_u[j] = DT_(0);
            }
          }
        }

        /**
         * \brief Copies the data arrays from a CSR input matrix
         *
         * \param[in] matrix
         * The input matrix.
         */
        void copy_data(const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix)
        {
          this->copy_data_csr(matrix.row_ptr(), matrix.col_ind(), matrix.val());
        }

        /**
         * \brief Performs the (I+L)*(D+U) numeric factorization
         *
         * This function performs the "classic" (I+L)*(D+U) numeric factorization of the
         * input matrix that has been set by a prior call to one of the #copy_data functions.
         */
        void factorize_numeric_il_du()
        {
          // get data arrays
          const IT_* rptr_l = this->_row_ptr_l.data();
          const IT_* rptr_u = this->_row_ptr_u.data();
          const IT_* cidx_l = (this->_col_idx_l.empty() ? nullptr : this->_col_idx_l.data());
          const IT_* cidx_u = (this->_col_idx_u.empty() ? nullptr : this->_col_idx_u.data());
          DT_* data_l = (this->_data_l.empty() ? nullptr : this->_data_l.data());
          DT_* data_u = (this->_data_u.empty() ? nullptr : this->_data_u.data());
          DT_* data_d = this->_data_d.data();

          // loop over all rows
          for(IT_ i(0); i < this->_n; ++i)
          {
            // get row-end pointers of L and U
            const IT_ ql = rptr_l[i+1];
            const IT_ qu = rptr_u[i+1];

            // loop over row i of L
            for(IT_ j(rptr_l[i]); j < rptr_l[i+1]; ++j)
            {
              // get column index of L_ij
              const IT_ cj = cidx_l[j];

              // get L/U pointers
              IT_ pl = j;
              IT_ pu = rptr_u[i];

              // update L_ij <- L_ij / D_jj
              data_l[j] *= data_d[cj];

              // loop over row j of U and process row i of L
              IT_ k(rptr_u[cj]);
              for(; k < rptr_u[cj+1]; ++k)
              {
                const IT_ ck = cidx_u[k];
                if(ck >= i)
                  break;
                for(; (pl < ql) && (cidx_l[pl] <= ck); ++pl)
                {
                  if(cidx_l[pl] == ck)
                    data_l[pl] -= data_l[j] * data_u[k];
                }
              }

              // process main diagonal entry
              if((k < rptr_u[cj+1]) && (cidx_u[k] == i))
              {
                data_d[i] -= data_l[j] * data_u[k];
                ++k;
              }
              // loop over row j of U and process row i of U
              for(; k < rptr_u[cj+1]; ++k)
              {
                const IT_ ck = cidx_u[k];
                for(; (pu < qu) && (cidx_u[pu] <= ck); ++pu)
                {
                  if(cidx_u[pu] == ck)
                    data_u[pu] -= data_l[j] * data_u[k];
                }
              }
            }

            // invert main diagonal entry
            data_d[i] = DT_(1) / data_d[i];
          }
        }

        /**
         * \brief Solves (I+L)*x = b
         *
         * \param[out] x
         * The solution vector
         *
         * \param[in] b
         * The right-hand-side vector
         *
         * \note
         * \p x and \p b are allowed to refer to the same array.
         */
        void solve_il(DT_* x, const DT_* b) const
        {

          const IT_* rptr = this->_row_ptr_l.data();
          const IT_* cidx = (this->_col_idx_l.empty() ? nullptr : this->_col_idx_l.data());
          const DT_* data_l = (this->_data_l.empty() ? nullptr : this->_data_l.data());

          for(IT_ i(0); i < this->_n; ++i)
          {
            DT_ r(b[i]);
            for(IT_ j(rptr[i]); j < rptr[i+1]; ++j)
            {
              r -= data_l[j] * x[cidx[j]];
            }
            x[i] = r;
          }
        }

        /**
         * \brief Solves (D+U)*x = b
         *
         * \param[out] x
         * The solution vector
         *
         * \param[in] b
         * The right-hand-side vector
         *
         * \note
         * \p x and \p b are allowed to refer to the same array.
         */
        void solve_du(DT_* x, const DT_* b) const
        {
          const IT_* rptr = this->_row_ptr_u.data();
          const IT_* cidx = (this->_col_idx_u.empty() ? nullptr : this->_col_idx_u.data());
          const DT_* data_u = (this->_data_u.empty() ? nullptr : this->_data_u.data());
          const DT_* data_d = this->_data_d.data();

          for(IT_ i(this->_n); i > IT_(0); )
          {
            --i;
            DT_ r(b[i]);
            for(IT_ j(rptr[i]); j < rptr[i+1]; ++j)
            {
              r -= data_u[j] * x[cidx[j]];
            }
            x[i] = data_d[i] * r;
          }
        }

        /**
         * \brief Solves (I+L^T)*x = b
         *
         * \param[out] x
         * The solution vector
         *
         * \param[in] b
         * The right-hand-side vector
         *
         * \note
         * \p x and \p b are allowed to refer to the same array.
         */
        void solve_ilt(DT_* x, const DT_* b) const
        {
          const IT_* rptr = this->_row_ptr_l.data();
          const IT_* cidx = (this->_col_idx_l.empty() ? nullptr : this->_col_idx_l.data());
          const DT_* data_l = (this->_data_l.empty() ? nullptr : this->_data_l.data());

          if(x != b)
          {
            for(IT_ i(0); i < this->_n; ++i)
              x[i] = b[i];
          }

          for(IT_ i(this->_n); i > IT_(0); )
          {
            --i;
            for(IT_ j(rptr[i]); j < rptr[i+1]; ++j)
            {
              // x_j -= L_ij * x[i]
              x[cidx[j]] -= data_l[j] * x[i];
            }
          }
        }

        /**
         * \brief Solves (D+U^T)*x = b
         *
         * \param[out] x
         * The solution vector
         *
         * \param[in] b
         * The right-hand-side vector
         *
         * \note
         * \p x and \p b are allowed to refer to the same array.
         */
        void solve_dut(DT_* x, const DT_* b) const
        {
          const IT_* rptr = this->_row_ptr_u.data();
          const IT_* cidx = (this->_col_idx_u.empty() ? nullptr : this->_col_idx_u.data());
          const DT_* data_u = (this->_data_u.empty() ? nullptr : this->_data_u.data());
          const DT_* data_d = this->_data_d.data();

          if(x != b)
          {
            for(IT_ i(0); i < this->_n; ++i)
              x[i] = b[i];
          }

          for(IT_ i(0); i < this->_n; ++i)
          {
            x[i] *= data_d[i];
            for(IT_ j(rptr[i]); j < rptr[i+1]; ++j)
            {
              // x_j -= U_ij * x_i
              x[cidx[j]] -= data_u[j] * x[i];
            }
          }
        }
      }; // class ILUCoreScalar

      /**
       * \brief ILU Core for Blocked Matrices
       *
       * This class is responsible for the factorization and management of an ILU(p) factorization.
       *
       * \tparam DT_
       * The data-type to be used.
       *
       * \tparam IT_
       * The index-type to be used.
       *
       * \tparam dim_
       * The row/column block dimension
       *
       * \author Peter Zajac
       */
      template<typename DT_, typename IT_, int dim_>
      class ILUCoreBlocked :
        public ILUCoreSymbolic<IT_>
      {
      protected:
        typedef ILUCoreSymbolic<IT_> BaseClass;

        static_assert(dim_ > 0, "invalid block dimension");

        typedef Tiny::Matrix<DT_, dim_, dim_> MatBlock;
        typedef Tiny::Vector<DT_, dim_> VecBlock;

        /// The data arrays of L, U and D.
        std::vector<MatBlock> _data_l, _data_u, _data_d;

      public:
        /// Clears all data arrays.
        void clear()
        {
          BaseClass::clear();
          _data_l.clear();
          _data_u.clear();
          _data_d.clear();
        }

        /// Allocates the data arrays for the numeric factorization
        void alloc_data()
        {
          // resize data arrays
          _data_d.resize(this->_n);
          _data_l.resize(this->_col_idx_l.size());
          _data_u.resize(this->_col_idx_u.size());
        }

        /// Returns the size of the numeric factorization in bytes
        std::size_t bytes_numeric() const
        {
          return sizeof(DT_) * std::size_t(Math::sqr(dim_)) * (_data_l.size() + _data_u.size() + _data_d.size());
        }

        /// Returns the total size of the factorization in bytes
        std::size_t bytes() const
        {
          return this->bytes_numeric() + this->bytes_symbolic();
        }

        /**
         * \brief Copies the data arrays from an CSR input matrix
         *
         * \param[in] row_ptr_a
         * The row-pointer array of the CSR input matrix.
         *
         * \param[in] col_idx_a
         * The column-index array of the CSR input matrix.
         *
         * \param[in] data_a
         * The data arrays of the CSR input matrix.
         */
        void copy_data_bcsr(const IT_* row_ptr_a, const IT_* col_idx_a, const MatBlock* data_a)
        {
          // loop over all rows
          for(IT_ i(0); i < this->_n; ++i)
          {
            // get row pointer of a
            IT_ ra = row_ptr_a[i];
            IT_ xa = row_ptr_a[i+1];

            // fetch row i of L
            for(IT_ j(this->_row_ptr_l[i]); j < this->_row_ptr_l[i+1]; ++j)
            {
              if(this->_col_idx_l[j] == col_idx_a[ra])
                _data_l[j] = data_a[ra++];
              else //if(this->_col_idx_l[j] < col_idx_a[ra])
                _data_l[j] = DT_(0);
            }

            // fetch diagonal of a
            _data_d[i] = data_a[ra++];

            // fetch row i of U
            for(IT_ j(this->_row_ptr_u[i]); j < this->_row_ptr_u[i+1]; ++j)
            {
              if((ra < xa) && (this->_col_idx_u[j] == col_idx_a[ra]))
                _data_u[j] = data_a[ra++];
              else// if(this->_col_idx_u[j] < col_idx_a[ra])
                _data_u[j] = DT_(0);
            }
          }
        }

        /**
         * \brief Copies the data arrays from a CSR input matrix
         *
         * \param[in] matrix
         * The input matrix.
         */
        void copy_data(const LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, dim_>& matrix)
        {
          this->copy_data_bcsr(matrix.row_ptr(), matrix.col_ind(), matrix.val());
        }

        /**
         * \brief Performs the (I+L)*(D+U) numeric factorization
         *
         * This function performs the "classic" (I+L)*(D+U) numeric factorization of the
         * input matrix that has been set by a prior call to one of the #copy_data functions.
         */
        void factorize_numeric_il_du()
        {
          // get data arrays
          const IT_* rptr_l = this->_row_ptr_l.data();
          const IT_* rptr_u = this->_row_ptr_u.data();
          const IT_* cidx_l = (this->_col_idx_l.empty() ? nullptr : this->_col_idx_l.data());
          const IT_* cidx_u = (this->_col_idx_u.empty() ? nullptr : this->_col_idx_u.data());
          MatBlock* data_l = (this->_data_l.empty() ? nullptr : this->_data_l.data());
          MatBlock* data_u = (this->_data_u.empty() ? nullptr : this->_data_u.data());
          MatBlock* data_d = this->_data_d.data();

          // loop over all rows
          for(IT_ i(0); i < this->_n; ++i)
          {
            // get row-end pointers of L and U
            const IT_ ql = rptr_l[i+1];
            const IT_ qu = rptr_u[i+1];

            // loop over row i of L
            for(IT_ j(rptr_l[i]); j < rptr_l[i+1]; ++j)
            {
              // get column index of L_ij
              const IT_ cj = cidx_l[j];

              // get L/U pointers
              IT_ pl = j;
              IT_ pu = rptr_u[i];

              // update L_ij <- D_jj^{-1} * L_ij
              {
                //data_l[j] *= data_d[cj];
                const MatBlock l_ij(data_l[j]);
                data_l[j].set_mat_mat_mult(data_d[cj], l_ij);
              }

              // loop over row j of U and process row i of L
              IT_ k(rptr_u[cj]);
              for(; k < rptr_u[cj+1]; ++k)
              {
                const IT_ ck = cidx_u[k];
                if(ck >= i)
                  break;
                for(; (pl < ql) && (cidx_l[pl] <= ck); ++pl)
                {
                  if(cidx_l[pl] == ck)
                    //data_l[pl] -= data_l[j] * data_u[k];
                    data_l[pl].add_mat_mat_mult(data_l[j], data_u[k], -DT_(1));
                }
              }

              // process main diagonal entry
              if((k < rptr_u[cj+1]) && (cidx_u[k] == i))
              {
                //data_d[i] -= data_l[j] * data_u[k];
                data_d[i].add_mat_mat_mult(data_l[j], data_u[k], -DT_(1));
                ++k;
              }

              // loop over row j of U and process row i of U
              for(; k < rptr_u[cj+1]; ++k)
              {
                const IT_ ck = cidx_u[k];
                for(; (pu < qu) && (cidx_u[pu] <= ck); ++pu)
                {
                  if(cidx_u[pu] == ck)
                    //data_u[pu] -= data_l[j] * data_u[k];
                    data_u[pu].add_mat_mat_mult(data_l[j], data_u[k], -DT_(1));
                }
              }
            }

            // invert main diagonal entry
            {
              //data_d[i] = DT_(1) / data_d[i];
              const MatBlock d_ii(data_d[i]);
              data_d[i].set_inverse(d_ii);
            }
          }
        }

        /**
         * \brief Solves (I+L)*x = b
         *
         * \param[out] x
         * The solution vector
         *
         * \param[in] b
         * The right-hand-side vector
         *
         * \note
         * \p x and \p b are allowed to refer to the same array.
         */
        void solve_il(VecBlock* x, const VecBlock* b) const
        {
          const IT_* rptr = this->_row_ptr_l.data();
          const IT_* cidx = (this->_col_idx_l.empty() ? nullptr : this->_col_idx_l.data());
          const MatBlock* data_l = (this->_data_l.empty() ? nullptr : this->_data_l.data());

          for(IT_ i(0); i < this->_n; ++i)
          {
            VecBlock r(b[i]);
            for(IT_ j(rptr[i]); j < rptr[i+1]; ++j)
            {
              //r -= data_l[j] * x[cidx[j]];
              r.add_mat_vec_mult(data_l[j], x[cidx[j]], -DT_(1));
            }
            x[i] = r;
          }
        }

        /**
         * \brief Solves (D+U)*x = b
         *
         * \param[out] x
         * The solution vector
         *
         * \param[in] b
         * The right-hand-side vector
         *
         * \note
         * \p x and \p b are allowed to refer to the same array.
         */
        void solve_du(VecBlock* x, const VecBlock* b) const
        {
          const IT_* rptr = this->_row_ptr_u.data();
          const IT_* cidx = (this->_col_idx_u.empty() ? nullptr : this->_col_idx_u.data());
          const MatBlock* data_u = (this->_data_u.empty() ? nullptr : this->_data_u.data());
          const MatBlock* data_d = this->_data_d.data();

          for(IT_ i(this->_n); i > IT_(0); )
          {
            --i;
            VecBlock r(b[i]);
            for(IT_ j(rptr[i]); j < rptr[i+1]; ++j)
            {
              //r -= data_u[j] * x[cidx[j]];
              r.add_mat_vec_mult(data_u[j], x[cidx[j]], -DT_(1));
            }
            //x[i] = data_d[i] * r;
            x[i].set_mat_vec_mult(data_d[i], r);
          }
        }
      }; // class ILUCoreBlocked
    } // namespace Intern
    /// \endcond

    /**
     * \brief Inheritances inside ilu_precond.hpp
     *
     * ILUPrecondBase <- ILUPrecondWithBackend.
     * Both for internal use only.
     *
     * SolverBase <- ILUPrecond.
     * ILUPrecond uses the internal classes.
     *
     * \author Dirk Ribbrock
     * \author Lisa-Marie Walter
     */

    /// ILU preconditioner base class for internal use
    template<typename Matrix_>
    class ILUPrecondBase
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::VectorTypeL VectorType;

      virtual void set_fill_in_param(int p) = 0;

      virtual void init_symbolic() = 0;

      virtual void done_symbolic() = 0;

      virtual void init_numeric() = 0;

      virtual String name() const = 0;

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) = 0;

      virtual ~ILUPrecondBase() {}
    }; // ILUPrecondBase class

    /**
     * \brief ILU(0) and ILU(p) preconditioner internal implementation
     *
     * This class implements a simple ILU(0) and ILU(p) preconditioner.
     *
     * \author Dirk Ribbrock
     * \author Peter Zajac
     */
    template<PreferredBackend backend_, typename Matrix_, typename Filter_>
    class ILUPrecondWithBackend;

    /**
     * \brief ILU(p) specialization for SparseMatrixCSR
     *
     * This specialization works for the LAFEM::SparseMatrixCSR matrix type
     *
     * \author Peter Zajac
     */
    template<typename DT_, typename IT_, typename Filter_>
    class ILUPrecondWithBackend<PreferredBackend::generic, LAFEM::SparseMatrixCSR<DT_, IT_>, Filter_> :
      public ILUPrecondBase<typename LAFEM::SparseMatrixCSR<DT_, IT_>>
    {
    public:
      typedef LAFEM::SparseMatrixCSR<DT_, IT_> MatrixType;
      typedef DT_ DataType;
      typedef IT_ IndexType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;

    protected:
      const MatrixType& _matrix;
      const FilterType& _filter;
      Intern::ILUCoreScalar<DataType, IndexType> _ilu;
      int _p;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] p
       * Maximum level of fillin.
       */
      explicit ILUPrecondWithBackend(const MatrixType& matrix, const FilterType& filter, const int p = 0) :
        _matrix(matrix),
        _filter(filter),
        _p(p)
      {
      }

      explicit ILUPrecondWithBackend(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter) :
        _matrix(matrix),
        _filter(filter),
        _p(-1)
      {
        auto fill_in_param_p = section->query("fill_in_param");
        if(fill_in_param_p.second && !fill_in_param_p.first.parse(_p))
          throw ParseError(section_name + ".fill_in_param", fill_in_param_p.first, "a non-negative integer");
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~ILUPrecondWithBackend()
      {
      }

      /**
       * \brief Sets the fill-in parameter
       *
       * \param[in] p
       * The new fill-in parameter
       */
      virtual void set_fill_in_param(int p) override
      {
        XASSERT(p > 0);
        _p = p;
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "ILU";
      }

      virtual void init_symbolic() override
      {
        if (_matrix.columns() != _matrix.rows())
        {
          XABORTM("Matrix is not square!");
        }

        // set matrix structure
        _ilu.set_struct(_matrix);

        // perform symbolic factorization
        _ilu.factorize_symbolic(_p);

        // allocate data arrays
        _ilu.alloc_data();
      }

      virtual void done_symbolic() override
      {
        _ilu.clear();
      }

      virtual void init_numeric() override
      {
        _ilu.copy_data(_matrix);
        _ilu.factorize_numeric_il_du();
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual Status apply(VectorType& out, const VectorType& in) override
      {
        TimeStamp ts_start;

        // get vector data arrays
        DataType* x = out.elements();
        const DataType* b = in.elements();

        // solve: (I+L)*y = b
        _ilu.solve_il(x, b);

        // solve: (D+U)*x = y
        _ilu.solve_du(x, x);

        // apply filter
        this->_filter.filter_cor(out);

        TimeStamp ts_stop;
        Statistics::add_time_precon(ts_stop.elapsed(ts_start));
        Statistics::add_flops(_ilu.get_nnze() * 2 + out.size());

        return Status::success;
      }
    }; // class ILUPrecondWithBackend<generic, SparseMatrixCSR,...>

    /**
     * \brief ILU(p) specialization for SparseMatrixBCSR<...>
     *
     * \author Peter Zajac
     */
    template<typename DT_, typename IT_, int dim_, typename Filter_>
    class ILUPrecondWithBackend<PreferredBackend::generic, LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, dim_>, Filter_> :
      public ILUPrecondBase<typename LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, dim_>>
    {
    public:
      typedef LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, dim_> MatrixType;
      typedef DT_ DataType;
      typedef IT_ IndexType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;

    protected:
      const MatrixType& _matrix;
      const FilterType& _filter;
      Intern::ILUCoreBlocked<DataType, IndexType, dim_> _ilu;
      int _p;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] p
       * Maximum level of fillin.
       */
      explicit ILUPrecondWithBackend(const MatrixType& matrix, const FilterType& filter, const int p = 0) :
        _matrix(matrix),
        _filter(filter),
        _p(p)
      {
      }

      explicit ILUPrecondWithBackend(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter) :
        _matrix(matrix),
        _filter(filter),
        _p(-1)
      {
        auto fill_in_param_p = section->query("fill_in_param");
        if(fill_in_param_p.second && !fill_in_param_p.first.parse(_p))
          throw ParseError(section_name + ".fill_in_param", fill_in_param_p.first, "a non-negative integer");
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~ILUPrecondWithBackend()
      {
      }

      /**
       * \brief Sets the fill-in parameter
       *
       * \param[in] p
       * The new fill-in parameter
       */
      virtual void set_fill_in_param(int p) override
      {
        XASSERT(p > 0);
        _p = p;
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "ILU";
      }

      virtual void init_symbolic() override
      {
        if (_matrix.columns() != _matrix.rows())
        {
          XABORTM("Matrix is not square!");
        }

        // set matrix structure
        _ilu.set_struct(_matrix);

        // perform symbolic factorization
        _ilu.factorize_symbolic(_p);

        // allocate data arrays
        _ilu.alloc_data();
      }

      virtual void done_symbolic() override
      {
        _ilu.clear();
      }

      virtual void init_numeric() override
      {
        _ilu.copy_data(_matrix);
        _ilu.factorize_numeric_il_du();
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual Status apply(VectorType& out, const VectorType& in) override
      {
        TimeStamp ts_start;

        // get vector data arrays
        auto* x = out.elements();
        const auto* b = in.elements();

        // solve: (I+L)*y = b
        _ilu.solve_il(x, b);

        // solve: (D+U)*x = y
        _ilu.solve_du(x, x);

        // apply filter
        this->_filter.filter_cor(out);

        TimeStamp ts_stop;
        Statistics::add_time_precon(ts_stop.elapsed(ts_start));
        Statistics::add_flops(_ilu.get_nnze() * 2 + out.size());

        return Status::success;
      }
    }; // class ILUPrecondWithBackend<generic, SparseMatrixBCSR<...>,...>

#ifdef FEAT_HAVE_CUDA
    /**
     * \brief ILU(0) preconditioner internal implementation
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
    class ILUPrecondWithBackend<PreferredBackend::cuda, LAFEM::SparseMatrixCSR<double, unsigned int>, Filter_> :
      public ILUPrecondBase<LAFEM::SparseMatrixCSR<double, unsigned int>>
    {
    public:
      typedef LAFEM::SparseMatrixCSR<double, unsigned int> MatrixType;
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
      explicit ILUPrecondWithBackend(const MatrixType& matrix, const FilterType& filter, const int = 0) :
        _matrix(matrix),
        _filter(filter)
      {
      }

      explicit ILUPrecondWithBackend(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter) :
        _matrix(matrix),
        _filter(filter)
      {
        // Check for _p
        auto fill_in_param_p = section->query("fill_in_param");
        if(fill_in_param_p.second)
        {
          XASSERTM(std::stoi(fill_in_param_p.first) == 0, "For PreferredBackend::cuda, the fill in parameter has to be == 0!");
        }
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~ILUPrecondWithBackend()
      {
      }

      /**
       * \brief Sets the fill-in parameter
       *
       * \param[in] p
       * The new fill-in parameter
       */
      virtual void set_fill_in_param(int p) override
      {
        XASSERTM(p == 0, "For PreferredBackend::cuda, the fill in parameter has to be == 0!");
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "ILU";
      }

      virtual void init_symbolic() override
      {
        if (_matrix.columns() != _matrix.rows())
        {
          XABORTM("Matrix is not square!");
        }

        _lu_matrix.clone(_matrix, LAFEM::CloneMode::Layout);

        cuda_info = Intern::cuda_ilu_init_symbolic(
          (int)_lu_matrix.rows(),
          (int)_lu_matrix.used_elements(),
          _lu_matrix.val(),
          (int*)_lu_matrix.row_ptr(),
          (int*)_lu_matrix.col_ind());
      }

      virtual void done_symbolic() override
      {
        Intern::cuda_ilu_done_symbolic(cuda_info);
      }

      virtual void init_numeric() override
      {
        _lu_matrix.copy(_matrix);

        Intern::cuda_ilu_init_numeric(
          _lu_matrix.val(),
          (int*)_lu_matrix.row_ptr(),
          (int*)_lu_matrix.col_ind(),
          cuda_info);
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        XASSERTM(_matrix.rows() == vec_cor.size(), "matrix / vector size mismatch!");
        XASSERTM(_matrix.rows() == vec_def.size(), "matrix / vector size mismatch!");

        TimeStamp ts_start;

        int status = Intern::cuda_ilu_apply(
          vec_cor.elements(),
          vec_def.elements(),
          _lu_matrix.val(),
          (int*)_lu_matrix.row_ptr(),
          (int*)_lu_matrix.col_ind(),
          cuda_info);

        this->_filter.filter_cor(vec_cor);

        TimeStamp ts_stop;
        Statistics::add_time_precon(ts_stop.elapsed(ts_start));
        Statistics::add_flops(_matrix.used_elements() * 2 + vec_cor.size());

        return (status == 0) ? Status::success :  Status::aborted;
      }
    }; // class ILUPrecondWithBackend<cuda, SparseMatrixCSR<...>,...>

    /*
    template<typename Filter_>
    class ILUPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned int>, Filter_> :
      public SolverBase<LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned int>::VectorTypeL>
    {
      public:
      typedef LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned int> MatrixType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef Filter_ FilterType;

      explicit ILUPrecond(const MatrixType&, const FilterType&, const int = 0)
      {
      }

      Status apply(VectorType &, const VectorType &) override
      {
          XABORTM("not implemented yet!");
      }

      String name() const override
      {
          XABORTM("not implemented yet!");
      }
    };

    template<typename Filter_>
    class ILUPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned long>, Filter_> :
      public SolverBase<LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned long>::VectorTypeL>
    {
      public:
      typedef LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned long> MatrixType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef Filter_ FilterType;

      explicit ILUPrecond(const MatrixType&, const FilterType&, const int = 0)
      {
      }

      Status apply(VectorType &, const VectorType &) override
      {
          XABORTM("not implemented yet!");
      }

      String name() const override
      {
          XABORTM("not implemented yet!");
      }
    };

    template<typename Filter_>
    class ILUPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned long>, Filter_> :
      public SolverBase<LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned long>::VectorTypeL>
    {
      public:
      typedef LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned long> MatrixType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef Filter_ FilterType;

      explicit ILUPrecond(const MatrixType&, const FilterType&, const int = 0)
      {
      }

      Status apply(VectorType &, const VectorType &) override
      {
          XABORTM("not implemented yet!");
      }

      String name() const override
      {
          XABORTM("not implemented yet!");
      }
    };*/

    /**
     * \brief ILU(0) preconditioner internal implementation
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
    class ILUPrecondWithBackend<PreferredBackend::cuda, LAFEM::SparseMatrixBCSR<double, unsigned int, blocksize_, blocksize_>, Filter_> :
      public ILUPrecondBase<typename LAFEM::SparseMatrixBCSR<double, unsigned int, blocksize_, blocksize_>>
    {
    public:
      typedef LAFEM::SparseMatrixBCSR<double, unsigned int, blocksize_, blocksize_> MatrixType;
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
      explicit ILUPrecondWithBackend(const MatrixType& matrix, const FilterType& filter, const int = 0) :
        _matrix(matrix),
        _filter(filter)
      {
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] matrix
       * The system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       */
      explicit ILUPrecondWithBackend(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter) :
        _matrix(matrix),
        _filter(filter)
        {
          // Check for _p
          auto fill_in_param_p = section->query("fill_in_param");
          if(fill_in_param_p.second)
          {
            XASSERTM(std::stoi(fill_in_param_p.first) == 0, "For PreferredBackend::cuda, the fill in parameter has to be == 0!");
          }
        }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~ILUPrecondWithBackend()
      {
      }

      /**
       * \brief Sets the fill-in parameter
       *
       * \param[in] p
       * The new fill-in parameter
       */
      virtual void set_fill_in_param(int p) override
      {
        XASSERTM(p == 0, "For PreferredBackend::cuda, the fill in parameter has to be == 0!");
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "ILU";
      }

      virtual void init_symbolic() override
      {
        if (_matrix.columns() != _matrix.rows())
        {
          XABORTM("Matrix is not square!");
        }

        _lu_matrix.clone(_matrix, LAFEM::CloneMode::Layout);

        cuda_info = Intern::cuda_ilub_init_symbolic(
          (int)_lu_matrix.rows(),
          (int)_lu_matrix.used_elements(),
          _lu_matrix.template val<LAFEM::Perspective::pod>(),
          (int*)_lu_matrix.row_ptr(),
          (int*)_lu_matrix.col_ind(),
          blocksize_);
      }

      virtual void done_symbolic() override
      {
        Intern::cuda_ilub_done_symbolic(cuda_info);
      }

      virtual void init_numeric() override
      {
        _lu_matrix.copy(_matrix);

        Intern::cuda_ilub_init_numeric(_lu_matrix.template val<LAFEM::Perspective::pod>(), (int*)_lu_matrix.row_ptr(), (int*)_lu_matrix.col_ind(), cuda_info);
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        XASSERTM(_matrix.rows() == vec_cor.size(), "matrix / vector size mismatch!");
        XASSERTM(_matrix.rows() == vec_def.size(), "matrix / vector size mismatch!");

        TimeStamp ts_start;

        int status = Intern::cuda_ilub_apply(
          vec_cor.template elements<LAFEM::Perspective::pod>(),
          vec_def.template elements<LAFEM::Perspective::pod>(),
          _lu_matrix.template val<LAFEM::Perspective::pod>(),
          (int*)_lu_matrix.row_ptr(),
          (int*)_lu_matrix.col_ind(),
          cuda_info);

        this->_filter.filter_cor(vec_cor);

        TimeStamp ts_stop;
        Statistics::add_time_precon(ts_stop.elapsed(ts_start));
        Statistics::add_flops(_matrix.used_elements() * 2 + vec_cor.size());

        return (status == 0) ? Status::success :  Status::aborted;
      }
    }; // class ILUPrecondWithBackend<cuda, SparseMatrixBCSR<...>,...>
#endif //FEAT_HAVE_CUDA

    /// Dummy class for not implemented specializations
    template<PreferredBackend backend_, typename Matrix_, typename Filter_>
    class ILUPrecondWithBackend :
      public ILUPrecondBase<Matrix_>
    {
      public:

      explicit ILUPrecondWithBackend(const Matrix_&, const Filter_&, const int = 0)
      {
      }

      explicit ILUPrecondWithBackend(const String& , PropertyMap*, const Matrix_& , const Filter_& )
      {
      }

      virtual Status apply(typename Matrix_::VectorTypeL &, const typename Matrix_::VectorTypeL &) override
      {
          XABORTM("not implemented yet!");
      }

      virtual String name() const override
      {
          XABORTM("not implemented yet!");
      }

      virtual void init_symbolic() override
      {
        XABORTM("not implemented yet!");
      }

      virtual void done_symbolic() override
      {
        XABORTM("not implemented yet!");
      }

      virtual void set_fill_in_param(int /*p*/) override
      {
        XABORTM("not implemented yet!");
      }

      virtual void init_numeric() override
      {
        XABORTM("not implemented yet!");
      }
    };

    /**
     * \brief ILU(0) and ILU(p) preconditioner implementation
     *
     * This class implements a simple ILU(0) and ILU(p) preconditioner.
     *
     * This implementation works for the following matrix types and combinations thereof:
     * - LAFEM::SparseMatrixCSR
     * - LAFEM::SparseMatrixBCSR
     *
     * \author Dirk Ribbrock
     */
    template<typename Matrix_, typename Filter_>
    class ILUPrecond : public SolverBase<typename Matrix_::VectorTypeL>
    {
    private:
      std::shared_ptr<ILUPrecondBase<Matrix_>> _impl;

    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef typename MatrixType::DataType DataType;

      /// Our base class
      typedef SolverBase<VectorType> BaseClass;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] backend
       * The backend to be preferred. This implementation works with generic and cuda.
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] filter
       * The filter to be used for the correction vector.
       *
       * \param[in] p
       * Maximum level of fillin.
       *
       */
      ILUPrecond(PreferredBackend backend, const MatrixType& matrix, const FilterType& filter, const int p = 0)
      {
        switch (backend)
        {
          case PreferredBackend::cuda:
            _impl = std::make_shared<ILUPrecondWithBackend<PreferredBackend::cuda, Matrix_, Filter_>>(matrix, filter);
            break;
          case PreferredBackend::mkl:
          case PreferredBackend::generic:
          default:
            _impl = std::make_shared<ILUPrecondWithBackend<PreferredBackend::generic, Matrix_, Filter_>>(matrix, filter, p);
        }
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] backend
       * The backend to be preferred. This implementation works with generic and cuda.
       *
       * \param[in] matrix
       * The system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       */
      ILUPrecond(const String& section_name, PropertyMap* section, PreferredBackend backend, const MatrixType& matrix, const FilterType& filter) :
        BaseClass(section_name, section)
      {
        switch (backend)
        {
          case PreferredBackend::cuda:
            _impl = std::make_shared<ILUPrecondWithBackend<PreferredBackend::cuda, Matrix_, Filter_>>(section_name, section, matrix, filter);
            break;
          case PreferredBackend::mkl:
          case PreferredBackend::generic:
          default:
            _impl = std::make_shared<ILUPrecondWithBackend<PreferredBackend::generic, Matrix_, Filter_>>(section_name, section, matrix, filter);
        }
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~ILUPrecond()
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return _impl->name();
      }

      /**
       * \brief Sets the fill-in parameter
       *
       * \param[in] p
       * The new fill-in parameter
       */
      virtual void set_fill_in_param(int p)
      {
        _impl->set_fill_in_param(p);
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        return _impl->apply(vec_cor, vec_def);
      }

      virtual void init_symbolic() override
      {
        _impl->init_symbolic();
      }

      virtual void done_symbolic() override
      {
        _impl->done_symbolic();
      }

      virtual void init_numeric() override
      {
        _impl->init_numeric();
      }
    }; // class ILUPrecond

    /**
     * \brief Creates a new ILUPrecond solver object
     *
     * \param[in] backend
     * The backend to be preferred. This implementation works with generic and cuda.
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] p
     * Maximum level of fill-in.
     *
     * \returns
     * A shared pointer to a new ILUPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<ILUPrecond<Matrix_, Filter_>> new_ilu_precond(PreferredBackend backend,
      const Matrix_& matrix, const Filter_& filter, const int p = 0)
    {
      return std::make_shared<ILUPrecond<Matrix_, Filter_>>(backend, matrix, filter, p);
    }

    /**
     * \brief Creates a new ILUPrecond solver object based on a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] backend
     * The backend to be preferred. This implementation works with generic and cuda.
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new ILUPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<ILUPrecond<Matrix_, Filter_>> new_ilu_precond(
      const String& section_name, PropertyMap* section, PreferredBackend backend,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<ILUPrecond<Matrix_, Filter_>>(section_name, section, backend, matrix, filter);
    }
  } // namespace Solver
} // namespace FEAT
