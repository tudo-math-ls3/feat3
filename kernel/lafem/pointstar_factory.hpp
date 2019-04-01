// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_POINTSTAR_FACTORY_HPP
#define KERNEL_LAFEM_POINTSTAR_FACTORY_HPP 1

#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/lafem/pointstar_structure.hpp>

#include <stack>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Pointstar factory base class
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_ = Index>
    class PointstarFactoryBase
    {
    protected:
      /// number of inner nodes per dimension
      Index _m;
      /// number of dimensions
      Index _d;

      PointstarFactoryBase(Index m, Index d) :
        _m(m), _d(d)
      {
        XASSERTM(m >= Index(3), "You need at least 3 inner nodes per dimension");
        XASSERTM(d >= Index(1), "You need at least 1 dimension");
      }

      virtual ~PointstarFactoryBase() {}

    public:
      /**
       * \brief Computes the pointstar matrix in CSR format.
       *
       * \returns
       * The m^d x m^d pointstar matrix
       */
      virtual SparseMatrixCSR<Mem::Main, DataType_, IndexType_> matrix_csr() const = 0;

      /**
       * \brief Computes the smallest eigenvalue of the pointstar matrix.
       *
       * \returns
       * The smallest eigenvalue of the matrix.
       */
      virtual DataType_ lambda_min() const = 0;

      /**
       * \brief Computes the largest eigenvalue of the pointstar matrix.
       *
       * \returns
       * The largest eigenvalue of the matrix.
       */
      virtual DataType_ lambda_max() const = 0;

      /**
       * \brief Computes the spectral condition number of the pointstar matrix.
       *
       * \returns
       * The spectral condition number of the matrix.
       */
      virtual DataType_ spectral_cond() const
      {
        return lambda_max() / lambda_min();
      }

      /**
       * \brief Computes the eigenvector of the smallest eigenvalue
       *
       * This function generates the eigenvector with respect to the smallest eigenvalue
       * of the corresponding pointstar matrix.
       *
       * \returns
       * The m^d x m^d eigenvector
       */
      virtual DenseVector<Mem::Main, DataType_, IndexType_> eigenvector_min() const
      {
        // compute vector length
        Index neq(1);
        for(Index i(0); i < _d; ++i)
          neq *= _m;

        // create vector
        DenseVector<Mem::Main, DataType_, IndexType_> vector(neq, DataType_(1));
        DataType_* v = vector.elements();

        // compute x scaling factor
        const DataType_ xsf = DataType_(Math::pi<DataType_>()) / DataType_(_m+1);

        for(Index i(0); i < neq; ++i)
        {
          for(Index quo(1); quo < neq; quo *= _m)
          {
            v[i] *= Math::sin(xsf * DataType_(((i / quo) % _m) + Index(1)));
          }
        }

        // return vector
        return std::move(vector);
      }

      /**
       * \brief Computes a Q2-bubble vector.
       *
       * This function can be used to generate a useful solution vector for the pointstar matrix.
       * The corresponding right hand side can be computed by multiplying this vector by the matrix.
       *
       * \returns
       * The m^d x m^d Q2-bubble vector
       */
      DenseVector<Mem::Main, DataType_, IndexType_> vector_q2_bubble() const
      {
        // compute vector length
        Index neq(1);
        for(Index i(0); i < _d; ++i)
          neq *= _m;

        // create vector
        DenseVector<Mem::Main, DataType_, IndexType_> vector(neq, DataType_(1));
        DataType_* v = vector.elements();

        // compute x scaling factor
        const DataType_ xsf = DataType_(1) / DataType_(_m+1);

        for(Index i(0); i < neq; ++i)
        {
          for(Index quo(1); quo < neq; quo *= _m)
          {
            // compute coordinate x in (0,1)
            DataType_ x(xsf * DataType_(((i / quo) % _m) + Index(1)));
            v[i] *= DataType_(4) * x * (DataType_(1) - x);
          }
        }

        // return vector
        return std::move(vector);
      }
    }; // PointstarFactoryBase

    /**
     * \brief Finite-Differences pointstar matrix factory.
     *
     * This class generates a CSR poinstar matrices in Finite-Differences style.
     * Moreover, this class can compute the largest and smallest eigenvalue of the
     * corresponding matrix as well as the eigenvector with respect to the smallest
     * eigenvalue.
     *
     * This class works for arbitrary dimensions; for instance, calling this class'
     * constructor with the parameter \p d = 2 will lead to the famous 2D 5-pointstar.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_ = Index>
    class PointstarFactoryFD :
      public PointstarFactoryBase<DataType_, IndexType_>
    {
    public:
      /**
       * \brief Creates the pointstar factory.
       *
       * \param[in] m
       * The number of inner nodes per dimension. Must be >= 3.
       *
       * \param[in] d
       * The dimension of the grid. Must be >= 1.
       */
      PointstarFactoryFD(Index m, Index d = Index(2)) :
        PointstarFactoryBase<DataType_, IndexType_>(m, d)
      {
      }

      virtual ~PointstarFactoryFD() {}

      /**
       * \brief Generates a FD-style pointstar CSR matrix
       *
       * \returns
       * The m^d x m^d FD-stye pointstar matrix.
       */
      virtual SparseMatrixCSR<Mem::Main, DataType_, IndexType_> matrix_csr() const override
      {
        // declare some variables
        Index i,j,k,l,n,neq,nnze,bpl,rpb,rps,b;
        IndexType_ *epr,*dor,*ior,*col_idx;
        DataType_ *data;
        const Index m(this->_m), d(this->_d);

        // calculate number of equations and non-zero entries
        neq = m;
        nnze = 3*m - 2;
        for(i = 1; i < d; ++i)
        {
          nnze *= m;
          nnze += 2*(m-1) * neq;
          neq *= m;
        }

        // create five vectors
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_row_ptr(neq+IndexType_(1));
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_col_idx(nnze);
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_dor(neq, IndexType_(0));
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_epr(neq, IndexType_(1));
        DenseVector<Mem::Main, DataType_, IndexType_> vec_data(nnze);

        // get data
        data = vec_data.elements();

        // get epr
        epr = vec_epr.elements();

        // get dor
        dor = vec_dor.elements();

        // get ior
        ior = vec_row_ptr.elements();

        // get column indices
        col_idx = vec_col_idx.elements();

        /* In this step we will recursively calculate the number of non-zero
         * elements in a row (= epr), and the number of non-zero elements left of
         * the main diagonal element in a row (= dor).
         * The indices of the first element of a row (= ior) and the indices of
         * the main diagonal elements (= dor) will later be recalculated from
         * epr and dor.
         */

        /* blocks per level */
        bpl = 1;

        /* rows per block */
        rpb = neq;

        /* go through all levels */
        for(l = 0; l < d; ++l)
        {
          /* rows per sub-block */
          rps = rpb / m;

          /* go through all blocks */
          for(b = 0; b < bpl; ++b)
          {
            /* index of b-th block */
            k = rpb * b;

            /* add the right diagonal to sub-block 1 */
            for(j = 0; j < rps; ++j)
              ++epr[k+j];

            /* add left and right diagonals to sub-blocks 2..m-1 */
            for(i = 1; i < m - 1; ++i)
            {
              /* index of i-th sub-block */
              n = k + rps*i;
              for(j = 0; j < rps; ++j)
              {
                epr[n+j] += 2;
                ++dor[n+j];
              }
            }

            /* add left diagonal to sub-block m */
            n = k + rps*(m - 1);
            for(j = 0; j < rps; ++j)
            {
              ++epr[n+j];
              ++dor[n+j];
            }
          }

          /* increase blocks per level */
          bpl *= m;

          /* decrease rows per block */
          rpb = rps;
        }

        /* initialize row and diagonal pointer arrays */
        ior[0] = 0;
        for(i = 0; i < neq; ++i)
        {
          ior[i+1] = ior[i] + epr[i];
          dor[i] += ior[i];
        }

        /* initialize data array */
        const DataType_ v1 = DataType_(2*d);
        const DataType_ v2 = DataType_(-1);
        for(i = 0; i < neq; ++i)
        {
          /* left of main diagonal */
          for(j = ior[i]; j < dor[i]; ++j)
            data[j] = v2;

          /* main diagonal */
          data[dor[i]] = v1;

          /* right of main diagonal */
          for(j = dor[i]+1; j < ior[i+1]; ++j)
            data[j] = v2;
        }

        /* We will now overwrite epr with the number of non-zero elements
         * right of the main diagonal element of a row.
         * At the same time, we can already set the indices of the main
         * diagonal elements.
         */
        for(i = 0; i < neq; ++i)
        {
          epr[i] = ior[i+1] - dor[i] - 1;
          col_idx[dor[i]] = IndexType_(i);
        }

        /* Now we need to calculate the indices for the off-main-diagonal
         * non-zero elements.
         * Once again, we will do this recursively ^_^
         *
         * We will go through every level of the matrix, and in each level,
         * we will set the indices of the diagonal-matrix blocks which to
         * not belong to the main diagonal.
         */

        /* blocks per level */
        bpl = 1;

        /* rows per block */
        rpb = neq;

        /* go through all levels */
        for(l = 0; l < d; ++l)
        {
          /* rows per sub-block */
          rps = rpb / m;

          /* go through all blocks */
          for(b = 0; b < bpl; ++b)
          {
            /* index of b-th block */
            k = rpb * b;

            /* go through the first m-1 sub-blocks */
            for(i = 0; i < rpb - rps; ++i)
            {
              /* and this is where the magic happens */
              col_idx[dor[k+i]+epr[k+i]] = IndexType_(k + i + rps);
              col_idx[dor[neq-k-i-1]-epr[k+i]] = IndexType_(neq - k - i - rps - 1);
            }

            for(i = 0; i < rpb - rps; ++i)
              --epr[k+i];
          }

          /* increase blocks per level */
          bpl *= m;

          /* decrease rows per block */
          rpb = rps;
        }

        // return the matrix
        return SparseMatrixCSR<Mem::Main, DataType_, IndexType_>(neq, neq, vec_col_idx, vec_data, vec_row_ptr);
      }

      /**
       * \brief Computes the smallest eigenvalue of the FD-style matrix.
       *
       * The smallest eigenvalue of the FD-style matrix wrt the sine-bubble eigenvector is given as
       *
       *   \f[ \lambda_{\min} = 2 \cdot d \cdot (1 - \cos\big(\pi/(m+1)\big) ) \f]
       *
       * \returns
       * The smallest eigenvalue of the matrix.
       */
      virtual DataType_ lambda_min() const override
      {
        // 4 * d * sin(pi/(2*m+2))^2  =  2 * d * (1 - cos(pi/(m+1)))
        return DataType_(2) * DataType_(this->_d) *
          (DataType_(1) - Math::cos(Math::pi<DataType_>() / DataType_(this->_m+1)));
      }

      /**
       * \brief Computes the largest eigenvalue of the FD-style matrix.
       *
       * The largest eigenvalue of the FD-style matrix is given as
       *
       * \f[ \lambda_{\max} = 2 \cdot d \cdot (1 + \cos\big(\pi/(m+1)\big) ) \f]
       *
       * \returns
       * The largest eigenvalue of the matrix.
       */
      virtual DataType_ lambda_max() const override
      {
        // 2 * d * (1 + cos(pi/(m+1)))
        return DataType_(2) * DataType_(this->_d) *
          (DataType_(1) + Math::cos(Math::pi<DataType_>() / DataType_(this->_m+1)));
      }
    }; // class PointstarFactoryFD

    /**
     * \brief Finite-Element pointstar matrix factory.
     *
     * This class generates a CSR poinstar matrices in Finite-Element style.
     * Moreover, this class can compute the largest and smallest eigenvalue of the
     * corresponding matrix as well as the eigenvector with respect to the smallest
     * eigenvalue.
     *
     * For now, this class can only compute 2D FE pointstars.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_ = Index>
    class PointstarFactoryFE :
      public PointstarFactoryBase<DataType_, IndexType_>
    {
    public:
      /**
       * \brief Creates the pointstar factory.
       *
       * \param[in] m
       * The number of inner nodes per dimension. Must be >= 3.
       */
      PointstarFactoryFE(Index m) :
        PointstarFactoryBase<DataType_, IndexType_>(m, Index(2))
      {
      }

      virtual ~PointstarFactoryFE() {}

      /**
       * \brief Generates a FE-style pointstar CSR matrix
       *
       * This function generates pointstar matrices in Finite-Element style.
       *
       * \returns
       * The m^d x m^d FE-stye pointstar matrix.
       */
      virtual SparseMatrixCSR<Mem::Main, DataType_, IndexType_> matrix_csr() const override
      {
        // declare some variables
        Index i,j,k,l,n,neq,nnze,bpl,rpb,rps,b;
        IndexType_ *epr,*dor,*ior,*col_idx;
        DataType_ *data;
        const Index m(this->_m), d(this->_d);

        /* calculate number of equations and non-zero entries */
        neq = m;
        nnze = j = 3*m - 2;
        for(i = 1; i < d; i++)
        {
          nnze *= j;
          neq *= m;
        }

        // create five vectors
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_row_ptr(neq+IndexType_(1));
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_col_idx(nnze, IndexType_(neq));
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_dor(neq, IndexType_(0));
        DenseVector<Mem::Main, IndexType_, IndexType_> vec_epr(neq, IndexType_(1));
        DenseVector<Mem::Main, DataType_, IndexType_> vec_data(nnze);

        // get data
        data = vec_data.elements();

        // get epr
        epr = vec_epr.elements();

        // get dor
        dor = vec_dor.elements();

        // get ior
        ior = vec_row_ptr.elements();

        // get column indices
        col_idx = vec_col_idx.elements();

        /* In this step we will recursively calculate the number of non-zero
         * elements in a row (= epr), and the number of non-zero elements left of
         * the main diagonal element in a row (= dor).
         * The indices of the first element of a row (= ior) and the indices of
         * the main diagonal elements (= dor) will later be recalculated from
         * epr and dor.
         */

        /* blocks per level */
        bpl = neq / m;

        /* rows per block */
        rpb = m;

        /* go through all levels */
        for(l = 0; l < d; ++l)
        {
          /* rows per sub-block */
          rps = rpb / m;

          /* go through all blocks */
          for(b = 0; b < bpl; ++b)
          {
            /* index of b-th block */
            k = rpb * b;

            /* add the right diagonal to sub-block 1 */
            for(j = 0; j < rps; ++j)
              epr[k+j] *= 2;

            /* add left and right diagonals to sub-blocks 2..m-1 */
            for(i = 1; i < m - 1; ++i)
            {
              /* index of i-th sub-block */
              n = k + rps*i;
              for(j = 0; j < rps; ++j)
              {
                dor[n+j] += epr[n+j];
                epr[n+j] *= 3;
              }
            }

            /* add left diagonal to sub-block m */
            n = k + rps*(m - 1);
            for(j = 0; j < rps; ++j)
            {
              dor[n+j] += epr[n+j];
              epr[n+j] *= 2;
            }
          }

          /* decrease blocks per level */
          bpl /= m;

          /* increase rows per block */
          rpb *= m;
        }

        /* initialize row and diagonal pointer arrays */
        /* calculate row and diagonal indices */
        ior[0] = 0;
        for(i = 0; i < neq; ++i)
        {
          ior[i+1] = ior[i] + epr[i];
          dor[i] += ior[i];
        }

        /* initialize data array */
        j = Index(3);
        for(i = 1; i < d; ++i)
          j *= Index (3);
        const DataType_ v1 = DataType_(j - Index(1));
        const DataType_ v2 = DataType_(-1);

        for(i = 0; i < neq; ++i)
        {
          /* left of main diagonal */
          for(j = ior[i]; j < dor[i]; ++j)
            data[j] = v2;

          /* main diagonal */
          data[dor[i]] = v1;

          /* right of main diagonal */
          for(j = dor[i]+1; j < ior[i+1]; ++j)
            data[j] = v2;
        }

        /* First, we initialize all indices to lie on the main diagonal.
         * At the same time, we will reset epr to 1.
         */
        for(i = 0; i < neq; ++i)
        {
          epr[i] = 1;
          for(j = ior[i]; j < ior[i+1]; ++j)
            col_idx[j] = IndexType_(neq+i);
        }

        /* blocks per level */
        bpl = neq / m;

        /* rows per block */
        rpb = m;

        /* go through all levels */
        for(l = d; l > 0; --l)
        {
          /* rows per sub-block */
          rps = rpb / m;

          /* go through all blocks */
          for(b = 0; b < bpl; ++b)
          {
            /* index of b-th block */
            k = rpb * b;

            /* go through the first sub-block */
            for(i = 0; i < rps; ++i)
            {
              for(n = 0; 2*n*epr[k+i] < ior[k+i+1]-ior[k+i]; ++n)
              {
                for(j = 0; j < epr[k+i]; ++j)
                {
                  col_idx[ior[k+i]+2*n*epr[k+i]+j+epr[k+i]] += IndexType_(rps);
                  col_idx[ior[neq-k-i]-2*n*epr[k+i]-j-epr[k+i]-1] -= IndexType_(rps);
                }
              }
              epr[k+i] *= 2;
              epr[k+rpb-i-1] *= 2;
            }

            /* go through the next m-2 sub-blocks */
            for(i = rps; i < rpb - rps; ++i)
            {
              for(n = 0; 3*n*epr[k+i] < ior[k+i+1]-ior[k+i]; ++n)
              {
                for(j = 0; j < epr[k+i]; ++j)
                {
                  col_idx[ior[k+i]+3*n*epr[k+i]+j+2*epr[k+i]] += IndexType_(rps);
                  col_idx[ior[neq-k-i]-3*n*epr[k+i]-j-2*epr[k+i]-1] -= IndexType_(rps);
                }
              }
              epr[k+i] *= 3;
            }
          }

          /* decrease blocks per level */
          bpl /= m;

          /* increase rows per block */
          rpb *= m;
        }

        // finally, update col_idx
        for(i = 0; i < nnze; ++i)
          col_idx[i] -= IndexType_(neq);

        // return the matrix
        return SparseMatrixCSR<Mem::Main, DataType_, IndexType_>(neq, neq, vec_col_idx, vec_data, vec_row_ptr);
      }

      /**
       * \brief Generates a FE-style pointstar CSR matrix with Neumann boundaries
       *
       * This function generates pointstar matrices in Finite-Element style using
       * Neumann boundaries.
       *
       * \returns
       * The m^d x m^d FE-stye pointstar matrix with Neumann boundaries.
       */
      virtual SparseMatrixCSR<Mem::Main, DataType_, IndexType_> matrix_csr_neumann() const
      {
        // create matrix with default values
        SparseMatrixCSR<Mem::Main, DataType_, IndexType_> matrix = this->matrix_csr();

        // get matrix arrays
        const IndexType_ m = IndexType_(this->_m);
        const IndexType_* row_ptr = matrix.row_ptr();
        const IndexType_* col_idx = matrix.col_ind();
        DataType_* val = matrix.val();

        // process first vertex line
        {
          // process first vertex in first line
          {
            const IndexType_ row = 0;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              const IndexType_ col = col_idx[k];
              if(col == row)
                val[k] = DataType_(4);
              else if(col == row+1)
                val[k] = DataType_(-1);
              else if(col == row+m)
                val[k] = DataType_(-1);
              else if(col == row+m+1)
                val[k] = DataType_(-2);
              else
                throw 0;
            }
          }

          // loop over all inner vertices in first line
          for(IndexType_ jj(1); jj + 1 < m; ++jj)
          {
            // loop over all non-zeros
            const IndexType_ row = jj;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              const IndexType_ col = col_idx[k];
              if(col == row)
                val[k] = DataType_(8);
              else if(col+1 == row)
                val[k] = DataType_(-1);
              else if(col == row+1)
                val[k] = DataType_(-1);
              else
                val[k] = DataType_(-2);
            }
          }

          // process last vertex in first line
          {
            const IndexType_ row =  m - 1;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              const IndexType_ col = col_idx[k];
              if(col == row)
                val[k] = DataType_(4);
              else if(col+1 == row)
                val[k] = DataType_(-1);
              else if(col == row+m)
                val[k] = DataType_(-1);
              else if(col+1 == row+m)
                val[k] = DataType_(-2);
              else
                throw 0;
            }
          }
        }

        // loop over all inner vertex lines
        for(IndexType_ ii(1); ii + 1 < m; ++ii)
        {
          // process first vertex in line ii
          {
            const IndexType_ row = ii*m;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              const IndexType_ col = col_idx[k];
              if(col == row)
                val[k] = DataType_(8);
              else if(col+m == row)
                val[k] = DataType_(-1);
              else if(col == row+m)
                val[k] = DataType_(-1);
              else
                val[k] = DataType_(-2);
            }
          }

          // loop over all inner vertices in line ii
          for(IndexType_ jj(1); jj + 1 < m; ++jj)
          {
            // loop over all non-zeros
            const IndexType_ row = ii*m + jj;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row+1]; ++k)
              val[k] = DataType_(col_idx[k] == row ? 16 : -2);
          }

          // process last vertex in line ii
          {
            const IndexType_ row = ii*m + m - 1;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              const IndexType_ col = col_idx[k];
              if(col == row)
                val[k] = DataType_(8);
              else if(col+m == row)
                val[k] = DataType_(-1);
              else if(col == row+m)
                val[k] = DataType_(-1);
              else
                val[k] = DataType_(-2);
            }
          }
        }


        // process last vertex line
        {
          // process first vertex in last line
          {
            const IndexType_ row = (m-1)*m;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              const IndexType_ col = col_idx[k];
              if(col == row)
                val[k] = DataType_(4);
              else if(col == row+1)
                val[k] = DataType_(-1);
              else if(col+m == row)
                val[k] = DataType_(-1);
              else if(col+m == row+1)
                val[k] = DataType_(-2);
              else
                throw 0;
            }
          }

          // loop over all inner vertices in last line
          for(IndexType_ jj(1); jj + 1 < m; ++jj)
          {
            // loop over all non-zeros
            const IndexType_ row = (m-1)*m + jj;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              const IndexType_ col = col_idx[k];
              if(col == row)
                val[k] = DataType_(8);
              else if(col+1 == row)
                val[k] = DataType_(-1);
              else if(col == row+1)
                val[k] = DataType_(-1);
              else
                val[k] = DataType_(-2);
            }
          }

          // process last vertex in last line
          {
            const IndexType_ row = (m-1)*m + m - 1;
            for(IndexType_ k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              const IndexType_ col = col_idx[k];
              if(col == row)
                val[k] = DataType_(4);
              else if(col+1 == row)
                val[k] = DataType_(-1);
              else if(col+m == row)
                val[k] = DataType_(-1);
              else if(col+m+1 == row)
                val[k] = DataType_(-2);
              else
                throw 0;
            }
          }
        }

        // return modified matrix
        return matrix;
      }


      /**
       * \brief Computes the smallest eigenvalue of the FE-style matrix.
       */
      virtual DataType_ lambda_min() const override
      {
        const DataType_ c(Math::cos(Math::pi<DataType_>() / DataType_(this->_m+1)));

        // lambda_min = 3^d - 1 - \sum_{k=1}^{d-1} 2^{d-k} {d \choose k} \cos\Big(\frac{\pi}{m+1}\Big)^k

        Index v(Index(1)); // v := 2^d
        Index a(Index(1)); // a := 3^d-1
        for(Index i(0); i < this->_d; ++i)
        {
          v *= Index(2);
          a *= Index(3);
        }

        DataType_ y = DataType_(v);
        for(Index k(1); k < this->_d; ++k)
        {
          (y *= c) += ((v *= this->_d - k + Index(1)) /= (Index(2) * k));
        }

        return DataType_(a - Index(1)) - c*y;
      }

      /**
       * \brief Computes the largest eigenvalue of the FE-style matrix.
       */
      virtual DataType_ lambda_max() const override
      {
        return DataType_(8) + DataType_(4)*Math::sqr(Math::cos(Math::pi<DataType_>() / DataType_(this->_m+1)));
      }
    }; // class PointstarFactory


    /**
     * \brief Pointstar factory base class
     *
     * \author Christoph Lohmann
     */
    template<typename DataType_, typename IndexType_ = Index>
    class PointstarFactoryBase2
    {
    protected:
      /// vector with number of subintervalls per dimension (plus a leading 2).
      const std::vector<IndexType_> & _num_of_subintervalls;
      /// vector with lengths of the n-dimensional hypercube
      const std::vector<DataType_> & _dimensions;
      /// number of dimensions
      const Index _d;

      PointstarFactoryBase2(const std::vector<IndexType_> & num_of_subintervalls,
                            const std::vector<DataType_> & dimensions) :
        _num_of_subintervalls(num_of_subintervalls),
        _dimensions(dimensions),
        _d(Index(_dimensions.size()))
      {
        XASSERTM(_d >= Index(1), "You need at least 1 dimension");
        XASSERTM(_d + 1 == _num_of_subintervalls.size(), "Vector-sizes do not match (n_1 = n_2 + 1)!");

        for (Index i(1); i <= _d; ++i)
        {
          XASSERTM(_num_of_subintervalls[i] >= Index(4), "You need at least 4 subintervalls per dimension");
        }
      }

      virtual ~PointstarFactoryBase2() {}

    public:
      /**
       * \brief Computes the pointstar matrix in banded format.
       *
       * \returns
       * The m^d x m^d pointstar matrix
       */
      virtual SparseMatrixBanded<Mem::Main, DataType_, IndexType_> matrix_banded() const = 0;

      /**
       * \brief Computes the smallest eigenvalue of the pointstar matrix.
       *
       * \returns
       * The smallest eigenvalue of the matrix.
       */
      virtual DataType_ lambda_min() const = 0;

      /**
       * \brief Computes the largest eigenvalue of the pointstar matrix.
       *
       * \returns
       * The largest eigenvalue of the matrix.
       */
      virtual DataType_ lambda_max() const = 0;

      /**
       * \brief Computes the spectral condition number of the pointstar matrix.
       *
       * \returns
       * The spectral condition number of the matrix.
       */
      virtual DataType_ spectral_cond() const
      {
        return lambda_max() / lambda_min();
      }

      /**
       * \brief Computes the eigenvector of the smallest eigenvalue
       *
       * This function generates the eigenvector with respect to the smallest eigenvalue
       * of the corresponding pointstar matrix.
       *
       * \returns
       * The m^d x m^d eigenvector
       */
      virtual DenseVector<Mem::Main, DataType_, IndexType_> eigenvector_min() const = 0;

      /**
       * \brief Computes a Q2-bubble vector.
       *
       * This function can be used to generate a useful solution vector for the pointstar matrix.
       * The corresponding right hand side can be computed by multiplying this vector by the matrix.
       *
       * \returns
       * The m^d x m^d Q2-bubble vector
       */
      virtual DenseVector<Mem::Main, DataType_, IndexType_> vector_q2_bubble() const = 0;
    }; // PointstarFactoryBase2

    /**
     * \brief Finite-Differences pointstar matrix factory.
     *
     * This class generates poinstar matrices in Finite-Differences style.
     * Moreover, this class can compute the largest and smallest eigenvalue of the
     * corresponding matrix as well as the eigenvector with respect to the smallest
     * eigenvalue.
     *
     * \author Christoph Lohmann
     */
    template<typename DataType_, typename IndexType_ = Index>
    class PointstarFactoryFD2 :
      public PointstarFactoryBase2<DataType_, IndexType_>
    {
    public:
      /**
       * \brief Creates the pointstar factory.
       *
       * \param[in] num_of_subintervalls
       * The vector with number of subintervalls per dimension plus a leading 2.
       * \param[in] dimensions
       * The vector with lengths of the n-dimensional hypercube
       */
      PointstarFactoryFD2(const std::vector<IndexType_> & num_of_subintervalls,
                          const std::vector<DataType_> & dimensions) :
        PointstarFactoryBase2<DataType_, IndexType_>(num_of_subintervalls, dimensions)
      {
      }

      virtual ~PointstarFactoryFD2() {}

      /**
       * \brief Generates a FD-style pointstar banded matrix
       *
       * \returns
       * The m^d x m^d FD-stye pointstar matrix.
       */
      virtual SparseMatrixBanded<Mem::Main, DataType_, IndexType_> matrix_banded() const override
      {
        const IndexType_ d(IndexType_(this->_d));
        const IndexType_ * const pnos(this->_num_of_subintervalls.data());
        const DataType_ * const pdim(this->_dimensions.data());

        /**
         * Create matrix-structure
         */
        SparseMatrixBanded<Mem::Main, DataType_, IndexType_> matrix(PointstarStructureFD::value<DataType_, IndexType_>(this->_num_of_subintervalls));
        const IndexType_ neq(IndexType_(matrix.rows()));
        DataType_ * const pval(matrix.val());

        /**
         * Fill matrix
         */
        // calculate diagonal entries
        DataType_ diagonal_entry(0);
        for (IndexType_ i(0); i < d; ++i)
        {
          diagonal_entry += DataType_(2.0) * Math::sqr(pdim[i] / DataType_(pnos[i + 1]));
        }

        // Fill diagonal with entries
        for (IndexType_ i(0); i < neq; ++i)
        {
          pval[neq * d + i] = diagonal_entry;
        }

        // Fill subdiagonals with entries
        for (IndexType_ i(0), k(pnos[0] - 1); i < d; ++i, k *= pnos[i] - 1)
        {
          const DataType_ subdiagonal_entry(-Math::sqr(pdim[i] / DataType_(pnos[i + 1])));

          for (IndexType_ j(0); j < neq / k / (pnos[i + 1] - 1); ++j)
          {
            const IndexType_ k1(neq * d + k * (pnos[i + 1] - 1) * j);
            const IndexType_ k2(neq * (i + 1));

            if (j != 0)
            {
              for (IndexType_ l(pnos[i] - 1); l > 0; --l)
              {
                pval[k1 + k2 - l] = DataType_(0);
                pval[k1 - k2 - l + k] = DataType_(0);
              }
            }

            for (IndexType_ l(0); l < k * (pnos[i + 1] - 2); ++l)
            {
              pval[k1 + l + k2] = subdiagonal_entry;
              pval[k1 + l - k2 + k] = subdiagonal_entry;
            }
          }
        }

        // return the matrix
        return matrix;
      }

      /**
       * \brief Computes the smallest eigenvalue of the FD-style matrix.
       *
       * The smallest eigenvalue of the FD-style matrix wrt the sine-bubble eigenvector is given as
       *
       *   \f[ \lambda_{\min} = 2 \cdot \sum^d_{k=1} h_k^2 \big(1 - \cos(\pi h_k) \big) \f]
       *
       * \returns
       * The smallest eigenvalue of the matrix.
       */
      virtual DataType_ lambda_min() const override
      {
        DataType_ x(DataType_(0.0));
        const IndexType_ * const pnos(this->_num_of_subintervalls.data());
        const DataType_ * const pdim(this->_dimensions.data());

        for (Index i(0); i < this->_d; ++i)
        {
          const DataType_ h(pdim[i] / DataType_(pnos[i + 1]));
          x += Math::sqr(h) * (DataType_(1.0) - Math::cos(Math::pi<DataType_>() * h / pdim[i]));
        }

        return 2 * x;
      }

      /**
       * \brief Computes the largest eigenvalue of the FD-style matrix.
       *
       * The largest eigenvalue of the FD-style matrix is given as
       *
       * \f[ \lambda_{\max} = 2 \cdot \sum^d_{k=1} h_k^2 \big(1 + \cos(\pi h_k) \big) \f]
       *
       * \returns
       * The largest eigenvalue of the matrix.
       */
      virtual DataType_ lambda_max() const override
      {
        DataType_ x(DataType_(0.0));
        const IndexType_ * const pnos(this->_num_of_subintervalls.data());
        const DataType_ * const pdim(this->_dimensions.data());

        for (Index i(0); i < this->_d; ++i)
        {
          const DataType_ h(pdim[i] / DataType_(pnos[i + 1]));
          x += Math::sqr(h) * (DataType_(1.0) + Math::cos(Math::pi<DataType_>() * h / pdim[i]));
        }

        return 2 * x;
      }

      virtual DenseVector<Mem::Main, DataType_, IndexType_> eigenvector_min() const override
      {
        const IndexType_ * const pnos(this->_num_of_subintervalls.data());
        const Index d(this->_d);

        // compute vector length
        Index size(1);
        for (Index i(1); i <= d; ++i)
        {
          size *= pnos[i] - 1;
        }

        // create vector
        DenseVector<Mem::Main, DataType_, IndexType_> vector(size, DataType_(1));
        DataType_* v = vector.elements();

        for(Index i(0); i < size; ++i)
        {
          for(Index j(0), quo(1); j < d; ++j, quo *= pnos[j] - 1)
          {
            // compute x scaling factor
            const DataType_ xsf(Math::pi<DataType_>() / DataType_(pnos[j + 1]));
            v[i] *= Math::sin(xsf * DataType_(((i / quo) % (pnos[j + 1] - 1)) + Index(1)));
          }
        }

        // return vector
        return std::move(vector);
      }

      /**
       * \brief Computes a Q2-bubble vector.
       *
       * This function can be used to generate a useful solution vector for the pointstar matrix.
       * The corresponding right hand side can be computed by multiplying this vector by the matrix.
       *
       * \returns
       * The m^d x m^d Q2-bubble vector
       */
      DenseVector<Mem::Main, DataType_, IndexType_> vector_q2_bubble() const override
      {
        const Index d(this->_d);
        const IndexType_ * const pnos(this->_num_of_subintervalls.data());
        const DataType_ * const pdim(this->_dimensions.data());

        // compute vector length
        Index size(1);
        for (Index i(1); i <= d; ++i)
        {
          size *= pnos[i] - 1;
        }

        // create vector
        DenseVector<Mem::Main, DataType_, IndexType_> vector(size, DataType_(1));
        DataType_* v = vector.elements();

        for(Index i(0); i < size; ++i)
        {
          for(Index j(0), quo(1); j < d; ++j, quo *= pnos[j] - 1)
          {
            // compute x scaling factor
            const DataType_ xsf(pdim[j] / DataType_(pnos[j + 1]));

            // compute coordinate x in (0,pdim[j])
            DataType_ x(xsf * DataType_(((i / quo) % (pnos[j + 1] - 1)) + Index(1)));

            v[i] *= DataType_(4) * x * (DataType_(1) - x);
          }
        }

        // return vector
        return std::move(vector);
      }
    }; // class PointstarFactoryFD2

  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_POINTSTAR_FACTORY_HPP
