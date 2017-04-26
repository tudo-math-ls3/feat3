#pragma once
#ifndef KERNEL_SOLVER_SPAI_PRECOND_HPP
#define KERNEL_SOLVER_SPAI_PRECOND_HPP 1

#define FEAT3_SPAI_DEBUG_OUTPUT 0

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/math.hpp>
#include <kernel/solver/legacy_preconditioners.hpp>

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
#ifdef FEAT_HAVE_CUDA
      // forward declaration of GPU functions
      unsigned int cuda_spai_maxrowlen(const unsigned int m, const unsigned int *rowptr, unsigned int *rowlens);

      template <typename DT_>
        void cuda_spai_construct_rows(const unsigned int this_chunk_size, const unsigned int first_row_in_chunk,
            const DT_ *A_val, const unsigned int *A_colind, const unsigned int *A_rowptr, DT_ *M_val,
            const unsigned int *rowlens, unsigned int maxrowlen);
#endif

#ifdef FEAT_HAVE_MKL
      // forward declaration of MKL functions
      // solve minimization problem using mkl qr decomposition
      // Ahat is stored in column-major order
      void spai_mkl_solve_minimization_problem(float *Ahat, unsigned long lenJk, unsigned long rowlen, unsigned long mrl, float *rhs, float *sol);
      void spai_mkl_solve_minimization_problem(double *Ahat, unsigned long lenJk, unsigned long rowlen, unsigned long mrl, double *rhs, double *sol);
#endif

      /**
       * \brief Factory Class for the creation of a SPAI-1 inverse
       *
       * Supports CSR Format only
       * Supports Men::Main in combination with DataTypes float and double and IndexType unsigned long
       * Supports Men::CUDA in combination with DataTypes float and double and IndexType unsigned int
       */
      template<class MatrixType>
        class SPAIFactory;

#ifdef FEAT_HAVE_CUDA
      // SPAI-1 Factory for float/double and Mem::CUDA
      template<typename DT_>
        class SPAIFactory<LAFEM::SparseMatrixCSR<Mem::CUDA, DT_, unsigned int>>
        {
          using MatrixType = LAFEM::SparseMatrixCSR<Mem::CUDA, DT_, unsigned int>;

          public:
          static MatrixType construct_M(const MatrixType& A, Index chunksize_in)
          {
            TimeStamp tstart;
            unsigned int chunksize = (unsigned int)chunksize_in;

            const unsigned int m = (unsigned int)A.rows();

            // inverse shares layout with A
            MatrixType M;
            M.clone(A, LAFEM::CloneMode::Layout);

            // determine maximum row length (globally!)
            LAFEM::DenseVector<Mem::CUDA, unsigned int, unsigned int> rowlens(m);
            unsigned int maxrowlen = cuda_spai_maxrowlen(m, A.row_ptr(), rowlens.elements());

            // distribution in chunks
            if (chunksize == 0 || chunksize > m) chunksize = m;
            unsigned int nchunks = 0, lastchunksize = 0;
            // all chunks are equal
            if (m % chunksize == 0)
            {
              nchunks = m / chunksize;
              lastchunksize = chunksize;
            }
            // last one is smaller than the rest
            else
            {
              nchunks = m / chunksize + 1;
              lastchunksize = m - (nchunks - 1)*chunksize;
            }

            // for each chunk
            for (unsigned i = 0; i < nchunks; i++)
            {
              const unsigned this_chunk_size = (i == nchunks - 1) ? lastchunksize : chunksize;
              const unsigned first_row_in_chunk = i*chunksize;

              // construct rows
              cuda_spai_construct_rows(this_chunk_size, first_row_in_chunk,
                  A.val(), A.col_ind(), A.row_ptr(), M.val(),
                  rowlens.elements(), maxrowlen);
            }

            TimeStamp tdone;
#       if FEAT3_SPAI_DEBUG_OUTPUT > 0
            std::cout << "Elapsed construction time: " << tdone.elapsed(tstart) << " s" << std::endl;
#       endif
#       if FEAT3_SPAI_DEBUG_OUTPUT > 2
            M.write_out_mtx(std::cout);
#       endif

            return M;
          }
        };
#endif

#ifdef FEAT_HAVE_MKL
      // SPAI-1 factory for float/double and Mem::Main - MKL
      template<typename DT_>
        class SPAIFactory<LAFEM::SparseMatrixCSR<Mem::Main, DT_, unsigned long>>
        {
          using MatrixType = LAFEM::SparseMatrixCSR<Mem::Main, DT_, unsigned long>;

          public:
          static MatrixType construct_M(const MatrixType& A, Index /*chunksize*/)
          {
            // measure construction time
            TimeStamp tstart;

            const unsigned long m = A.rows();

            // inverse shares layout with A
            MatrixType M;
            M.clone(A, LAFEM::CloneMode::Layout);

            // raw data
            const unsigned long * const A_rowptr = A.row_ptr();
            const unsigned long * const A_colind = A.col_ind();
            const DT_ * const A_val = A.val();
            DT_ * const M_val = M.val();

            // calculate maxrowlength
            unsigned long mrl = 0;
            for(unsigned long i=0; i<m; i++)
            {
              unsigned long tmp = A_rowptr[i+1] - A_rowptr[i];
              if(tmp > mrl) mrl = tmp;
            }

            for(unsigned long k=0; k<m; k++)
            {
              // length of the current row
              const unsigned long rowlen = A_rowptr[k+1] - A_rowptr[k];

              // Jk is the set of rows that are connected to row k via its nnz
              LAFEM::DenseVector<Mem::Main, unsigned long, unsigned long> Jk(rowlen*mrl, (unsigned long)(0));
              unsigned long * const Jk_raw = Jk.elements();
              // Ahat is the reduced matrix A to Ik, Jk
              // Storage is column major
              LAFEM::DenseVector<Mem::Main, DT_, unsigned long> Ahat(rowlen*rowlen*mrl, DT_(0));
              DT_ * const Ahat_raw = Ahat.elements();

              // construct Jk (without duplicates) and fill Ahat at the same time
              unsigned long lenJk = 0;
              for(unsigned long i=0; i < rowlen; i++)
              {
                const unsigned long ii = A_colind[A_rowptr[k] + i];
                const unsigned long lenii = A_rowptr[ii+1] - A_rowptr[ii];
                for(unsigned long j=0; j < lenii; j++)
                {
                  const unsigned long jj = A_colind[A_rowptr[ii] + j];
                  // already in matrix?
                  for (unsigned long l=0; l < lenJk; l++)
                  {
                    if (Jk_raw[l] == jj)
                    {
                      // insert at this position. length of a full column is rowlen*mrl
                      Ahat_raw[i*rowlen*mrl + l] = A_val[A_rowptr[ii] + j];
                      goto nextj; // break inner loop and skip rest of outer loop
                    }
                  }
                  // not yet in matrix -> new line
                  Jk_raw[lenJk] = jj;
                  Ahat_raw[i*rowlen*mrl + lenJk] = A_val[A_rowptr[ii] + j];
                  lenJk++;
nextj:
                  continue;
                }
              }

              // construct rhs of the minimization problem
              LAFEM::DenseVector<Mem::Main, DT_, unsigned long> rhs(lenJk, DT_(0));
              DT_ * const rhs_raw = rhs.elements();
              for(unsigned long i=0; i < lenJk; i++)
              {
                if(Jk_raw[i] == k)
                {
                  // there are no duplicates in Jk
                  rhs_raw[i] = DT_(1);
                  break;
                }
              }

              // perform actual minimization
              spai_mkl_solve_minimization_problem(Ahat_raw, lenJk, rowlen, mrl, rhs_raw, M_val + A_rowptr[k]);
            }

            TimeStamp tdone;
#       if FEAT3_SPAI_DEBUG_OUTPUT > 0
            std::cout << "Elapsed construction time: " << tdone.elapsed(tstart) << " s" << std::endl;
#       endif
#       if FEAT3_SPAI_DEBUG_OUTPUT > 2
            M.write_out_mtx(std::cout);
#       endif

            return M;
          }
        };
#endif

    } // namespace Intern


    template<typename Matrix_, typename Filter_>
    class SPAIPrecond;

    /// general spai implementation without any optimisations
    template<typename DT_, typename IT_, typename Filter_>
      class SPAIPrecond<LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>, Filter_> : public SolverBase<typename LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>::VectorTypeL>
    {
      public:
        typedef LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_> MatrixType;
        typedef typename MatrixType::VectorTypeR VectorType;

        /// the filter object
        const Filter_& _filter;
        /// the actual preconditioner object
        SPAIPreconditioner<MatrixType, VectorType> _precond;

        template<typename... Args_>
          explicit SPAIPrecond(const MatrixType& matrix, const Filter_& filter) :
            _filter(filter),
            _precond(matrix)
      {
      }

        /// Returns the name of the solver.
        virtual String name() const override
        {
          return _precond.name();
        }

        /// Applies the preconditioner.
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          TimeStamp ts_start;
          _precond.apply(vec_cor, vec_def);
          this->_filter.filter_cor(vec_cor);
          TimeStamp ts_stop;
          Statistics::add_time_precon(ts_stop.elapsed(ts_start));
          return Status::success;
        }
    };

    /// general spai implementation without any optimisations
    template<typename DT_, typename IT_, typename Filter_>
      class SPAIPrecond<LAFEM::SparseMatrixELL<Mem::Main, DT_, IT_>, Filter_> : public SolverBase<typename LAFEM::SparseMatrixELL<Mem::Main, DT_, IT_>::VectorTypeL>
    {
      public:
        typedef LAFEM::SparseMatrixELL<Mem::Main, DT_, IT_> MatrixType;
        typedef typename MatrixType::VectorTypeR VectorType;

        /// the filter object
        const Filter_& _filter;
        /// the actual preconditioner object
        SPAIPreconditioner<MatrixType, VectorType> _precond;

        template<typename... Args_>
          explicit SPAIPrecond(const MatrixType& matrix, const Filter_& filter) :
            _filter(filter),
            _precond(matrix)
      {
      }

        /// Returns the name of the solver.
        virtual String name() const override
        {
          return _precond.name();
        }

        /// Applies the preconditioner.
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          TimeStamp ts_start;
          _precond.apply(vec_cor, vec_def);
          this->_filter.filter_cor(vec_cor);
          TimeStamp ts_stop;
          Statistics::add_time_precon(ts_stop.elapsed(ts_start));
          return Status::success;
        }
    };

#ifdef FEAT_HAVE_MKL
    /**
     * \brief SPAI(1) preconditioner implementation
     *
     * This class implements a sparse approximate inverse preconditioner
     * which approximates the inverse by a matrix with the same layout.
     *
     * So far, this implementation only works for the matrix type
     * LAFEM::SparseMatrixCSR.
     *
     * This implementation supports Mem::CUDA and Mem::Main architecture.
     *
     * MKL support is required in order to use this class.
     *
     * \author David Schneider
     */
    template<typename DT_,  typename Filter_>
      class SPAIPrecond<LAFEM::SparseMatrixCSR<Mem::Main, DT_, unsigned long>, Filter_> :
      public SolverBase<typename LAFEM::SparseMatrixCSR<Mem::Main, DT_, unsigned long>::VectorTypeL>
      {
        public:
          using Matrix_ = LAFEM::SparseMatrixCSR<Mem::Main, DT_, unsigned long>;
          typedef Matrix_ MatrixType;
          typedef Filter_ FilterType;
          typedef typename MatrixType::VectorTypeL VectorType;
          typedef typename MatrixType::DataType DataType;
          typedef typename MatrixType::IndexType IndexType;

        protected:
          const MatrixType& _matrix;
          const FilterType& _filter;
          Index _chunksize;
          MatrixType _M;

        public:
          /**
           * \brief Constructor.
           *
           * \param[in] matrix
           * The system matrix.
           *
           * \param[in] filter
           * The system filter.
           *
           * \praram[in] chunksize
           * For Mem::CUDA, maximum number of rows of M that are constructed in parallel.
           * The optimal value of this parameter depends on the available GPU RAM.
           * When set to zero (default), all rows are computed at once.
           * For Mem::Main, this parameter has no effect.
           */
          SPAIPrecond(
              const MatrixType& matrix, const FilterType& filter, Index chunksize = 0) :
            _matrix(matrix),
            _filter(filter),
            _chunksize(chunksize)
        {
        }

          /// Returns the name of the solver.
          virtual String name() const override
          {
            return "SPAI_1";
          }

          virtual void init_symbolic() override
          {
            // check if construction of SPAI is possible
            if (_matrix.columns() != _matrix.rows())
            {
              throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square.");
            }
          }

          virtual void done_symbolic() override
          {
            _M.clear();
          }

          virtual void init_numeric() override
          {
            _M = Intern::SPAIFactory<MatrixType>::construct_M(_matrix, _chunksize);
          }

          virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
          {
            // multiply with the SPAI
            _M.apply(vec_cor, vec_def);
            // apply filter
            this->_filter.filter_cor(vec_cor);

            return Status::success;
          }

          const MatrixType& get_spai() const
          {
            return _M;
          }

      }; // class SPAIPrecond<...>
#endif

#ifdef FEAT_HAVE_CUDA
    template<typename DT_, typename Filter_>
      class SPAIPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, DT_, unsigned int>, Filter_> :
      public SolverBase<typename LAFEM::SparseMatrixCSR<Mem::CUDA, DT_, unsigned int>::VectorTypeL>
      {
        public:
          using Matrix_ = LAFEM::SparseMatrixCSR<Mem::CUDA, DT_, unsigned int>;
          typedef Matrix_ MatrixType;
          typedef Filter_ FilterType;
          typedef typename MatrixType::VectorTypeL VectorType;
          typedef typename MatrixType::DataType DataType;
          typedef typename MatrixType::IndexType IndexType;

        protected:
          const MatrixType& _matrix;
          const FilterType& _filter;
          Index _chunksize;
          MatrixType _M;

        public:
          /**
           * \brief Constructor.
           *
           * \param[in] matrix
           * The system matrix.
           *
           * \param[in] filter
           * The system filter.
           *
           * \praram[in] chunksize
           * For Mem::CUDA, maximum number of rows of M that are constructed in parallel.
           * The optimal value of this parameter depends on the available GPU RAM.
           * When set to zero (default), all rows are computed at once.
           * For Mem::Main, this parameter has no effect.
           */
          SPAIPrecond(
              const MatrixType& matrix, const FilterType& filter, Index chunksize = 0) :
            _matrix(matrix),
            _filter(filter),
            _chunksize(chunksize)
        {
        }

          /// Returns the name of the solver.
          virtual String name() const override
          {
            return "SPAI_1";
          }

          virtual void init_symbolic() override
          {
            // check if construction of SPAI is possible
            if (_matrix.columns() != _matrix.rows())
            {
              throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square.");
            }
          }

          virtual void done_symbolic() override
          {
            _M.clear();
          }

          virtual void init_numeric() override
          {
            _M = Intern::SPAIFactory<MatrixType>::construct_M(_matrix, _chunksize);
          }

          virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
          {
            // multiply with the SPAI
            _M.apply(vec_cor, vec_def);
            // apply filter
            this->_filter.filter_cor(vec_cor);

            return Status::success;
          }

          const MatrixType& get_spai() const
          {
            return _M;
          }

      }; // class SPAIPrecond<...>
#endif

    /// SPAIPrecond specialisation for saddle point matrices
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename Filter_>
    class SPAIPrecond<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>, Filter_> :
      public SolverBase<LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL>>
      {
        public:
        /// our matrix type
        typedef LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_> MatrixType;
        /// our filter type
        typedef Filter_ FilterType;

        /// our vector type
        typedef typename MatrixType::VectorTypeR VectorType;
        /// our data type
        typedef typename MatrixType::DataType DataType;
        /// our index type
        typedef typename MatrixType::IndexType IndexType;
        /// our index type
        typedef typename MatrixType::MemType MemType;

        protected:
        const MatrixType& _matrix;
        const FilterType& _filter;
        MatrixType _M;

        public:

        SPAIPrecond(
            const MatrixType& matrix, const FilterType& filter) :
          _matrix(matrix),
          _filter(filter)
        {
        }

        /// Returns the name of the solver.
        virtual String name() const override
        {
          return "SPAI_1";
        }

        virtual void init_symbolic() override
        {
        }

        virtual void done_symbolic() override
        {
          _M.clear();
        }

        virtual void init_numeric() override
        {
          LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> csr;
          csr.convert(_matrix);
          auto csr_result = Intern::SPAIFactory<decltype(csr)>::construct_M(csr, 0);
          _M = _matrix.clone(LAFEM::CloneMode::Layout);
          csr_result.convert_reverse(_M);
        }

        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // multiply with the SPAI
          _M.apply(vec_cor, vec_def);
          // apply filter
          this->_filter.filter_cor(vec_cor);

          return Status::success;
        }

        const MatrixType& get_spai() const
        {
          return _M;
        }

      }; // class SPAIPrecond<...>

    /// Dummy class for not implemented specialisations
    template<typename Matrix_, typename Filter_>
    class SPAIPrecond :
      public SolverBase<typename Matrix_::VectorTypeL>
    {
      public:

      explicit SPAIPrecond(const Matrix_&, const Filter_&)
      {
      }

      Status apply(typename Matrix_::VectorTypeL &, const typename Matrix_::VectorTypeL &) override
      {
          throw InternalError(__func__, __FILE__, __LINE__, "not implemented yet!");
      }

      String name() const override
      {
          throw InternalError(__func__, __FILE__, __LINE__, "not implemented yet!");
      }
    };

    /**
     * \brief Creates a new SPAIPrecond solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new SPAIPrecond object.
     */
    template<typename Matrix_, typename Filter_, typename... Args_>
      inline std::shared_ptr<SPAIPrecond<Matrix_, Filter_>> new_spai_precond(
          const Matrix_& matrix, const Filter_& filter, Args_&&... args)
      {
        return std::make_shared<SPAIPrecond<Matrix_, Filter_>>
          (matrix, filter, std::forward<Args_>(args)...);
      }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_SPAI_PRECOND_HPP
