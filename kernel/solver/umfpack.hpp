// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_UMFPACK_HPP
#define KERNEL_SOLVER_UMFPACK_HPP 1

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
#if defined(FEAT_HAVE_UMFPACK) || defined(DOXYGEN)
    /**
     * \brief UMFPACK solver class
     *
     * This class provides an implementation of the SolverBase interface using the
     * direct solver UMFPACK for doing the actual dirty work.
     *
     * \note
     * This solver can only be applied onto SparseMatrixCSR<Mem::Main,double,Index> matrices.
     * If you want to apply UMFPACK on other matrix types, use the GenericUmfpack solver instead.
     *
     * \attention
     * This class is only declared if FEAT was configured to build and link against
     * the \c UMFPACK third-party library.
     *
     * \author Peter Zajac
     */
    class Umfpack :
      public SolverBase<LAFEM::DenseVector<Mem::Main, double, Index>>
    {
    public:
      /// compatible matrix type
      typedef LAFEM::SparseMatrixCSR<Mem::Main, double, Index> MatrixType;
      /// compatible vector type
      typedef LAFEM::DenseVector<Mem::Main, double, Index> VectorType;

      /// our base class
      typedef SolverBase<VectorType> BaseClass;

    private:
      /// system matrix
      const MatrixType& _system_matrix;
      /// umfpack control array
      double* _umf_control;
      /// umfpack symbolic factorization pointer
      void* _umf_symbolic;
      /// umfpack numeric factorization pointer
      void* _umf_numeric;

      /// symbolic peak memory size
      std::size_t _sym_peak_size;
      /// symbolic factorization memory size
      std::size_t _sym_mem_size;
      /// numeric factorization memory size
      std::size_t _num_mem_size;
      /// total peak memory size
      std::size_t _umf_peak_size;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] system_matrix
       * A reference to the system matrix to be factorized.
       */
      explicit Umfpack(const MatrixType& system_matrix);

      /// virtual destructor
      virtual ~Umfpack();

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Umfpack";
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
    }; // class Umfpack

    /**
     * \brief Creates a new Umfpack solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \returns
     * A shared pointer to a new Umfpack object.
     */
    inline std::shared_ptr<Umfpack> new_umfpack(const LAFEM::SparseMatrixCSR<Mem::Main, double, Index>& matrix)
    {
      return std::make_shared<Umfpack>(matrix);
    }

    /**
     * \brief UMFPACK Mean solver class
     *
     * This class implements a variant of the Umfpack solver, which is capable of solving linear
     * systems including the integral mean constraint or any other constraint that can be expressed
     * as a scalar Lagrange multiplier to the original system matrix.
     *
     * As Umfpack is a direct solver, it cannot directly utilise the MeanFilter class, which is used
     * by iterative solvers to enforce the integral mean constraint. Therefore, this class implements
     * an algorithm, which extends the original system matrix by a Lagrange multiplier vector and
     * applied the Umfpack solver onto this extended linear system.
     *
     * This class offers two constructors:
     * - A CTOR which takes a SparseMatrixCSR and a MeanFilter as input
     * - A CTOR which takes a SparseMatrixCSR and a DenseVector representing the Lagrange multiplier
     *   as input.
     *
     * \author Peter Zajac
     */
    class UmfpackMean :
      public SolverBase<LAFEM::DenseVector<Mem::Main, double, Index>>
    {
    public:
      /// compatible matrix type
      typedef LAFEM::SparseMatrixCSR<Mem::Main, double, Index> MatrixType;
      /// compatible vector type
      typedef LAFEM::DenseVector<Mem::Main, double, Index> VectorType;

      /// our base class
      typedef SolverBase<VectorType> BaseClass;

    private:
      /// system matrix
      const MatrixType& _system_matrix;
      /// weight vector
      const VectorType& _weight_vector;
      /// our extended system matrix
      MatrixType _solver_matrix;
      /// two temporary extended vectors
      VectorType _vec_x, _vec_b;
      /// our internal Umfpack solver object
      Umfpack _umfpack;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] system_matrix
       * A reference to the system matrix to be factorized.
       *
       * \param[in] weight_vector
       * The weight vector to be used as a Lagrange multiplier.
       */
      explicit UmfpackMean(const MatrixType& system_matrix, const VectorType& weight_vector);

      /**
       * \brief Constructor
       *
       * \param[in] system_matrix
       * A reference to the system matrix to be factorized.
       *
       * \param[in] mean_filter
       * A reference to the mean filter containing the weight vector.
       */
      explicit UmfpackMean(
        const MatrixType& system_matrix,
        const LAFEM::MeanFilter<Mem::Main, double, Index>& mean_filter) :
        UmfpackMean(system_matrix, mean_filter.get_vec_dual())
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "UmfpackMean";
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
    }; // class UmfpackMean

    /**
     * \brief Creates a new UmfpackMean solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] weight_vector
     * The weight vector to be used as a Lagrange multiplier.
     *
     * \returns
     * A shared pointer to a new UmfpackMean object.
     */
    inline std::shared_ptr<UmfpackMean> new_umfpack_mean(
      const LAFEM::SparseMatrixCSR<Mem::Main, double, Index>& matrix,
      const LAFEM::DenseVector<Mem::Main, double, Index>& weight_vector)
    {
      return std::make_shared<UmfpackMean>(matrix, weight_vector);
    }

    /**
     * \brief Creates a new UmfpackMean solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * A reference to the mean filter containing the weight vector.
     *
     * \returns
     * A shared pointer to a new UmfpackMean object.
     */
    inline std::shared_ptr<UmfpackMean> new_umfpack_mean(
      const LAFEM::SparseMatrixCSR<Mem::Main, double, Index>& matrix,
      const LAFEM::MeanFilter<Mem::Main, double, Index>& filter)
    {
      return std::make_shared<UmfpackMean>(matrix, filter);
    }

    /**
     * \brief UMFPACK Saddle-Point Mean solver class
     *
     * This class implements a variant of the Umfpack solver, which is capable of solving linear
     * systems on SaddlePointMatrix objects including the integral mean constraint or any other
     * constraint that can be expressed as a scalar Lagrange multiplier to the second solution
     * component (usually the pressure) of the original system.
     *
     * This class can be used to solve (Navier-)Stokes systems with full Dirichlet boundary
     * conditions for the velocity, which in turn requires a MeanFilter constraint for the
     * pressure to make the system regular. In this case, the sub-matrix blocks A and B
     * are required to be filtered by the corresponding UnitFilterBlocked by the caller,
     * whereas the MeanFilter is passed to this solver class.
     *
     * This class offers two constructors:
     * - A CTOR which takes a SaddlePointMatrix and a MeanFilter for the second component as input
     * - A CTOR which takes a SaddlePointMatrix and a DenseVector representing the Lagrange multiplier
     *   as input.
     *
     * \author Peter Zajac
     */
    template<typename DT_, typename IT_, int dim_>
    class SaddleUmfpackMean :
      public SolverBase<LAFEM::TupleVector<LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_>, LAFEM::DenseVector<Mem::Main, DT_, IT_>>>
    {
    public:
      /// compatible matrix type
      typedef LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, dim_> MatrixTypeA;
      typedef LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, 1> MatrixTypeB;
      typedef LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, 1, dim_> MatrixTypeD;
      typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> MatrixType;

      /// compatible vector type
      typedef LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_> VectorTypeV;
      typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorTypeP;
      typedef LAFEM::TupleVector<VectorTypeV, VectorTypeP> VectorType;

      typedef LAFEM::SparseMatrixCSR<Mem::Main, double, Index> UmfMatrixType;
      typedef LAFEM::DenseVector<Mem::Main, double, Index> UmfVectorType;

      /// our base class
      typedef SolverBase<VectorType> BaseClass;

    private:
      /// system matrix
      const MatrixType& _system_matrix;
      /// weight vector from mean filter
      const VectorTypeP& _weight_vector;
      /// our extended system matrix
      UmfMatrixType _solver_matrix;
      /// two temporary extended vectors
      UmfVectorType _vec_x, _vec_b;
      /// our internal Umfpack solver object
      Umfpack _umfpack;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] system_matrix
       * A reference to the system matrix to be factorized.
       *
       * \param[in] weight_vector
       * The weight vector to be used as a Lagrange multiplier.
       */
      explicit SaddleUmfpackMean(const MatrixType& system_matrix, const VectorTypeP& weight_vector) :
        BaseClass(),
        _system_matrix(system_matrix),
        _weight_vector(weight_vector),
        _umfpack(_solver_matrix)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] system_matrix
       * A reference to the system matrix to be factorized.
       *
       * \param[in] mean_filter
       * A reference to the pressure mean filter containing the weight vector.
       */
      explicit SaddleUmfpackMean(
        const MatrixType& system_matrix,
        const LAFEM::MeanFilter<Mem::Main, DT_, IT_>& mean_filter) :
        SaddleUmfpackMean(system_matrix, mean_filter.get_vec_dual())
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "SaddleUmfpackMean";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // get sub-matrices
        const MatrixTypeA& matrix_a = _system_matrix.block_a();
        const MatrixTypeB& matrix_b = _system_matrix.block_b();
        const MatrixTypeD& matrix_d = _system_matrix.block_d();

        // get velocity/pressure space dimensions
        const Index n_v = matrix_b.rows();
        const Index n_p = matrix_b.columns();

        // verify matrix dimensions
        XASSERT(n_v == matrix_a.rows());
        XASSERT(n_v == matrix_a.columns());
        XASSERT(n_p == matrix_d.rows());
        XASSERT(n_v == matrix_d.columns());
        XASSERTM(n_p == _weight_vector.size(), "invalid weight vector/mean filter size");

        // get the number of non-zeroes in sub-matrices
        const Index nnze_a = matrix_a.used_elements();
        const Index nnze_b = matrix_b.used_elements();
        const Index nnze_d = matrix_d.used_elements();

        // allocate our solver matrix
        _solver_matrix = UmfMatrixType(dim_*n_v+n_p+1, dim_*n_v+n_p+1, dim_*dim_*nnze_a + dim_*(nnze_b+nnze_d) + 2*n_p);

        // get our input matrix arrays
        const IT_* irow_ptr_a = matrix_a.row_ptr();
        const IT_* icol_idx_a = matrix_a.col_ind();
        const IT_* irow_ptr_b = matrix_b.row_ptr();
        const IT_* icol_idx_b = matrix_b.col_ind();
        const IT_* irow_ptr_d = matrix_d.row_ptr();
        const IT_* icol_idx_d = matrix_d.col_ind();

        // get our output matrix arrays
        Index* orow_ptr = _solver_matrix.row_ptr();
        Index* ocol_idx = _solver_matrix.col_ind();

        // assemble the solver matrix structure
        orow_ptr[0] = Index(0);
        // loop over all rows of A and B
        for(Index i(0); i < n_v; ++i)
        {
          Index id = i*Index(dim_);
          // loop over all block dimensions
          for(Index j(0); j < Index(dim_); ++j, ++id)
          {
            Index op = orow_ptr[id];
            // copy input matrix A row pattern
            for(Index ip(irow_ptr_a[i]); ip < irow_ptr_a[i+1]; ++ip)
            {
              // convert block to scalar indices
              for(Index k(0); k < Index(dim_); ++k, ++op)
                ocol_idx[op] = Index(icol_idx_a[ip])*Index(dim_) + k;
            }
            // copy input matrix B row pattern
            for(Index ip(irow_ptr_b[i]); ip < irow_ptr_b[i+1]; ++ip, ++op)
              ocol_idx[op] = n_v*Index(dim_) + Index(icol_idx_b[ip]);
            // set next row pointer
            orow_ptr[id+1] = op;
          }
        }
        // loop over all rows of D
        for(Index i(0); i < n_p; ++i)
        {
          Index op = orow_ptr[Index(dim_)*n_v + i];
          // copy input matrix D row pattern
          for(Index ip(irow_ptr_d[i]); ip < irow_ptr_d[i+1]; ++ip)
          {
            // convert block to scalar value
            for(Index k(0); k < Index(dim_); ++k, ++op)
              ocol_idx[op] = Index(icol_idx_d[ip])*Index(dim_) + k;
          }
          // append a single entry in the last column
          ocol_idx[op] = Index(dim_)*n_v + n_p;
          // set next row pointer
          orow_ptr[Index(dim_)*n_v+i+1] = ++op;
        }

        // append an (almost) dense row for pressure
        Index op(orow_ptr[Index(dim_)*n_v + n_p]);
        for(Index j(0); j < n_p; ++j, ++op)
        {
          ocol_idx[op] = Index(dim_)*n_v + j;
        }
        // set last row pointer
        orow_ptr[Index(dim_)*n_v+n_p+1] = op;

        // Okay, solver matrix structure assembly complete

        // create temporary vectors
        _vec_x = _solver_matrix.create_vector_r();
        _vec_b = _solver_matrix.create_vector_r();

        // initialize umfpack
        _umfpack.init_symbolic();
      }

      virtual void done_symbolic() override
      {
        _umfpack.done_symbolic();
        _vec_b.clear();
        _vec_x.clear();
        _solver_matrix.clear();

        BaseClass::done_symbolic();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // get sub-matrices
        const MatrixTypeA& matrix_a = _system_matrix.block_a();
        const MatrixTypeB& matrix_b = _system_matrix.block_b();
        const MatrixTypeD& matrix_d = _system_matrix.block_d();

        // get velocity/pressure space dimensions
        const Index n_v = matrix_b.rows();
        const Index n_p = matrix_b.columns();

        // get our input matrix arrays
        const IT_* irow_ptr_a = matrix_a.row_ptr();
        const auto* idata_a   = matrix_a.val();
        const IT_* irow_ptr_b = matrix_b.row_ptr();
        const auto* idata_b   = matrix_b.val();
        const IT_* irow_ptr_d = matrix_d.row_ptr();
        const auto* idata_d   = matrix_d.val();

        // get input vector array
        const DT_* weight = _weight_vector.elements();

        // get output matrix arrays
        const Index* orow_ptr = _solver_matrix.row_ptr();
        double* odata = _solver_matrix.val();

        // loop over all rows of A and B
        for(Index i(0); i < n_v; ++i)
        {
          Index id = i*Index(dim_);
          // loop over all block dimensions
          for(int j(0); j < dim_; ++j, ++id)
          {
            Index op = orow_ptr[id];
            // copy input matrix A row data
            for(Index ip(irow_ptr_a[i]); ip < irow_ptr_a[i+1]; ++ip)
            {
              for(int k(0); k < dim_; ++k, ++op)
                odata[op] = double(idata_a[ip][j][k]);
            }
            // copy input matrix B row data
            for(Index ip(irow_ptr_b[i]); ip < irow_ptr_b[i+1]; ++ip, ++op)
              odata[op] = double(idata_b[ip][j][0]);
          }
        }
        // loop over all rows of D
        for(Index i(0); i < n_p; ++i)
        {
          Index op = orow_ptr[Index(dim_)*n_v + i];
          // copy input matrix D row data
          for(Index ip(irow_ptr_d[i]); ip < irow_ptr_d[i+1]; ++ip)
          {
            // convert block to scalar value
            for(Index k(0); k < Index(dim_); ++k, ++op)
              odata[op] = double(idata_d[ip][0][k]);
          }
          // copy weight vector entry as last column entry
          odata[op] = double(weight[i]);
        }

        // copy weight vector into last row
        Index op(orow_ptr[Index(dim_)*n_v + n_p]);
        for(Index j(0); j < n_p; ++j, ++op)
        {
          odata[op] = double(weight[j]);
        }

        // okay, solver matrix data assembly completed

        // initialize umfpack
        _umfpack.init_numeric();
      }

      virtual void done_numeric() override
      {
        _umfpack.done_numeric();

        BaseClass::done_numeric();
      }

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
      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        const Index n_v = vec_sol.template at<0>().size() * Index(dim_);
        const Index n_p = vec_sol.template at<1>().size();

        // copy RHS vector
        const DT_* vr_v = vec_rhs.template at<0>().template elements<LAFEM::Perspective::pod>();
        const DT_* vr_p = vec_rhs.template at<1>().elements();
        double* vb = _vec_b.elements();

        // copy velocity RHS
        for(Index i(0); i < n_v; ++i, ++vb)
          *vb = double(vr_v[i]);
        // copy pressure RHS
        for(Index i(0); i < n_p; ++i, ++vb)
          *vb = double(vr_p[i]);
        // append Lagrange multiplier RHS
        *vb = 0.0;

        // solve system
        Status status = _umfpack.apply(_vec_x, _vec_b);

        // copy sol vector
        const double* vx = _vec_x.elements();
        DT_* vs_v = vec_sol.template at<0>().template elements<LAFEM::Perspective::pod>();
        DT_* vs_p = vec_sol.template at<1>().elements();
        // copy velocity solution
        for(Index i(0); i < n_v; ++i, ++vx)
          vs_v[i] = DT_(*vx);
        // copy pressure solution
        for(Index i(0); i < n_p; ++i, ++vx)
          vs_p[i] = DT_(*vx);

        // okay
        return status;
      }
    }; // class SaddleUmfpackMean

    /**
     * \brief Creates a new UmfpackMean solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] weight_vector
     * The weight vector to be used as a Lagrange multiplier.
     *
     * \returns
     * A shared pointer to a new UmfpackMean object.
     */
    template<typename DT_, typename IT_, int dim_>
    inline std::shared_ptr<SaddleUmfpackMean<DT_, IT_, dim_>> new_saddle_umfpack_mean(
      const LAFEM::SaddlePointMatrix<
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, dim_>,
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, 1>,
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, 1, dim_>>& matrix,
      const LAFEM::DenseVector<Mem::Main, DT_, IT_>& weight_vector)
    {
      return std::make_shared<SaddleUmfpackMean<DT_, IT_, dim_>>(matrix, weight_vector);
    }

    /**
     * \brief Creates a new UmfpackMean solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * A reference to the mean filter containing the weight vector.
     *
     * \returns
     * A shared pointer to a new UmfpackMean object.
     */
    template<typename DT_, typename IT_, int dim_>
    inline std::shared_ptr<SaddleUmfpackMean<DT_, IT_, dim_>> new_saddle_umfpack_mean(
      const LAFEM::SaddlePointMatrix<
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, dim_>,
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, 1>,
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, 1, dim_>>& matrix,
      const LAFEM::MeanFilter<Mem::Main, DT_, IT_>& filter)
    {
      return std::make_shared<SaddleUmfpackMean<DT_, IT_, dim_>>(matrix, filter);
    }

    /// \cond internal
    namespace Intern
    {
      // This helper is class is required if the data/index type pair of
      // the system vector is not double/Index. In this case, the system
      // vector has to be converted to a DenseVector with the same data/index
      // type pair first and then converted to a double/Index dense vector
      // in a second step.
      template<typename DT_, typename IT_>
      class GenericUmfpackVectorHelper
      {
      private:
        LAFEM::DenseVector<Mem::Main, DT_, IT_> _vec_tmp;

      public:
        void init(const LAFEM::DenseVector<Mem::Main, double, Index>& v)
        {
          _vec_tmp.convert(v);
        }

        void done()
        {
          _vec_tmp.clear();
        }

        template<typename VT_>
        void download(LAFEM::DenseVector<Mem::Main, double, Index>& vo, const VT_& vi)
        {
          _vec_tmp.copy(vi);
          vo.convert(_vec_tmp);
        }

        template<typename VT_>
        void upload(VT_& vo, const LAFEM::DenseVector<Mem::Main, double, Index>& vi)
        {
          _vec_tmp.convert(vi);
          _vec_tmp.copy_inv(vo);
        }
      };

      // specialization for double/Index vectors; no intermediate conversion is necessary here
      template<>
      class GenericUmfpackVectorHelper<double, Index>
      {
      public:
        void init(const LAFEM::DenseVector<Mem::Main, double, Index>&)
        {
        }

        void done()
        {
        }

        template<typename VT_>
        void download(LAFEM::DenseVector<Mem::Main, double, Index>& vo, const VT_& vi)
        {
          vo.copy(vi);
        }

        template<typename VT_>
        void upload(VT_& vo, const LAFEM::DenseVector<Mem::Main, double, Index>& vi)
        {
          vi.copy_inv(vo);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Generic UMFPACK solver class
     *
     * This class effectively wraps around an Umfpack solver object
     *
     * \attention
     * This class is only declared if FEAT was configured to build and link against
     * the \c UMFPACK third-party library.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_>
    class GenericUmfpack :
      public SolverBase<typename Matrix_::VectorTypeL>
    {
    public:
      /// our base class
      typedef SolverBase<typename Matrix_::VectorTypeL> BaseClass;
      /// our matrix type
      typedef Matrix_ MatrixType;
      /// our vector type
      typedef typename MatrixType::VectorTypeL VectorType;

    protected:
      /// our matrix
      const MatrixType& _matrix;
      /// the matrix for our Umfpack solver (SparseMatrixCSR<Mem::Main, double, Index>)
      typename Umfpack::MatrixType _umf_matrix;
      /// the vectors for our Umfpack solver (DenseVector<Mem::Main, double, Index>)
      typename Umfpack::VectorType _umf_vsol, _umf_vrhs;
      /// vector helper
      Intern::GenericUmfpackVectorHelper<typename VectorType::DataType, typename VectorType::IndexType> _vec_helper;
      /// the actual Umfpack solver object
      Umfpack _umfpack;

    public:
      explicit GenericUmfpack(const MatrixType& matrix) :
        _matrix(matrix),
        _umfpack(_umf_matrix)
      {
      }

      /// virtual destructor
      virtual ~GenericUmfpack()
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "GenericUmfpack";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // convert matrix to obtain the structure
        _umf_matrix.convert(_matrix);

        // create vectors
        _umf_vsol = _umf_matrix.create_vector_l();
        _umf_vrhs = _umf_matrix.create_vector_l();

        // initialiase vector helper
        _vec_helper.init(_umf_vsol);

        // factorize symbolic
        _umfpack.init_symbolic();
      }

      virtual void done_symbolic() override
      {
        _umfpack.done_symbolic();

        _vec_helper.done();

        _umf_vrhs.clear();
        _umf_vsol.clear();
        _umf_matrix.clear();

        BaseClass::done_symbolic();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // convert our system matrix to <Mem::Main, double, Index> (if necessary)
        typename MatrixType::template ContainerType<Mem::Main, double, Index> mat_main;
        mat_main.convert(_matrix);

        // get the array of our CSR matrix
        const Index* row_ptr = _umf_matrix.row_ptr();
        /*const*/ Index* col_idx = _umf_matrix.col_ind();
        double* val = _umf_matrix.val();

        // copy entries into our CSR matrix
        for(Index i(0); i < _umf_matrix.rows(); ++i)
          mat_main.set_line(i, val + row_ptr[i], col_idx + row_ptr[i], 0);

        // factorize
        _umfpack.init_numeric();
      }

      virtual void done_numeric() override
      {
        _umfpack.done_numeric();
        BaseClass::done_numeric();
      }

      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // convert RHS vector
        _vec_helper.download(_umf_vrhs, vec_rhs);

        // solve
        Status status = _umfpack.apply(_umf_vsol, _umf_vrhs);

        // convert sol vector
        _vec_helper.upload(vec_sol, _umf_vsol);

        return status;
      }
    }; // class GenericUmfpack<...>

    /**
     * \brief Creates a new GenericUmfpack solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \returns
     * A shared pointer to a new GenericUmfpack object.
     */
    template<typename Matrix_>
    inline std::shared_ptr<GenericUmfpack<Matrix_>> new_generic_umfpack(const Matrix_& matrix)
    {
      return std::make_shared<GenericUmfpack<Matrix_>> (matrix);
    }

#endif // defined(FEAT_HAVE_UMFPACK) || defined(DOXYGEN)
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_UMFPACK_HPP
