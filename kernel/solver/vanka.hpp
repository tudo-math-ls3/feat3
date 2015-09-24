#pragma once
#ifndef KERNEL_SOLVER_VANKA_HPP
#define KERNEL_SOLVER_VANKA_HPP 1

// includes, FEAST
#include <kernel/solver/base.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_csr_blocked.hpp>
#include <kernel/lafem/power_diag_matrix.hpp>
#include <kernel/lafem/power_full_matrix.hpp>
#include <kernel/lafem/power_row_matrix.hpp>
#include <kernel/lafem/power_col_matrix.hpp>
#include <kernel/lafem/power_full_matrix.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>

// includes, system
#include <map>
#include <set>
#include <vector>

namespace FEAST
{
  namespace Solver
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Internal Vanka Vector helper class template
       *
       * This class template is used by Vanka to perform 2 tasks:
       * - Gather the local defect vector entries.
       * - Scatter the local correction vector entries.
       *
       * Currently, the following specialisations exist:
       * - LAFEM::DenseVector<Mem::Main,...>
       * - LAFEM::PowerVector<...>
       *
       * \author Peter Zajac
       */
      template<typename Vector_>
      class VankaVector;

      /**
       * \brief Internal Vanka Matrix helper class template
       *
       * this class template is used by Vanka to perform 2 tasks:
       * - Gather the local matrix entries (full or diagonal only)
       * - Compute the local defect vector entries.
       *
       * Currently, the following specialisations exist:
       * - LAFEM::SparseMatrixCSR<Mem::Main,...>
       * - LAFEM::PowerDiagMatrix<...>
       * - LAFEM::PowerFullMatrix<...>
       * - LAFEM::PowerRowMatrix<...>
       * - LAFEM::PowerColMatrix<...>
       *
       * \author Peter Zajac
       */
      template<typename Matrix_>
      class VankaMatrix;

      template<typename DT_, typename IT_>
      class VankaVector<LAFEM::DenseVector<Mem::Main, DT_, IT_>>
      {
      public:
        typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorType;
        static constexpr int dim = 1;

      protected:
        DT_* _vec_cor;

      public:
        explicit VankaVector(VectorType& vec_cor) :
          _vec_cor(vec_cor.elements())
        {
        }

        IT_ gather_def(DT_* x, const VectorType& vec, const IT_* idx, const IT_ n, const IT_ off) const
        {
          const DT_* vdef = vec.elements();

          for(IT_ i(0); i < n; ++i)
          {
            x[off+i] = vdef[idx[i]];
          }
          return off+n;
        }

        IT_ scatter_cor(const DT_ omega, const DT_* x, const IT_* idx, const IT_ n, const IT_ off)
        {
          for(IT_ i(0); i < n; ++i)
          {
            _vec_cor[idx[i]] += omega * x[off+i];
          }
          return off+n;
        }
      };

      template<typename SubVector_, int dim_>
      class VankaVector<LAFEM::PowerVector<SubVector_, dim_>>
      {
      public:
        typedef LAFEM::PowerVector<SubVector_, dim_> VectorType;
        typedef VankaVector<SubVector_> FirstClass;
        typedef VankaVector<LAFEM::PowerVector<SubVector_, dim_-1>> RestClass;

        static constexpr int dim = FirstClass::dim + RestClass::dim;

      protected:
        FirstClass _first;
        RestClass _rest;

      public:
        explicit VankaVector(VectorType& vec_cor) :
          _first(vec_cor.first()),
          _rest(vec_cor.rest())
        {
        }

        template<typename DT_, typename IT_>
        IT_ gather_def(DT_* x, const VectorType& vec, const IT_* idx, const IT_ n, const IT_ off) const
        {
          IT_ noff = _first.gather_def(x, vec.first(), idx, n, off);
          return _rest.gather_def(x, vec.rest(), idx, n, noff);
        }

        template<typename DT_, typename IT_>
        IT_ scatter_cor(const DT_ omega, const DT_* x, const IT_* idx, const IT_ n, const IT_ off)
        {
          IT_ noff = _first.scatter_cor(omega, x, idx, n, off);
          return _rest.scatter_cor(omega, x, idx, n, noff);
        }
      };

      template<typename SubVector_>
      class VankaVector<LAFEM::PowerVector<SubVector_, 1>>
      {
      public:
        typedef LAFEM::PowerVector<SubVector_, 1> VectorType;
        typedef VankaVector<SubVector_> FirstClass;

        static constexpr int dim = FirstClass::dim;

      protected:
        FirstClass _first;

      public:
        explicit VankaVector(VectorType& vec_cor) :
          _first(vec_cor.first())
        {
        }

        template<typename DT_, typename IT_>
        IT_ gather_def(DT_* x, const VectorType& vec, const IT_* idx, const IT_ n, const IT_ off) const
        {
          return _first.gather_def(x, vec.first(), idx, n, off);
        }

        template<typename DT_, typename IT_>
        IT_ scatter_cor(const DT_ omega, const DT_* x, const IT_* idx, const IT_ n, const IT_ off)
        {
          return _first.scatter_cor(omega, x, idx, n, off);
        }
      };

      template<typename DT_, typename IT_>
      class VankaMatrix<LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>>
      {
      public:
        typedef LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_> MatrixType;
        typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorTypeR;

        static constexpr int row_dim = 1;
        static constexpr int col_dim = 1;

      protected:
        const IT_* _row_ptr;
        const IT_* _col_idx;
        const DT_* _mat_val;

      public:
        explicit VankaMatrix(const MatrixType& matrix) :
          _row_ptr(matrix.row_ptr()),
          _col_idx(matrix.col_ind()),
          _mat_val(matrix.val())
        {
        }

        std::pair<IT_, IT_> gather_full(
          DT_* data, const IT_* ridx, const IT_* cidx,
          const IT_ m, const IT_ n, const IT_ stride,
          const IT_ mo = IT_(0), const IT_ no = IT_(0)) const
        {
          // empty matrix?
          if(_mat_val == nullptr)
            return std::make_pair(mo+m, no+n);

          // loop over all local rows
          for(IT_ i(0); i < m; ++i)
          {
            // get row index
            const IT_ ri = ridx[i];

            // initialise loop variable for local columns
            IT_ j = IT_(0);

            // loop over the ri row of our matrix
            for(IT_ ai = _row_ptr[ri]; ai < _row_ptr[ri+1]; ++ai)
            {
              // get column index
              const IT_ rj = _col_idx[ai];

              // now search its position in our local matrix
              while((j < n) && (cidx[j] < rj))
              {
                ++j;
              }
              // did we find our local entry?
              if((j < n) && (cidx[j] == rj))
              {
                // found, so store it in our local matrix
                data[(mo + i) * stride + no + j] = _mat_val[ai];
              }
            }
          }
          return std::make_pair(mo+m, no+n);
        }

        IT_ gather_diag(DT_* data, const IT_* idx, const IT_ m, const IT_ mo = IT_(0)) const
        {
          // loop over all local rows
          for(IT_ i(0); i < m; ++i)
          {
            // get row index
            const IT_ ri = idx[i];

            // find diagonal entry
            for(IT_ ai = _row_ptr[ri]; ai < _row_ptr[ri+1]; ++ai)
            {
              // is it our diagonal?
              if(_col_idx[ai] == ri)
              {
                data[mo+i] = _mat_val[ai];
                break;
              }
            }
          }
          return mo + m;
        }

        IT_ mult_cor(DT_* x, const DT_ alpha, const VectorTypeR& vec_cor, const IT_* idx, const IT_ m, const IT_ off) const
        {
          if(_mat_val == nullptr)
            return off + m;

          const DT_* v = vec_cor.elements();

          // loop over all local rows
          for(IT_ i(0); i < m; ++i)
          {
            const IT_ vi = idx[i];

            // loop over all columns
            DT_ r = DT_(0);
            for(IT_ k = _row_ptr[vi]; k < _row_ptr[vi+1]; ++k)
            {
              r += _mat_val[k] * v[_col_idx[k]];
            }
            x[off+i] += alpha * r;
          }
          return off + m;
        }
      };

      template<typename SubMatrix_, int dim_>
      class VankaMatrix<LAFEM::PowerDiagMatrix<SubMatrix_, dim_>>
      {
      public:
        typedef LAFEM::PowerDiagMatrix<SubMatrix_, dim_> MatrixType;
        typedef typename MatrixType::VectorTypeR VectorTypeR;
        typedef VankaMatrix<SubMatrix_> FirstClass;
        typedef VankaMatrix<LAFEM::PowerDiagMatrix<SubMatrix_, dim_-1>> RestClass;

        static constexpr int row_dim = FirstClass::row_dim + RestClass::row_dim;
        static constexpr int col_dim = FirstClass::col_dim + RestClass::col_dim;

      protected:
        FirstClass _first;
        RestClass _rest;

      public:
        explicit VankaMatrix(const MatrixType& matrix) :
          _first(matrix.first()),
          _rest(matrix.rest())
        {
        }

        template<typename DT_, typename IT_>
        std::pair<IT_, IT_> gather_full(
          DT_* data, const IT_* ridx, const IT_* cidx,
          const IT_ m, const IT_ n, const IT_ stride,
          const IT_ mo = IT_(0), const IT_ no = IT_(0)) const
        {
          std::pair<IT_, IT_> mn = _first.gather_full(data, ridx, cidx, m, n, stride, mo, no);
          return _rest.gather_full(data, ridx, cidx, m, n, stride, mn.first, mn.second);
        }

        template<typename DT_, typename IT_>
        IT_ gather_diag(DT_* data, const IT_* idx, const IT_ m, const IT_ mo = IT_(0)) const
        {
          IT_ mno = _first.gather_diag(data, idx, m, mo);
          return _rest.gather_diag(data, idx, m, mno);
        }

        template<typename DT_, typename IT_>
        IT_ mult_cor(DT_* x, const DT_ alpha, const VectorTypeR& vec_cor, const IT_* idx, const IT_ m, const IT_ off) const
        {
          IT_ noff = _first.mult_cor(x, alpha, vec_cor.first(), idx, m, off);
          return _rest.mult_cor(x, alpha, vec_cor.rest(), idx, m, noff);
        }
      };

      template<typename SubMatrix_>
      class VankaMatrix<LAFEM::PowerDiagMatrix<SubMatrix_, 1>>
      {
      public:
        typedef LAFEM::PowerDiagMatrix<SubMatrix_, 1> MatrixType;
        typedef typename MatrixType::VectorTypeR VectorTypeR;
        typedef VankaMatrix<SubMatrix_> FirstClass;

        static constexpr int row_dim = FirstClass::row_dim;
        static constexpr int col_dim = FirstClass::col_dim;

      protected:
        FirstClass _first;

      public:
        explicit VankaMatrix(const MatrixType& matrix) :
          _first(matrix.first())
        {
        }

        template<typename DT_, typename IT_>
        std::pair<IT_, IT_> gather_full(
          DT_* data, const IT_* ridx, const IT_* cidx,
          const IT_ m, const IT_ n, const IT_ stride,
          const IT_ mo = IT_(0), const IT_ no = IT_(0)) const
        {
          return _first.gather_full(data, ridx, cidx, m, n, stride, mo, no);
        }

        template<typename DT_, typename IT_>
        IT_ gather_diag(DT_* data, const IT_* idx, const IT_ m, const IT_ mo = IT_(0)) const
        {
          return _first.gather_diag(data, idx, m, mo);
        }

        template<typename DT_, typename IT_>
        IT_ mult_cor(DT_* x, const DT_ alpha, const VectorTypeR& vec_cor, const IT_* idx, const IT_ m, const IT_ off) const
        {
          return _first.mult_cor(x, alpha, vec_cor.first(), idx, m, off);
        }
      };

      template<typename SubMatrix_, int dim_>
      class VankaMatrix<LAFEM::PowerColMatrix<SubMatrix_, dim_>>
      {
      public:
        typedef LAFEM::PowerColMatrix<SubMatrix_, dim_> MatrixType;
        typedef typename MatrixType::VectorTypeR VectorTypeR;
        typedef VankaMatrix<SubMatrix_> FirstClass;
        typedef VankaMatrix<LAFEM::PowerColMatrix<SubMatrix_, dim_-1>> RestClass;

        static constexpr int row_dim = FirstClass::row_dim + RestClass::row_dim;
        static constexpr int col_dim = FirstClass::col_dim;

      protected:
        FirstClass _first;
        RestClass _rest;

      public:
        explicit VankaMatrix(const MatrixType& matrix) :
          _first(matrix.first()),
          _rest(matrix.rest())
        {
        }

        template<typename DT_, typename IT_>
        std::pair<IT_, IT_> gather_full(
          DT_* data, const IT_* ridx, const IT_* cidx,
          const IT_ m, const IT_ n, const IT_ stride,
          const IT_ mo = IT_(0), const IT_ no = IT_(0)) const
        {
          std::pair<IT_, IT_> mno = _first.gather_full(data, ridx, cidx, m, n, stride, mo, no);
          return _rest.gather_full(data, ridx, cidx, m, n, stride, mno.first, no);
        }

        template<int ro_, int co_, typename DT_, typename IT_>
        IT_ gather_diag_pfm(DT_* data, const IT_* idx, const IT_ m, const IT_ mo = IT_(0)) const
        {
          IT_ mno = _first.template gather_diag_pfm<ro_, co_>(data, idx, m, mo);
          return _rest.template gather_diag_pfm<ro_+1, co_>(data, idx, m, mno);
        }

        template<typename DT_, typename IT_>
        IT_ mult_cor(DT_* x, const DT_ alpha, const VectorTypeR& vec_cor, const IT_* idx, const IT_ m, const IT_ off) const
        {
          IT_ noff = _first.mult_cor(x, alpha, vec_cor, idx, m, off);
          return _rest.mult_cor(x, alpha, vec_cor, idx, m, noff);
        }
      };

      template<typename SubMatrix_>
      class VankaMatrix<LAFEM::PowerColMatrix<SubMatrix_, 1>>
      {
      public:
        typedef LAFEM::PowerColMatrix<SubMatrix_, 1> MatrixType;
        typedef typename MatrixType::VectorTypeR VectorTypeR;
        typedef VankaMatrix<SubMatrix_> FirstClass;

        static constexpr int row_dim = FirstClass::row_dim;
        static constexpr int col_dim = FirstClass::col_dim;

      protected:
        VankaMatrix<SubMatrix_> _first;

      public:
        explicit VankaMatrix(const MatrixType& matrix) :
          _first(matrix.first())
        {
        }

        template<typename DT_, typename IT_>
        std::pair<IT_, IT_> gather_full(
          DT_* data, const IT_* ridx, const IT_* cidx,
          const IT_ m, const IT_ n, const IT_ stride,
          const IT_ mo = IT_(0), const IT_ no = IT_(0)) const
        {
          return _first.gather_full(data, ridx, cidx, m, n, stride, mo, no);
        }

        template<int ro_, int co_, typename DT_, typename IT_>
        IT_ gather_diag_pfm(DT_* data, const IT_* idx, const IT_ m, const IT_ mo = IT_(0)) const
        {
          return _first.template gather_diag_pfm<ro_, co_>(data, idx, m, mo);
        }

        template<typename DT_, typename IT_>
        IT_ mult_cor(DT_* x, const DT_ alpha, const VectorTypeR& vec_cor, const IT_* idx, const IT_ m, const IT_ off) const
        {
          return _first.mult_cor(x, alpha, vec_cor, idx, m, off);
        }
      };

      template<typename SubMatrix_, int dim_>
      class VankaMatrix<LAFEM::PowerRowMatrix<SubMatrix_, dim_>>
      {
      public:
        typedef LAFEM::PowerRowMatrix<SubMatrix_, dim_> MatrixType;
        typedef typename MatrixType::VectorTypeR VectorTypeR;
        typedef VankaMatrix<SubMatrix_> FirstClass;
        typedef VankaMatrix<LAFEM::PowerRowMatrix<SubMatrix_, dim_-1>> RestClass;

        static constexpr int row_dim = FirstClass::row_dim;
        static constexpr int col_dim = FirstClass::col_dim + RestClass::col_dim;

      protected:
        FirstClass _first;
        RestClass _rest;

      public:
        explicit VankaMatrix(const MatrixType& matrix) :
          _first(matrix.first()),
          _rest(matrix.rest())
        {
        }

        template<typename DT_, typename IT_>
        std::pair<IT_, IT_> gather_full(
          DT_* data, const IT_* ridx, const IT_* cidx,
          const IT_ m, const IT_ n, const IT_ stride,
          const IT_ mo = IT_(0), const IT_ no = IT_(0)) const
        {
          std::pair<IT_, IT_> mno = _first.gather_full(data, ridx, cidx, m, n, stride, mo, no);
          return _rest.gather_full(data, ridx, cidx, m, n, stride, mo, mno.second);
        }

        template<int ro_, int co_, typename DT_, typename IT_>
        IT_ gather_diag_pfm(DT_* data, const IT_* idx, const IT_ m, const IT_ mo = IT_(0)) const
        {
          if(ro_ == co_)
            return _first.gather_diag(data, idx, m, mo);
          else
            return _rest.template gather_diag_pfm<ro_, co_+1>(data, idx, m, mo);
        }

        template<typename DT_, typename IT_>
        IT_ mult_cor(DT_* x, const DT_ alpha, const VectorTypeR& vec_cor, const IT_* idx, const IT_ m, const IT_ off) const
        {
          _first.mult_cor(x, alpha, vec_cor.first(), idx, m, off);
          return _rest.mult_cor(x, alpha, vec_cor.rest(), idx, m, off);
        }
      };

      template<typename SubMatrix_>
      class VankaMatrix<LAFEM::PowerRowMatrix<SubMatrix_, 1>>
      {
      public:
        typedef LAFEM::PowerRowMatrix<SubMatrix_, 1> MatrixType;
        typedef VankaMatrix<SubMatrix_> FirstClass;
        typedef typename MatrixType::VectorTypeR VectorTypeR;

        static constexpr int row_dim = FirstClass::row_dim;
        static constexpr int col_dim = FirstClass::col_dim;

      protected:
        FirstClass _first;

      public:
        explicit VankaMatrix(const MatrixType& matrix) :
          _first(matrix.first())
        {
        }

        template<typename DT_, typename IT_>
        std::pair<IT_, IT_> gather_full(
          DT_* data, const IT_* ridx, const IT_* cidx,
          const IT_ m, const IT_ n, const IT_ stride,
          const IT_ mo = IT_(0), const IT_ no = IT_(0)) const
        {
          return _first.gather_full(data, ridx, cidx, m, n, stride, mo, no);
        }

        template<int ro_, int co_, typename DT_, typename IT_>
        IT_ gather_diag_pfm(DT_* data, const IT_* idx, const IT_ m, const IT_ mo = IT_(0)) const
        {
          if(ro_ == co_)
            return _first.gather_diag(data, idx, m, mo);
          else
            return mo;
        }

        template<typename DT_, typename IT_>
        IT_ mult_cor(DT_* x, const DT_ alpha, const VectorTypeR& vec_cor, const IT_* idx, const IT_ m, const IT_ off) const
        {
          return _first.mult_cor(x, alpha, vec_cor.first(), idx, m, off);
        }
      };

      template<typename SubMatrix_, int dim_w_, int dim_h_>
      class VankaMatrix<LAFEM::PowerFullMatrix<SubMatrix_, dim_w_, dim_h_>>
      {
      public:
        typedef LAFEM::PowerFullMatrix<SubMatrix_, dim_w_, dim_h_> MatrixType;
        typedef typename MatrixType::VectorTypeR VectorTypeR;
        typedef VankaMatrix<typename MatrixType::ContClass> ContClass;

        static constexpr int row_dim = ContClass::row_dim;
        static constexpr int col_dim = ContClass::col_dim;

      protected:
        ContClass _cont;

      public:
        explicit VankaMatrix(const MatrixType& matrix) :
          _cont(matrix.get_container())
        {
        }

        template<typename DT_, typename IT_>
        std::pair<IT_, IT_> gather_full(
          DT_* data, const IT_* ridx, const IT_* cidx,
          const IT_ m, const IT_ n, const IT_ stride,
          const IT_ mo = IT_(0), const IT_ no = IT_(0)) const
        {
          return _cont.gather_full(data, ridx, cidx, m, n, stride, mo, no);
        }

        template<typename DT_, typename IT_>
        IT_ gather_diag(DT_* data, const IT_* idx, const IT_ m, const IT_ mo = IT_(0)) const
        {
          return _cont.template gather_diag_pfm<0,0>(data, idx, m, mo);
        }

        template<typename DT_, typename IT_>
        IT_ mult_cor(DT_* x, const DT_ alpha, const VectorTypeR& vec_cor, const IT_* idx, const IT_ m, const IT_ off) const
        {
          return _cont.mult_cor(x, alpha, vec_cor, idx, m, off);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Vanka type enumeration
     *
     * This enumeration specifies the various types of supported Vanka variants.
     * See the documentation of the Vanka solver class template for details.
     */
    enum class VankaType
    {
      /// Nodal diagonal Vanka type (multiplicative)
      nodal_diag_mult = 0x000,
      /// Nodal full Vanka type (multiplicative)
      nodal_full_mult = 0x001,
      /// Blocked diagonal Vanka type (multiplicative)
      block_diag_mult = 0x010,
      /// Blocked full Vanka type (multiplicative)
      block_full_mult = 0x011,
      /// Nodal diagonal Vanka type (additive)
      nodal_diag_add  = 0x100,
      /// Nodal full Vanka type (additive)
      nodal_full_add  = 0x101,
      /// Blocked diagonal Vanka type (additive)
      block_diag_add  = 0x110,
      /// Blocked full Vanka type (additive)
      block_full_add  = 0x111
    };

    /**
     * \brief Vanka Factorisation Error
     *
     * This exception is thrown by the Vanka preconditioner if the inversion
     * of a local system failed. This indicates that the system matrix may
     * (but does not necessarily need to) be singular.
     */
    class VankaFactorError :
      public SolverException
    {
    public:
      VankaFactorError() : SolverException("Vanka Factorisation Error") {}
    };

    /**
     * \brief Vanka preconditioner/smoother class template.
     *
     * This class template is only implemented for LAFEM::SaddlePointMatrix;
     * see the documentation of the corresponding specialisation for information.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class Vanka;

    /**
     * \brief Vanka preconditioner/smoother class template implementation
     *
     * This class template implements the Vanka preconditioner/smoother, which is a special
     * Block-SOR/Gauss-Seidel variant suitable for working with saddle-point systems of the form
     * \f[\begin{bmatrix} A & B\\D & 0\end{bmatrix} \cdot \begin{bmatrix} v\\p\end{bmatrix} = \begin{bmatrix} f\\g\end{bmatrix}\f]
     *
     * In total, this class template implements eight different variants of the Vanka, which can
     * be chosen at run-time by specifying the corresponding VankaType in the constructor of this
     * object. The eight variants are the cross-product of three binary variants:
     *
     * <b>1.a: Block Variants</b> (<c>VankaType::block_*_*</c>)\n
     * The block variants collect all pressure DOFs, which are adjacent to each other via the
     * connectivity graph of <em>D*B</em> with the same degree. If the \e B and \e D matrices
     * are build by using a <em>discontiunuous</em> pressure element, then these block variants
     * collect all DOFs associated with a single element of the underlying mesh, thus leading
     * to the 'true' element-based Vanka.\n
     * - Block variants are recommended for \e discontinuous pressure spaces.
     *
     * <b>1.b: Nodal Variants</b> (<c>VankaType::nodal_*_*</c>)\n
     * The nodal variants process exactly one pressure DOF per block.
     * - Nodal variants are recommended for \e continuous pressure spaces.
     *
     * <b>2.a: Diagonal Variants</b> (<c>VankaType::*_diag_*</c>)\n
     * These variants apply the Schur-complement approach onto each local system and
     * approximate the inverse of \e A by the inverse of its main diagonal. If there is
     * more than one local pressure DOF, the local Schur-Complement matrix itself is inverted
     * by using Gaussian elimination.
     *
     * <b>2.b: Full Variants</b> (<c>VankaType::*_full_*</c>)\n
     * These variants invert the full dense local systems by using Gaussian elimination.
     *
     * <b>3.a: Multiplicative Variants</b> (<c>VankaType::*_*_mult</c>)\n
     * These variants work successively within one Vanka iteration, just as the (scalar)
     * Gauss-Seidel- or SOR-iterations do.
     *
     * <b>3.b: Additive Variants</b> (<c>VankaType::*_*_add</c>)\n
     * These variants work in an additive way, i.e. they ignore the 'progess' made in the
     * current iteration, just as the (scalar) Jacobi-iterations do.
     *
     * Moreover, each of the above variants is parameterised in a relaxation/damping parameter
     * \e omega, which can be specified in the constructor. For multiplicative variants, \e omega
     * is a relaxation paramter (as in SOR), whereas for additive variants \e omega is a damping
     * parameter (as in damped Jacobi).
     *
     * <b>Guidelines, Remarks, Warnings & Tips</b>\n
     * - In many cases, setting the relaxation/damping parameter \e omega to a value less than 1
     *   can 'cure' a diverging Vanka variant. So before trying out a more expensive variant
     *   to achieve convergence, you may want to play around with \e omega first, e.g. by setting
     *   it to 0.8.
     * - In many cases and for any variant, choosing \e omega > 1 leads to divergence, so don't do that.
     * - If the pressure space is discontinuous, one should always choose block variants
     *   instead of nodal variants, as these require less memory, are less computationally
     *   expensive and lead to faster convergence in most cases.
     * - If the pressure space is continuous, one should stick to the nodal variants.\n
     *   Moreover, never use full variants for continuous pressure spaces, as these may require
     *   astronomical amounts of memory for the factorisation!
     * - Diagonal variants usually require (significantly) less memory than full variants.
     *   Therefore, when trying to configure a Vanka for a new type of problem and/or discretisation,
     *   one should first try whether a diagonal variant works before trying a full variant.
     * - For more complex problems (like convection-dominant Navier-Stokes), the full variants
     *   are often more stable than the diagonal variants.
     * - Multiplicative variants are roughly twice as computationally expensive as additive variants.
     * - Surprisingly, additive variants seem to converge \e faster than their multiplicative
     *   counterparts in many cases (at least for linear Stokes equations).
     *
     * \see \cite{Vanka86}
     *
     * \todo Find out how to use Vanka for Global::Matrix systems...
     *
     * \author Peter Zajac
     */
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename Filter_>
    class Vanka<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>, Filter_> :
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

      /// velocity vector type
      typedef typename MatrixD_::VectorTypeR VectorV;
      /// pressure vector type
      typedef typename MatrixB_::VectorTypeR VectorP;

    protected:
      // our Vanka matrix types
      typedef Intern::VankaMatrix<MatrixA_> VankaMatrixA;
      typedef Intern::VankaMatrix<MatrixB_> VankaMatrixB;
      typedef Intern::VankaMatrix<MatrixD_> VankaMatrixD;

      // dimension sanity checks
      static_assert(VankaMatrixA::row_dim == VankaMatrixA::col_dim, "Matrix A has invalid dimensions");
      static_assert(VankaMatrixA::row_dim == VankaMatrixB::row_dim, "Matrices A and B have incompatible dimensions");
      static_assert(VankaMatrixA::col_dim == VankaMatrixD::col_dim, "Matrices A and D have incompatible dimensions");
      static_assert(VankaMatrixB::col_dim == VankaMatrixD::row_dim, "Matrices B and D have incompatible dimensions");

      // for now, B and D cannot have more than 1 pressure dimension...
      static_assert(VankaMatrixB::col_dim == 1, "Invalid pressure space dimension");

      /// the system matrix
      const MatrixType& _matrix;
      /// the system filter
      const FilterType& _filter;
      /// the Vanka type
      VankaType _type;
      /// relaxation parameter
      DataType _omega;
      /// desired number of iterations
      Index _num_iter;
      /// maximum velocity block degree
      IndexType _degree_v;
      /// maximum pressure block degree
      IndexType _degree_p;
      /// velocity block structure
      std::vector<IndexType> _block_v_ptr, _block_v_idx;
      /// pressure block structure
      std::vector<IndexType> _block_p_ptr, _block_p_idx;
      /// factorisation data
      std::vector<DataType> _data;
      /// local defect/correction vectors
      std::vector<DataType> _vdef, _vcor;
      /// temporary vectors (additive types only)
      VectorType _vec_scale, _vec_tmp1, _vec_tmp2;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The saddle-point system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] type
       * The desired type of Vanka.
       *
       * \param[in] omega
       * The relaxation parameter.
       *
       * \param[in] num_iter
       * The number of iterations to be performed.
       */
      explicit Vanka(const MatrixType& matrix, const FilterType& filter, VankaType type,
        DataType omega = DataType(1), Index num_iter = Index(1)) :
        _matrix(matrix),
        _filter(filter),
        _type(type),
        _omega(omega),
        _num_iter(num_iter)
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Vanka";
      }

      /// Performs symbolic factorisation
      virtual void init_symbolic() override
      {
        bool block = ((int(_type) & 0x010) != 0);
        bool multi = ((int(_type) & 0x100) == 0);

        // compute pressure block graph
        if(block)
        {
          this->_build_p_block();
        }
        else
        {
          this->_build_p_nodal();
        }

        // compute velocity block graph
        this->_build_v_block();

        // allocate memory for numerical factorisation
        this->_alloc_data();

        // allocate temporary vector for additive
        if(!multi)
        {
          this->_vec_scale = this->_matrix.create_vector_r();
          this->_vec_tmp1 = this->_matrix.create_vector_r();
          this->_vec_tmp2 = this->_matrix.create_vector_r();
        }
      }

      /// Releases the symbolic factorisation data
      virtual void done_symbolic() override
      {
        _vec_tmp2.clear();
        _vec_tmp1.clear();
        _vec_scale.clear();
        _vdef.clear();
        _vcor.clear();
        _data.clear();
        _block_v_ptr.clear();
        _block_v_idx.clear();
        _block_p_ptr.clear();
        _block_p_idx.clear();
      }

      /// Performs numeric factorisation
      virtual void init_numeric() override
      {
        bool full = ((int(_type) & 0x001) != 0);
        bool multi = ((int(_type) & 0x100) == 0);

        if(full)
        {
          this->_factor_full();
        }
        else
        {
          this->_factor_diag();
        }

        if(!multi)
        {
          this->_calc_scale();
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        bool full = ((int(_type) & 0x001) != 0);
        if(full)
        {
          this->_apply_full(vec_cor, vec_def);
        }
        else
        {
          this->_apply_diag(vec_cor, vec_def);
        }

        return Status::success;
      }

      std::size_t data_size() const
      {
        return _data.size();
      }

    protected:
      /**
       * \brief Builds the 'blocked' pressure graph
       *
       * This function build the pressure graph by trying to re-create patched 'blocks'
       * by conjugating the D and B matrices of the system matrix.
       */
      void _build_p_block()
      {
        const auto& mat_b = _matrix.block_b().first();
        const auto& mat_d = _matrix.block_d().first();

        // fetch matrix dimensions
        const IndexType m = IndexType(mat_d.rows());

        // fetch the matrix arrays
        const IndexType* row_ptr_b = mat_b.row_ptr();
        const IndexType* col_idx_b = mat_b.col_ind();
        const IndexType* row_ptr_d = mat_d.row_ptr();
        const IndexType* col_idx_d = mat_d.col_ind();

        // clear block arrays
        _block_p_ptr.clear();
        _block_p_idx.clear();
        _block_p_ptr.push_back(IndexType(0));

        // local map
        std::map<IndexType, int> map_s;

        // allocate mask vector
        std::vector<int> mask(std::size_t(m), 0);

        // loop over all pressure nodes
        for(Index i(0); i < m; ++i)
        {
          // is this node already processed?
          if(mask[i] != 0)
            continue;

          // clear the local map
          map_s.clear();

          // okay, loop over all entries of D
          for(IndexType kd = row_ptr_d[i]; kd < row_ptr_d[i+1]; ++kd)
          {
            // fetch the column index of D
            const IndexType col_d = col_idx_d[kd];

            // loop over the row of B
            for(IndexType kb = row_ptr_b[col_d]; kb < row_ptr_b[col_d+1]; ++kb)
            {
              // fetch the column index of B
              const IndexType col_b = col_idx_b[kb];

              // insert into map
              auto ib = map_s.emplace(col_b, 1);
              if(!ib.second)
              {
                // node already exists, increment counter then
                ++(ib.first->second);
              }
            }
          }

          // compute the map degree
          int ideg = 0;
          for(auto it = map_s.begin(); it != map_s.end(); ++it)
            ideg = Math::max(ideg, it->second);

          // loop over all map entries with maximum degree
          for(auto it = map_s.begin(); it != map_s.end(); ++it)
          {
            if(ideg == it->second)
            {
              // insert into map
              _block_p_idx.push_back(it->first);
              mask[it->first] = 1;
            }
          }

          // push row
          _block_p_ptr.push_back(IndexType(_block_p_idx.size()));
        }
      }

      /**
       * \brief Builds the 'nodal' pressure graph
       */
      void _build_p_nodal()
      {
        const auto& mat_d = _matrix.block_d().first();
        const IndexType m = IndexType(mat_d.rows());

        // clear block arrays
        _block_p_ptr.clear();
        _block_p_idx.clear();
        _block_p_ptr.reserve(m+1);
        _block_p_idx.reserve(m);
        for(IndexType i(0); i < m; ++i)
        {
          _block_p_ptr.push_back(i);
          _block_p_idx.push_back(i);
        }
        _block_p_ptr.push_back(m);
      }

      /// Builds the velocity graph
      void _build_v_block()
      {
        // fetch the matrix arrays
        const auto& mat_d = _matrix.block_d().first();
        const IndexType* row_ptr_d = mat_d.row_ptr();
        const IndexType* col_idx_d = mat_d.col_ind();

        // fetch number of pressure blocks
        const IndexType m = IndexType(_block_p_ptr.size()-1);

        // clear block arrays
        _block_v_ptr.clear();
        _block_v_idx.clear();
        _block_v_ptr.reserve(m+1);
        _block_v_ptr.push_back(IndexType(0));

        std::set<IndexType> set_v;

        // loop over all pressure blocks
        for(IndexType i(0); i < m; ++i)
        {
          // loop over all pressure dofs in the current block
          for(IndexType j = _block_p_ptr[i]; j < _block_p_ptr[i+1]; ++j)
          {
            // get the pressure dof index
            const IndexType pix = _block_p_idx[j];

            // now loop over the corresponding row of D
            for(IndexType k = row_ptr_d[pix]; k < row_ptr_d[pix+1]; ++k)
            {
              // insert velocity dof into block set
              set_v.insert(col_idx_d[k]);
            }
          }

          // push velocity block
          for(auto it = set_v.begin(); it != set_v.end(); ++it)
          {
            _block_v_idx.push_back(*it);
          }

          // update pointer
          _block_v_ptr.push_back(IndexType(_block_v_idx.size()));

          // clear local set
          set_v.clear();
        }
      }

      /// Allocates the data arrays for numerical factorisation
      void _alloc_data()
      {
        // get the dimension of our system
        const IndexType dim = IndexType(Intern::VankaMatrix<MatrixA_>::row_dim);

        // use diagonal?
        bool diag = ((int(_type) & 1) == 0);

        // number of non-zero entries
        IndexType nze = IndexType(0);

        // reset degrees
        _degree_v = _degree_p = IndexType(0);

        // loop over all blocks
        const IndexType m = IndexType(_block_p_ptr.size()-1);
        for(IndexType i(0); i < m; ++i)
        {
          // get number of velocity and pressure dofs
          IndexType nv = _block_v_ptr[i+1] - _block_v_ptr[i];
          IndexType np = _block_p_ptr[i+1] - _block_p_ptr[i];

          // update degrees
          _degree_v = Math::max(nv, _degree_v);
          _degree_p = Math::max(np, _degree_p);

          // update count
          if(diag)
          {
            // main diagonal of A
            nze += dim * nv;
            // matrices B and D
            nze += IndexType(2) * dim * nv * np;
            // Schur-Complement matrix S
            nze += np * np;
          }
          else
          {
            // dense local systems
            nze += Math::sqr(dim*nv + np);
          }
        }

        // finally, allocate our memory
        _data.resize(nze);
        _vdef.resize(dim*_degree_v + _degree_p);
        _vcor.resize(dim*_degree_v + _degree_p);
      }

      /// Calculate the scaling vector for additive Vanka
      void _calc_scale()
      {
        // create a Vanka Vector
        this->_vec_scale.format();
        Intern::VankaVector<VectorV> vanka_v(this->_vec_scale.template at<0>());
        Intern::VankaVector<VectorP> vanka_p(this->_vec_scale.template at<1>());

        // get our data arrays
        const IndexType* vptr = _block_v_ptr.data();
        const IndexType* vidx = _block_v_idx.data();
        const IndexType* pptr = _block_p_ptr.data();
        const IndexType* pidx = _block_p_idx.data();

        // create a vector of ones
        std::vector<DataType> vec_one(vanka_v.dim*_degree_v + _degree_p, DataType(1));
        const DataType* vone = vec_one.data();

        // loop over all blocks
        const IndexType nblocks = IndexType(_block_v_ptr.size() - 1);
        for(IndexType iblock(0); iblock < nblocks; ++iblock)
        {
          // get number of velocity and pressure dofs
          const IndexType nv = vptr[iblock+1] - vptr[iblock];
          const IndexType np = pptr[iblock+1] - pptr[iblock];

          // scatter ones
          vanka_v.scatter_cor(DataType(1), vone, &vidx[vptr[iblock]], nv, IndexType(0));
          vanka_p.scatter_cor(DataType(1), vone, &pidx[pptr[iblock]], np, IndexType(0));
        }

        // invert components
        this->_vec_scale.component_invert(this->_vec_scale);
      }

      /// Performs the 'full' numerical factorisation
      void _factor_full()
      {
        Intern::VankaMatrix<MatrixA_> vanka_a(_matrix.block_a());
        Intern::VankaMatrix<MatrixB_> vanka_b(_matrix.block_b());
        Intern::VankaMatrix<MatrixD_> vanka_d(_matrix.block_d());

        // get the dimension of our system
        const IndexType dim = IndexType(vanka_a.row_dim);

        // get our data arrays
        DataType* data = _data.data();
        const IndexType* vptr = _block_v_ptr.data();
        const IndexType* vidx = _block_v_idx.data();
        const IndexType* pptr = _block_p_ptr.data();
        const IndexType* pidx = _block_p_idx.data();

        // format block data
        ::memset(data, 0, sizeof(DataType) * _data.size());

        // allocate pivot array
        std::vector<IndexType> pivot(3*(dim*_degree_v + _degree_p));

        // current block offset
        IndexType block_offset = IndexType(0);

        // loop over all blocks
        const IndexType nblocks = IndexType(_block_v_ptr.size() - 1);
        for(IndexType iblock(0); iblock < nblocks; ++iblock)
        {
          // get number of velocity and pressure dofs
          const IndexType nv = vptr[iblock+1] - vptr[iblock];
          const IndexType np = pptr[iblock+1] - pptr[iblock];

          // get our local indices
          const IndexType* loc_vidx = &vidx[vptr[iblock]];
          const IndexType* loc_pidx = &pidx[pptr[iblock]];

          // compute matrix stride
          const IndexType n = dim*nv + np;
          const IndexType block_size = n*n;

          // get our block data array pointer
          DataType* block_data = &data[block_offset];

          // gather matrix a
          std::pair<IndexType,IndexType> ao =
            vanka_a.gather_full(block_data, loc_vidx, loc_vidx, nv, nv, n);
          vanka_b.gather_full(block_data, loc_vidx, loc_pidx, nv, np, n, IndexType(0), ao.second);
          vanka_d.gather_full(block_data, loc_pidx, loc_vidx, np, nv, n, ao.first, IndexType(0));

          // invert local matrix block
          Math::invert_matrix(n, n, block_data, pivot.data());

          // make sure the matrix inversion did not fail
          for(IndexType i(0); i < block_size; ++i)
          {
            // make sure we have a finite value
            if(!Math::isfinite(block_data[i]))
            {
              throw VankaFactorError();
            }
          }

          // increment block data offset
          block_offset += block_size;
        }
      }

      /// Applies the 'full' Vanka iteration
      void _apply_full(VectorType& vec_cor, const VectorType& vec_def)
      {
        // format correction vector
        vec_cor.format();

        // do we use the multiplicative variant?
        const bool multi = ((int(_type) & 0x100) == 0);

        // additive?
        if(!multi)
        {
          this->_vec_tmp1.copy(vec_def);
          this->_vec_tmp2.format();
        }

        // get our sub-vectors
        VectorV& vec_cv = (multi ? vec_cor.template at<0>() : this->_vec_tmp2.template at<0>());
        VectorP& vec_cp = (multi ? vec_cor.template at<1>() : this->_vec_tmp2.template at<1>());
        const VectorV& vec_dv = (multi ? vec_def.template at<0>() : this->_vec_tmp1.template at<0>());
        const VectorP& vec_dp = (multi ? vec_def.template at<1>() : this->_vec_tmp1.template at<1>());

        // create our vanka matrix and vector objects
        Intern::VankaMatrix<MatrixA_> vanka_a(_matrix.block_a());
        Intern::VankaMatrix<MatrixB_> vanka_b(_matrix.block_b());
        Intern::VankaMatrix<MatrixD_> vanka_d(_matrix.block_d());
        Intern::VankaVector<VectorV> vanka_v(vec_cv);
        Intern::VankaVector<VectorP> vanka_p(vec_cp);

        // get velocity vector dimension
        const IndexType velo_dim = IndexType(vanka_v.dim);

        // get block count
        const IndexType num_blocks = IndexType(_block_v_ptr.size()-1);

        // get block data arrays
        const IndexType* vptr = _block_v_ptr.data();
        const IndexType* vidx = _block_v_idx.data();
        const IndexType* pptr = _block_p_ptr.data();
        const IndexType* pidx = _block_p_idx.data();

        // get local data and vectors
        const DataType* lmat_data = _data.data();
        DataType* lcor = _vcor.data();
        DataType* ldef = _vdef.data();

        // iterate
        for(IndexType iter(0); iter < _num_iter; ++iter)
        {
          // additive variant?
          if((!multi) && (iter > IndexType(0)))
          {
            // compute current defect
            this->_matrix.apply(_vec_tmp1, vec_cor, vec_def, -DataType(1));
            _vec_tmp2.format();
          }

          // reset block offset
          IndexType block_offset = IndexType(0);

          // loop over all blocks
          for(IndexType iblock(0); iblock < num_blocks; ++iblock)
          {
            // get local sizes
            const IndexType nv = vptr[iblock+1] - vptr[iblock];
            const IndexType np = pptr[iblock+1] - pptr[iblock];

            // get our local indices
            const IndexType* loc_vidx = &vidx[vptr[iblock]];
            const IndexType* loc_pidx = &pidx[pptr[iblock]];

            // compute local degree
            const IndexType n = velo_dim * nv + np;

            // get our local matrix data
            const DataType* lmat = &lmat_data[block_offset];

            // get local vectors
            DataType* lcor_v = lcor;
            DataType* lcor_p = &lcor[velo_dim * nv];
            DataType* ldef_v = ldef;
            DataType* ldef_p = &ldef[velo_dim * nv];

            // gather local defect
            vanka_v.gather_def(ldef_v, vec_dv, loc_vidx, nv, IndexType(0));
            vanka_p.gather_def(ldef_p, vec_dp, loc_pidx, np, IndexType(0));

            // multiplicative variants only:
            if(multi)
            {
              // subtract A*x
              vanka_a.mult_cor(ldef_v, -DataType(1), vec_cv, loc_vidx, nv, IndexType(0));
              vanka_b.mult_cor(ldef_v, -DataType(1), vec_cp, loc_vidx, nv, IndexType(0));
              vanka_d.mult_cor(ldef_p, -DataType(1), vec_cv, loc_pidx, np, IndexType(0));
            }

            // solve local system
            for(IndexType i(0); i < n; ++i)
            {
              lcor[i] = DataType(0);
              for(IndexType j(0); j < n; ++j)
              {
                lcor[i] += lmat[i*n + j] * ldef[j];
              }
            }

            // scatter result
            vanka_v.scatter_cor(_omega, lcor_v, loc_vidx, nv, IndexType(0));
            vanka_p.scatter_cor(_omega, lcor_p, loc_pidx, np, IndexType(0));

            // update block offset
            block_offset += n*n;
          }

          // additive variant?
          if(!multi)
          {
            // multiply by scaling vector
            this->_vec_tmp2.component_product(this->_vec_tmp2, this->_vec_scale);

            // update correction vector
            vec_cor.axpy(vec_cor, this->_vec_tmp2);
          }

          // apply filter
          _filter.filter_cor(vec_cor);
        }
      }

      /// Performs the 'diagonal' numerical factorisation
      void _factor_diag()
      {
        Intern::VankaMatrix<MatrixA_> vanka_a(_matrix.block_a());
        Intern::VankaMatrix<MatrixB_> vanka_b(_matrix.block_b());
        Intern::VankaMatrix<MatrixD_> vanka_d(_matrix.block_d());

        // get the dimension of our system
        const IndexType dim = IndexType(vanka_a.row_dim);

        // get our data arrays
        DataType* data = _data.data();
        const IndexType* vptr = _block_v_ptr.data();
        const IndexType* vidx = _block_v_idx.data();
        const IndexType* pptr = _block_p_ptr.data();
        const IndexType* pidx = _block_p_idx.data();

        // format block data
        ::memset(data, 0, sizeof(DataType) * _data.size());

        // allocate pivot array (pressure only)
        std::vector<IndexType> pivot(3*_degree_p);

        // current block offset
        IndexType block_offset = IndexType(0);

        // loop over all blocks
        const IndexType nblocks = IndexType(_block_v_ptr.size() - 1);
        for(IndexType iblock(0); iblock < nblocks; ++iblock)
        {
          // get number of velocity and pressure dofs
          const IndexType nv = vptr[iblock+1] - vptr[iblock];
          const IndexType np = pptr[iblock+1] - pptr[iblock];
          const IndexType dnv = dim*nv;

          // get our local indices
          const IndexType* loc_vidx = &vidx[vptr[iblock]];
          const IndexType* loc_pidx = &pidx[pptr[iblock]];

          // compute local block size
          const IndexType block_size = dnv + IndexType(2) * dnv * np + np*np;

          // get our block data array pointer
          DataType* block_data = &data[block_offset];

          // get our local matrices
          DataType* loc_a =  block_data;
          DataType* loc_b = &block_data[dnv];
          DataType* loc_d = &block_data[dnv + dnv*np];
          DataType* loc_s = &block_data[dnv + IndexType(2)*dnv*np];

          // gather diag(A)
          vanka_a.gather_diag(loc_a, loc_vidx, nv);
          // gather B and D
          vanka_b.gather_full(loc_b, loc_vidx, loc_pidx, nv, np, np);
          vanka_d.gather_full(loc_d, loc_pidx, loc_vidx, np, nv, dnv);

          // invert diag(A) and pre-multiply D by diag(A)^{-1}
          for(IndexType i(0); i < dnv; ++i)
          {
            // invert a_ii
            loc_a[i] = DataType(1) / loc_a[i];

            // make sure we have a finite value
            if(!Math::isfinite(loc_a[i]))
            {
              throw VankaFactorError();
            }

            // pre-multiply D by diag(A)^{-1}
            for(IndexType j(0); j < np; ++j)
            {
              loc_d[j*dnv + i] *= loc_a[i];
            }
          }

          // calculate Schur-complement of A:
          // S := -D * diag(A)^{-1} * B
          for(IndexType i(0); i < np; ++i)
          {
            for(IndexType j(0); j < np; ++j)
            {
              DataType s = DataType(0);
              for(IndexType k(0); k < dnv; ++k)
              {
                // Note: D is already pre-multiplied by diag(A)^{-1}
                s += loc_d[i*dnv + k] * loc_b[k*np + j];
              }
              loc_s[i*np + j] = -s;
            }
          }

          // invert local matrix S
          if(np == IndexType(1))
          {
            loc_s[0] = DataType(1) / loc_s[0];
          }
          else
          {
            Math::invert_matrix(np, np, loc_s, pivot.data());
          }

          // ensure that S^{-1} is not bogus
          for(IndexType i(0); i < np*np; ++i)
          {
            // make sure we have a finite value
            if(!Math::isfinite(loc_s[i]))
            {
              throw VankaFactorError();
            }
          }

          // increment block data offset
          block_offset += block_size;
        }
      }

      /// Applies the 'diagonal' Vanka iteration
      void _apply_diag(VectorType& vec_cor, const VectorType& vec_def)
      {
        // format correction vector
        vec_cor.format();

        // do we use the multiplicative variant?
        const bool multi = ((int(_type) & 0x100) == 0);

        // additive?
        if(!multi)
        {
          this->_vec_tmp1.copy(vec_def);
          this->_vec_tmp2.format();
        }

        // get our sub-vectors
        VectorV& vec_cv = (multi ? vec_cor.template at<0>() : this->_vec_tmp2.template at<0>());
        VectorP& vec_cp = (multi ? vec_cor.template at<1>() : this->_vec_tmp2.template at<1>());
        const VectorV& vec_dv = (multi ? vec_def.template at<0>() : this->_vec_tmp1.template at<0>());
        const VectorP& vec_dp = (multi ? vec_def.template at<1>() : this->_vec_tmp1.template at<1>());

        // create our vanka matrix and vector objects
        Intern::VankaMatrix<MatrixA_> vanka_a(_matrix.block_a());
        Intern::VankaMatrix<MatrixB_> vanka_b(_matrix.block_b());
        Intern::VankaMatrix<MatrixD_> vanka_d(_matrix.block_d());
        Intern::VankaVector<VectorV> vanka_v(vec_cv);
        Intern::VankaVector<VectorP> vanka_p(vec_cp);

        // get velocity vector dimension
        const IndexType velo_dim = IndexType(vanka_v.dim);

        // get block count
        const IndexType num_blocks = IndexType(_block_v_ptr.size()-1);

        // get block data arrays
        const IndexType* vptr = _block_v_ptr.data();
        const IndexType* vidx = _block_v_idx.data();
        const IndexType* pptr = _block_p_ptr.data();
        const IndexType* pidx = _block_p_idx.data();

        // get local data and vectors
        const DataType* data = _data.data();
        DataType* lcor = _vcor.data();
        DataType* ldef = _vdef.data();

        // iterate
        for(IndexType iter(0); iter < _num_iter; ++iter)
        {
          // additive variant?
          if((!multi) && (iter > IndexType(0)))
          {
            // compute current defect
            this->_matrix.apply(_vec_tmp1, vec_cor, vec_def, -DataType(1));
            _vec_tmp2.format();
          }

          // reset block offset
          IndexType block_offset = IndexType(0);

          // loop over all blocks
          for(IndexType iblock(0); iblock < num_blocks; ++iblock)
          {
            // get local sizes
            const IndexType nv = vptr[iblock+1] - vptr[iblock];
            const IndexType np = pptr[iblock+1] - pptr[iblock];
            const IndexType dnv = velo_dim*nv;

            // get our local indices
            const IndexType* loc_vidx = &vidx[vptr[iblock]];
            const IndexType* loc_pidx = &pidx[pptr[iblock]];

            // compute local block size
            const IndexType block_size = dnv + IndexType(2) * dnv * np + np*np;

            // get our block data array pointer
            const DataType* block_data = &data[block_offset];

            // get our local matrices
            const DataType* loc_a =  block_data;
            const DataType* loc_b = &block_data[dnv];
            const DataType* loc_d = &block_data[dnv + dnv*np];
            const DataType* loc_s = &block_data[dnv + IndexType(2)*dnv*np];

            // get local vectors
            DataType* lcor_v = lcor;
            DataType* lcor_p = &lcor[dnv];
            DataType* ldef_v = ldef;
            DataType* ldef_p = &ldef[dnv];

            // gather local defect
            vanka_v.gather_def(ldef_v, vec_dv, loc_vidx, nv, IndexType(0));
            vanka_p.gather_def(ldef_p, vec_dp, loc_pidx, np, IndexType(0));

            // multiplicative variants only:
            if(multi)
            {
              // subtract A*x
              vanka_a.mult_cor(ldef_v, -DataType(1), vec_cv, loc_vidx, nv, IndexType(0));
              vanka_b.mult_cor(ldef_v, -DataType(1), vec_cp, loc_vidx, nv, IndexType(0));
              vanka_d.mult_cor(ldef_p, -DataType(1), vec_cv, loc_pidx, np, IndexType(0));
            }

            // update pressure RHS:
            // g_p := f_p - D * diag(A)^{-1} * f_u
            for(IndexType i(0); i < np; ++i)
            {
              DataType r = DataType(0);
              for(IndexType j(0); j < dnv; ++j)
              {
                // Note: D is already pre-multiplied by diag(A)^{-1}
                r += loc_d[i*dnv + j] * ldef_v[j];
              }
              ldef_p[i] -= r;
            }

            // solve pressure:
            // p := S^{-1} * g_p
            for(IndexType i(0); i < np; ++i)
            {
              DataType r = DataType(0);
              for(IndexType j(0); j < np; ++j)
              {
                r += loc_s[i*np + j] * ldef_p[j];
              }
              lcor_p[i] = r;
            }

            // update velocity RHS and solve velocity
            for(IndexType i(0); i < dnv; ++i)
            {
              DataType xb = DataType(0);
              for(IndexType j(0); j < np; ++j)
              {
                xb += loc_b[i*np + j] * lcor_p[j];
              }
              // solve: u := diag(A)^{-1} * (f_u - B*p)
              lcor_v[i] = loc_a[i] * (ldef_v[i] - xb);
            }

            // scatter result
            vanka_v.scatter_cor(_omega, lcor_v, loc_vidx, nv, IndexType(0));
            vanka_p.scatter_cor(_omega, lcor_p, loc_pidx, np, IndexType(0));

            // update block offset
            block_offset += block_size;
          }

          // additive variant?
          if(!multi)
          {
            // multiply by scaling vector
            this->_vec_tmp2.component_product(this->_vec_tmp2, this->_vec_scale);

            // update correction vector
            vec_cor.axpy(vec_cor, this->_vec_tmp2);
          }

          // apply filter
          _filter.filter_cor(vec_cor);
        }
      }
    }; // class Vanka<...>

    /**
     * \brief Creates a new Vanka solver object
     *
     * \param[in] matrix
     * The saddle-point system matrix.
     *
     * \param[in] filter
     * The system filter
     *
     * \param[in] type
     * Specifies the Vanka type.
     *
     * \param[in] omega
     * The relaxation parameter.
     *
     * \param[in] num_iter
     * The number of Vanka iterations to be performed.
     *
     * \returns
     * A shared pointer to a new Vanka object.
     */
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename Filter_>
    std::shared_ptr<Vanka<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>, Filter_>> new_vanka(
      const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& matrix,
      const Filter_& filter,
      VankaType type,
      typename MatrixA_::DataType omega = typename MatrixA_::DataType(1),
      Index num_iter = Index(1))
    {
      return std::make_shared<Vanka<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>, Filter_>>
        (matrix, filter, type, omega, num_iter);
    }
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_VANKA_HPP
