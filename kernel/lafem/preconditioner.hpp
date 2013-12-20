#pragma once
#ifndef KERNEL_LAFEM_PRECONDITIONER_HPP
#define KERNEL_LAFEM_PRECONDITIONER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/algorithm.hpp>
#include <kernel/lafem/scale.hpp>
#include <kernel/lafem/algorithm.hpp>
#include <kernel/lafem/axpy.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_, typename MT_, typename VT_>
    class Preconditioner
    {
    public:
      virtual void apply(VT_ & out, const VT_ & in) = 0;
    };



    /**
     * \brief No preconditioner.
     *
     * This class represents a dummy for a preconditioner.
     *
     * \author Dirk Ribbrock
     */
    template <typename Algo_, typename MT_, typename VT_>
    class NonePreconditioner : public virtual Preconditioner<Algo_, MT_, VT_>
    {
    public:
      NonePreconditioner()
      {
      }

      static String type_name()
      {
        return "None_Preconditioner";
      }

      virtual void apply(VT_ & out, const VT_ & in)
      {
        copy(out, in);
      }
    };



    /**
     * \brief Jacobi-Preconditioner.
     *
     * This class represents the Jacobi-Preconditioner.
     *
     * \author Dirk Ribbrock
     */
    template <typename Algo_, typename MT_, typename VT_>
    class JacobiPreconditioner : public virtual Preconditioner<Algo_, MT_, VT_>
    {
    private:
      VT_ _jac;

    public:
      JacobiPreconditioner(const MT_ & A, typename VT_::DataType damping) :
        _jac(A.rows())
      {
        if (A.columns() != A.rows())
          throw InternalError("Matrix is not square!");

        for (Index i(0) ; i < _jac.size() ; ++i)
          _jac(i, damping / A(i, i));
      }

      static String type_name()
      {
        return "Jacobi_Preconditioner";
      }

      virtual void apply(VT_ & out, const VT_ & in)
      {
        out.template component_product<Algo_>(_jac, in);
      }
    };



    /**
     * \brief Gauss-Seidel-Preconditioner.
     *
     * This class represents a dummy class for the Gauss-Seidel-Preconditioner.
     *
     * \author Christoph Lohmann
     */
    template <typename Algo_, typename MT_, typename VT_>
    class GaussSeidelPreconditioner : public virtual Preconditioner<Algo_, MT_, VT_>
    {
    };


    /**
     * \brief Gauss-Seidel-Preconditioner for CSR-matrices.
     *
     * This class specializes the Gauss-Seidel-Preconditioner for CSR-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Arch_, typename DT_>
    class GaussSeidelPreconditioner<Algo::Generic, SparseMatrixCSR<Arch_, DT_>,
                                    DenseVector<Arch_, DT_> >
      : public virtual Preconditioner<Algo::Generic, SparseMatrixCSR<Arch_, DT_>,
                                      DenseVector<Arch_, DT_> >
    {
    private:
      DT_ damping;
      const SparseMatrixCSR<Arch_, DT_> & A;

    public:
      GaussSeidelPreconditioner(const SparseMatrixCSR<Arch_, DT_> & A,
                                DT_ damping) :
        damping(damping),
        A(A)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError("Matrix is not square!");
        }
      }

      static String type_name()
      {
        return "GaussSeidel_Preconditioner";
      }

      virtual void apply(DenseVector<Arch_, DT_> & out,
                         const DenseVector<Arch_, DT_> & in)
      {

        // copy in-vector to out-vector
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pA_val(A.val());
        const Index * pA_col_ind(A.col_ind());
        const Index * pA_row(A.row_ptr());


        Index row;

        // iteration over all columns
        for (Index i(0); i < A.columns(); ++i)
        {

          // iteration over all elements on the left side of the main-diagonal
          row = pA_row[i];

          while(pA_col_ind[row] < i)
          {
            pout[i] -= pA_val[row] * pout[pA_col_ind[row]];
            ++row;
          }

          // divide by the element of the main diagonal
          pout[i] /= pA_val[row];
        }

        // damping of solution
        out.template scale<Algo_>(out, damping);

      }
    };


    /**
     * \brief Gauss-Seidel-Preconditioner for COO-matrices.
     *
     * This class specializes the Gauss-Seidel-Preconditioner for COO-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Arch_, typename DT_>
    class GaussSeidelPreconditioner<Algo::Generic, SparseMatrixCOO<Arch_, DT_>,
                                    DenseVector<Arch_, DT_> >
      : public virtual Preconditioner<Algo::Generic, SparseMatrixCOO<Arch_, DT_>,
                                      DenseVector<Arch_, DT_> >
    {
    private:
      DT_ damping;
      const SparseMatrixCOO<Arch_, DT_> & A;

    public:
      GaussSeidelPreconditioner(const SparseMatrixCOO<Arch_, DT_> & A,
                                DT_ damping) :
        damping(damping),
        A(A)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError("Matrix is not square!");
        }
      }

      static String type_name()
      {
        return "GaussSeidel_Preconditioner";
      }

      virtual void apply(DenseVector<Arch_, DT_> & out,
                         const DenseVector<Arch_, DT_> & in)
      {

        // copy in-vector to out-vector
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pA_val(A.val());
        const Index * pA_column(A.column());
        const Index * pA_row(A.row());


        Index row = 0;

        // Notice, that the elements of A are sorted by the column.
        // So we can solve the system with only one iteration over the data-array

        // iteration over all columns
        for (Index i(0); i < A.columns(); ++i)
        {

          // go to the i-th column of A
          while (pA_row[row] < i)
          {
            ++row;
          }

          // iteration over all elements on the left side of the main diagonal
          while (pA_column[row] < i)
          {
            pout[i] -= pA_val[row] * out(pA_column[row]);
            ++row;
          }

          // divide by the element of the main diagonal
          pout[i] /= pA_val[row];
        }


        // damping of solution
        out.template scale<Algo_>(out, damping);

      }
    };


    /**
     * \brief Gauss-Seidel-Preconditioner for ELL-matrices.
     *
     * This class specializes the Gauss-Seidel-Preconditioner for ELL-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Arch_, typename DT_>
    class GaussSeidelPreconditioner<Algo::Generic, SparseMatrixELL<Arch_, DT_>,
                                    DenseVector<Arch_, DT_> >
      : public virtual Preconditioner<Algo::Generic, SparseMatrixELL<Arch_, DT_>,
                                      DenseVector<Arch_, DT_> >
    {
    private:
      DT_ damping;
      const SparseMatrixELL<Arch_, DT_> & A;

    public:
      GaussSeidelPreconditioner(const SparseMatrixELL<Arch_, DT_> & A,
                                DT_ damping) :
        damping(damping),
        A(A)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError("Matrix is not square!");
        }
      }

      static String type_name()
      {
        return "GaussSeidel_Preconditioner";
      }

      virtual void apply(DenseVector<Arch_, DT_> & out,
                         const DenseVector<Arch_, DT_> & in)
      {

        // copy in-vector to out-vector
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pA_val(A.Ax());
        const Index * pAj(A.Aj());


        Index col;
        const Index stride(A.stride());

        // iteration over all rows
        for (Index i(0); i < A.rows(); ++i)
        {
          // search entry on the main diagonal
          col = i;

          while(pAj[col] < i)
          {
            col += stride;
          }

          // divide by the element of the main diagonal
          pout[i] /= pA_val[col];

          // iteration over all cells under the main-diagonal
          for (Index j(i+1); j < A.rows(); ++j)
          {
            col = j;
            while(pAj[col] < i)
            {
              col += stride;
            }

            if (pAj[col] == i)
            {
               pout[j] -= pA_val[col] * pout[i];
            }
          }
        }

        // damping of solution
        out.template scale<Algo_>(out, damping);

      }
    };



    /**
     * \brief Polynomial-Preconditioner.
     *
     * This class represents the Neumann-Polynomial-Preconditioner
     *
     * \author Christoph Lohmann
     */
    template <typename Algo_, typename MT_, typename VT_>
    class PolynomialPreconditioner : public virtual Preconditioner<Algo_, MT_, VT_>
    {
    private:
      const MT_ & A;        // system-matrix
      Index m;              // order m of preconditioner
      VT_ aux;              // auxilary-vector
      VT_ * pscale;         // scaling-vector (if check is true)

    public:

    /**
     * \brief Constuctor of Neumann-Polynomial-Preconditioner.
     *
     * param[in] A system-matrix
     * param[in] m order of the polynom
     * param[in] bscale optional paramter: true convergence of neumann-series
     *                                     for strict diagonaldominant matrices
     *                                     false (default), dont scale at all
     *
     * \author Christoph Lohmann
     */
      PolynomialPreconditioner(const MT_ & A, Index m, bool bscale = false) :
        A(A),
        m(m),
        aux(A.rows()),
        pscale(nullptr)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError("Matrix is not square!");
        }

        typedef typename VT_::DataType DT_;

        if (bscale == true)
        {
          pscale = new VT_(A.rows());
          DT_ *pscale_elem(pscale->elements());
          for (Index i = 0; i < A.rows(); ++i)
          {
            pscale_elem[i] = -DT_(1.0) / A(i,i);
          }
        }
      }

      ~PolynomialPreconditioner()
      {
        if (pscale != nullptr)
        {
          delete pscale;
        }
      }

      static String type_name()
      {
        return "Polynomial_Preconditioner";
      }

      virtual void apply(VT_ & out, const VT_ & in)
      {
        typedef typename VT_::DataType DT_;

        /*
        * preconditioner is given by
        *   M^-1 = I + (I - A) + ... + (I - A)^m
        *
        * the preconditioner only works, if
        *   ||I - A||_2 < 1.
        */

        VT_ * pauxs[2];

        if (m%2 == 0)
        {
          pauxs[0] = &out;
          pauxs[1] = &aux;
        }
        else
        {
          pauxs[1] = &out;
          pauxs[0] = &aux;

        }

        if (pscale == nullptr)
        {
          copy(*pauxs[0], in);

          for (Index i = 1; i <= m; ++i)
          {
            *pauxs[i%2].template axpy<Algo_>(DT_(1.0), in, *pauxs[(i-1)%2]);
            *pauxs[i%2].template axpy<Algo_>(DT_(-1.0), A, *pauxs[(i-1)%2], *pauxs[i%2]);
          }
        }
        else
        {
          *pauxs[0].template component_product<Algo_>(*pscale, in);
          out.template scale<Algo_>(out, DT_(-1.0));

          for (Index i = 1; i <= m; ++i)
          {
            *pauxs[i%2].template axpy<Algo_>(DT_(1.0), in, *pauxs[(i-1)%2]);
            *pauxs[i%2].template axpy<Algo_>(*pscale, A, *pauxs[(i-1)%2], *pauxs[i%2]);
          }
        }


      }
    };

  } // namespace LAFEM
} // namespace FEAST




#endif // KERNEL_LAFEM_PRECONDITIONER_HPP
