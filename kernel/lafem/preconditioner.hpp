#pragma once
#ifndef KERNEL_LAFEM_PRECONDITIONER_HPP
#define KERNEL_LAFEM_PRECONDITIONER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <vector>

namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_, typename MT_, typename VT_>
    class Preconditioner
    {
    public:
      virtual ~Preconditioner()
      {
      }

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
    class NonePreconditioner : public Preconditioner<Algo_, MT_, VT_>
    {
    public:
      virtual ~NonePreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * Creates a dummy preconditioner
       */
      NonePreconditioner()
      {
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "None_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector, which is applied to the preconditioning
       */
      virtual void apply(VT_ & out, const VT_ & in)
      {
        out.copy(in);
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
    class JacobiPreconditioner : public Preconditioner<Algo_, MT_, VT_>
    {
    private:
      VT_ _jac;

    public:
      virtual ~JacobiPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] damping A damping-parameter
       *
       * Creates a Jacobi preconditioner to the given matrix and damping-parameter
       */
      JacobiPreconditioner(const MT_ & A, const typename VT_::DataType damping) :
        _jac(A.rows())
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        const Index n(A.rows());

        for (Index i(0) ; i < n ; ++i)
        {
          _jac(i, damping / A(i, i));
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "Jacobi_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
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
    class GaussSeidelPreconditioner : public Preconditioner<Algo_, MT_, VT_>
    {
    };


    /**
     * \brief Gauss-Seidel-Preconditioner for CSR-matrices.
     *
     * This class specializes the Gauss-Seidel-Preconditioner for CSR-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class GaussSeidelPreconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                                    DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const DT_ damping;
      const SparseMatrixCSR<Mem_, DT_> & A;

    public:
      virtual ~GaussSeidelPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] damping A damping-parameter
       *
       * Creates a Gauss-Seidel preconditioner to the given matrix and damping-parameter
       */
      GaussSeidelPreconditioner(const SparseMatrixCSR<Mem_, DT_> & A,
                                DT_ damping) :
        damping(damping),
        A(A)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "GaussSeidel_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        out.copy(in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(A.val());
        const Index * pcol_ind(A.col_ind());
        const Index * prow(A.row_ptr());
        const Index n(A.rows());

        Index col;

        // iteration over all rows
        for (Index i(n); i > 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          col = prow[i+1]-1;

          while(pcol_ind[col] > i)
          {
            pout[i] -= pval[col] * pout[pcol_ind[col]];
            --col;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        // damping of solution
        out.template scale<Algo::Generic>(out, damping);
      }
    };


    /**
     * \brief Gauss-Seidel-Preconditioner for COO-matrices.
     *
     * This class specializes the Gauss-Seidel-Preconditioner for COO-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class GaussSeidelPreconditioner<Algo::Generic, SparseMatrixCOO<Mem_, DT_>,
                                    DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCOO<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const DT_ damping;
      const SparseMatrixCOO<Mem_, DT_> & A;

    public:
      virtual ~GaussSeidelPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] damping A damping-parameter
       *
       * Creates a Gauss-Seidel preconditioner to the given matrix and damping-parameter
       */
      GaussSeidelPreconditioner(const SparseMatrixCOO<Mem_, DT_> & A,
                                DT_ damping) :
        damping(damping),
        A(A)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "GaussSeidel_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        out.copy(in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(A.val());
        const Index * pcol(A.column());
        const Index * prow(A.row());
        const Index n(A.rows());

        Index col(A.used_elements() - 1);

        // iteration over all rows
        for (Index i(n); i > 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          while (prow[col] > i)
          {
            --col;
          }

          while(pcol[col] > i)
          {
            pout[i] -= pval[col] * pout[pcol[col]];
            --col;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        // damping of solution
        out.template scale<Algo::Generic>(out, damping);
      }
    };


    /**
     * \brief Gauss-Seidel-Preconditioner for ELL-matrices.
     *
     * This class specializes the Gauss-Seidel-Preconditioner for ELL-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class GaussSeidelPreconditioner<Algo::Generic, SparseMatrixELL<Mem_, DT_>,
                                    DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixELL<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const DT_ damping;
      const SparseMatrixELL<Mem_, DT_> & A;

    public:
      virtual ~GaussSeidelPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] damping A damping-parameter
       *
       * Creates a Gauss-Seidel preconditioner to the given matrix and damping-parameter
       */
      GaussSeidelPreconditioner(const SparseMatrixELL<Mem_, DT_> & A,
                                DT_ damping) :
        damping(damping),
        A(A)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "GaussSeidel_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        out.copy(in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(A.Ax());
        const Index * paj(A.Aj());
        const Index stride(A.stride());
        const Index * parl(A.Arl());
        const Index n(A.rows());


        Index col;

        // iteration over all rows
        for (Index i(n); i > 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          col = i + stride * (parl[i] - 1);

          while (paj[col] > i)
          {
            pout[i] -= pval[col] * pout[paj[col]];
            col -= stride;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        // damping of solution
        out.template scale<Algo::Generic>(out, damping);
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
    class PolynomialPreconditioner : public Preconditioner<Algo_, MT_, VT_>
    {
    private:
      const MT_ & A;        // system-matrix
      Index m;              // order m of preconditioner
      VT_ aux, tmp;         // auxilary-vector
      VT_ * pscale;         // scaling-vector (if check is true)

    public:
      /**
       * \brief Destructor
       *
       * deletes a vector, which was only used if bscale was true
       */
      ~PolynomialPreconditioner()
      {
        if (pscale != nullptr)
        {
          delete pscale;
        }
      }

      /**
       * \brief Constuctor of Neumann-Polynomial-Preconditioner.
       *
       * param[in] A system-matrix
       * param[in] m order of the polynom
       * param[in] bscale optional paramter: true convergence of neumann-series
       *                                     for strict diagonaldominant matrices
       *                                     false (default), dont scale at all
       *
       * Creates a Polynomial preconditioner to the given matrix and the given order
       */
      PolynomialPreconditioner(const MT_ & A, Index m, bool bscale = false) :
        A(A),
        m(m),
        aux(A.rows()),
        tmp(A.rows()),
        pscale(nullptr)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
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

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "Polynomial_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
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
          pauxs[0]->copy(in);

          for (Index i = 1; i <= m; ++i)
          {
            pauxs[i%2]->template axpy<Algo::Generic>(in, *pauxs[(i-1)%2]);
            //pauxs[i%2]->template axpy<Algo::Generic>(DT_(-1.0), A, *pauxs[(i-1)%2], *pauxs[i%2]);
            A.template apply<Algo::Generic>(*pauxs[i%2], *pauxs[(i-1)%2], *pauxs[i%2], -DT_(1));
          }
        }
        else
        {
          pauxs[0]->template component_product<Algo::Generic>(*pscale, in);
          out.template scale<Algo::Generic>(out, DT_(-1.0));

          for (Index i = 1; i <= m; ++i)
          {
            pauxs[i%2]->template axpy<Algo::Generic>(in, *pauxs[(i-1)%2]);
            //pauxs[i%2]->template axpy<Algo::Generic>(*pscale, A, *pauxs[(i-1)%2], *pauxs[i%2]);
            A.template apply<Algo::Generic>(tmp, *pauxs[(i-1)%2]);
            pauxs[i%2]->template component_product<Algo::Generic>(tmp, *pscale, *pauxs[i%2]);
          }
        } // function apply

      }
    };


    /**
     * \brief ILU(p)-Preconditioner.
     *
     * This class represents a dummy class for the ILU(p)-Preconditioner.
     *
     * \author Christoph Lohmann
     */
    template <typename Algo_, typename MT_, typename VT_>
    class ILUPreconditioner : public Preconditioner<Algo_, MT_, VT_>
    {
    };


    /**
     * \brief ILU(p)-Preconditioner for CSR-matrices.
     *
     * This class specializes the ILU(p)-Preconditioner for CSR-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class ILUPreconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                            DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const SparseMatrixCSR<Mem_, DT_> & A;
      SparseMatrixCSR<Mem_, DT_> LU;

    public:
      virtual ~ILUPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] p level of fillin
       *           if p = 0, the layout of A is used for the ILU-decomposition
       *
       * Creates a ILU preconditioner to the given matrix and level of fillin
       */
      ILUPreconditioner(const SparseMatrixCSR<Mem_, DT_> & A, const int p) :
        A(A)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        if (p == 0)
        {
          SparseMatrixCSR<Mem_, DT_> tLU(A.layout());
          LU = tLU;

          copy_entries(false);
        }
        else
        {
          symbolic_lu_factorisation(p);
          copy_entries();
        }

        create_lu();
      }

      /**
       * \brief Constructor
       *
       * param[in] LU external LU-matrix
       *           the LU-decomposition is not calculated internally;
       *           so this preconditioner only solves \f$y \leftarrow L^{-1} U^{-1} x\f$
       *
       * Creates a ILU preconditioner to the given LU-decomposition
       */
      ILUPreconditioner(const SparseMatrixCSR<Mem_, DT_> & LU) :
        A(LU),
        LU(LU)
      {
        if (LU.columns() != LU.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] layout An external layout for the LU-decomposition
       *
       * Creates a ILU preconditioner to the given matrix and layout
       */
      ILUPreconditioner(const SparseMatrixCSR<Mem_, DT_> & A,
                        const SparseLayout<Mem_, SparseMatrixCSR<Mem_, DT_>::LayoutType> & layout) :
        A(A),
        LU(layout)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        if (LU.columns() != LU.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        if (A.columns() != LU.columns())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrices have different sizes!");
        }

        copy_entries();
        create_lu();
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "ILU_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        out.copy(in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(LU.val());
        const Index * pcol_ind(LU.col_ind());
        const Index * prow(LU.row_ptr());
        const Index n(LU.rows());

        Index col;

        // __forward-insertion__
        // iteration over all rows
        for (Index i(0); i < n; ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          col = prow[i];

          while(pcol_ind[col] < i)
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
          col = prow[i+1]-1;

          while(pcol_ind[col] > i)
          {
            pout[i] -= pval[col] * pout[pcol_ind[col]];
            --col;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }
      } // function apply

    private:
      void create_lu()
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/1014
         *    section 5.5.3; Algo 5.13.; page 132
         */

        DT_ * plu(LU.val());
        const Index * pcol(LU.col_ind());
        const Index * prow(LU.row_ptr());
        const Index n(LU.rows());

        // integer work array of length n
        //   for saving the position of the diagonal-entries
        Index * pw = new Index[n];

        Index row_start;
        Index row_end;

        // iteration over all columns
        for (Index i(0); i < n; ++i)
        {
          row_start = prow[i];
          row_end = prow[i + 1];

          // iteration over all elements on the left side of the main-diagonal
          // k -> \tilde a_{ik}
          Index k = row_start;
          while(pcol[k] < i)
          {
            plu[k] /= plu[pw[pcol[k]]];
            Index m(pw[pcol[k]] + 1);
            // m -> \tilde a_{kj}
            for (Index j(k+1); j < row_end; ++j)
            {
              // j -> \tilde a_{ij}
              while (m < prow[pcol[k] + 1])
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
      } // function create_lu


      void symbolic_lu_factorisation(int p)
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/1014
         *    section 5.5.6; Algo 5.19.; page 142
         */

        const Index n(A.rows());
        const Index * pacol(A.col_ind());
        const Index * parow(A.row_ptr());

        // type of list-entries
        typedef std::pair<int, Index> PAIR_;
        // lists for saving the non-zero entries of L per row
        std::list<PAIR_> * ll = new std::list<PAIR_>[n];

        // vector for saving the iterators to the diag-entries
        std::list<PAIR_>::iterator * pldiag = new std::list<PAIR_>::iterator[n];
        // auxillary-iterators
        std::list<PAIR_>::iterator it, it1, it2;

        Index col_begin, col_end, col, col2;
        int l, l2, neues_level;

        // fill list with non-zero entries of A
        // and save iterators to the diag- and last entry of each row
        for (Index row(0); row < n; ++row)
        {
          std::list<PAIR_> & ll_row (ll[row]);
          col_begin = parow[row];
          col_end = parow[row + 1] - 1;

          for (Index k(col_begin); k <= col_end; ++k)
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
        // calculate number of non-zero-elements
        Index nnz(0);
        for (Index i(0); i < n; ++i)
        {
          nnz += Index(ll[i].size());
        }

        DenseVector<Mem_, DT_> val(nnz);
        DenseVector<Mem_, Index> col_ind(nnz);
        DenseVector<Mem_, Index> row_ptr(n+1);
        Index * pcol_ind(col_ind.elements());
        Index * prow_ptr(row_ptr.elements());

        Index k1(0);
        row_ptr(0,0);

        for (Index i(0); i < n; ++i)
        {
          for (it = ll[i].begin(); it != ll[i].end(); ++it)
          {
            pcol_ind[k1] = it->second;
            ++k1;
          }
          prow_ptr[i+1] = k1;
        }

        SparseMatrixCSR<Mem_, DT_> tLU(n, n, col_ind, val, row_ptr);
        LU = tLU;

        delete[] ll;
        delete[] pldiag;
      } // symbolic_lu_factorisation

      void copy_entries(bool check = true)
      {
        if (check == false)
        {
          DT_ * plu(LU.val());
          const DT_ * pa(A.val());
          const Index used_elements(LU.used_elements());

          // initialize LU array to A
          for (Index i(0); i < used_elements; ++i)
          {
            plu[i] = pa[i];
          }
        }
        else
        {
          DT_ * plu(LU.val());
          const Index * plucol(LU.col_ind());
          const Index * plurow(LU.row_ptr());

          const DT_ * pa(A.val());
          const Index * pacol(A.col_ind());
          const Index * parow(A.row_ptr());

          const Index n(LU.rows());
          Index k;
          // initialize LU array to A
          for (Index i(0); i < n; ++i)
          {
            k = parow[i];
            for (Index j(plurow[i]); j < plurow[i + 1]; ++j)
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
      } // function copy-entries
    };


    /**
     * \brief ILU(p)-Preconditioner for ELL-matrices.
     *
     * This class specializes the ILU(p)-Preconditioner for ELL-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class ILUPreconditioner<Algo::Generic, SparseMatrixELL<Mem_, DT_>,
                            DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixELL<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const SparseMatrixELL<Mem_, DT_> & A;
      SparseMatrixELL<Mem_, DT_> LU;

    public:
      virtual ~ILUPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] p level of fillin
       *           if p = 0, the layout of A is used for the ILU-decomposition
       *
       * Creates a ILU preconditioner to the given matrix and level of fillin
       */
      ILUPreconditioner(const SparseMatrixELL<Mem_, DT_> & A, const int p) :
        A(A)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        if (p == 0)
        {
          SparseMatrixELL<Mem_, DT_> tLU(A.layout());
          LU = tLU;

          copy_entries(false);
        }
        else
        {
          symbolic_lu_factorisation(p);
          copy_entries();
        }

        create_lu();
      }

      /**
       * \brief Constructor
       *
       * param[in] LU external LU-matrix
       *           the LU-decomposition is not calculated internally;
       *           so this preconditioner only solves \f$y \leftarrow L^{-1} U^{-1} x\f$
       *
       * Creates a ILU preconditioner to the given LU-decomposition
       */
      ILUPreconditioner(const SparseMatrixELL<Mem_, DT_> & LU) :
        A(LU),
        LU(LU)
      {
        if (LU.columns() != LU.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] layout An external layout for the LU-decomposition
       *
       * Creates a ILU preconditioner to the given matrix and layout
       */
      ILUPreconditioner(const SparseMatrixELL<Mem_, DT_> & A,
                        const SparseLayout<Mem_, SparseMatrixELL<Mem_, DT_>::LayoutType> & layout) :
        A(A),
        LU(layout)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        if (LU.columns() != LU.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        if (A.columns() != LU.columns())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrices have different sizes!");
        }

        copy_entries();

        create_lu();
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "ILU_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        out.copy(in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(LU.Ax());
        const Index * paj(LU.Aj());
        const Index stride(LU.stride());
        const Index * parl(LU.Arl());
        const Index n(A.rows());


        Index col;

        // __forward-insertion__
        // iteration over all rows
        for (Index i(0); i < n; ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          col = i;

          while (paj[col] < i)
          {
            pout[i] -= pval[col] * pout[paj[col]];
            col += stride;
          }
        }

        // __backward-insertion__
        // iteration over all rows
        for (Index i(n); i > 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          col = i + stride * (parl[i] - 1);

          while (paj[col] > i)
          {
            pout[i] -= pval[col] * pout[paj[col]];
            col -= stride;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }
      } // function apply

    private:
      void create_lu()
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/1014
         *    section 5.5.3; Algo 5.13.; page 132
         */

        const Index n(A.rows());
        DT_ * plu(LU.Ax());
        const Index * pj(LU.Aj());
        const Index * prl(LU.Arl());
        const Index stride(LU.stride());

        // integer work array of length n
        //   for saving the position of the diagonal-entries
        Index * pw = new Index[n];

        // iteration over all columns
        for (Index i(0); i < n; ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          // k -> \tilde a_{ik}
          Index k = i;
          while (pj[k] < i)
          {
            plu[k] /= plu[pw[pj[k]]];
            Index m(pw[pj[k]] + stride);
            // m -> \tilde a_{kj}
            for (Index j(k + stride); j < i + prl[i] * stride; j += stride)
            {
              // j -> \tilde a_{ij}
              while (m < pj[k] + prl[pj[k]] * stride)
              {
                if (pj[m] == pj[j])
                {
                  plu[j] -= plu[k] * plu[m];
                  m += stride;
                  break;
                }
                else if(pj[m] > pj[j])
                {
                  break;
                }
                m += stride;
              }
            }
            k += stride;
          }
          // save the position of the diagonal-entry
          pw[i] = k;
        }

        delete[] pw;
      } // function create_lu

      void symbolic_lu_factorisation(int p)
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/1014
         *    section 5.5.6; Algo 5.19.; page 142
         */

        const Index n(A.rows());
        const Index * paj(A.Aj());
        const Index * parl(A.Arl());
        const Index stride(A.stride());

        // type of list-entries
        typedef std::pair<int, Index> PAIR_;
        // list for saving the non-zero entries of L
        std::list<PAIR_> * ll = new std::list<PAIR_>[n];

        // vector for saving the iterators to the diag-entries
        std::list<PAIR_>::iterator * pldiag = new std::list<PAIR_>::iterator[n];
        // auxillary-iterators
        std::list<PAIR_>::iterator it, it1, it2;

        Index col, col2;
        int l, l2, neues_level;

        // fill list with non-zero entries of A
        // and save iterators to the diag- and last entry of each row
        for (Index row(0); row < n; ++row)
        {
          std::list<PAIR_> & ll_row (ll[row]);
          for (Index k(0); k < parl[row]; ++k)
          {
            col = paj[row  + stride * k];

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
        // calculate maximal number of columns per row
        Index num_cols_per_row(0);
        Index nnz(0);

        for (Index i(0); i < n; ++i)
        {
          num_cols_per_row = Math::max(num_cols_per_row, Index(ll[i].size()));
          nnz += Index(ll[i].size());
        }

        DenseVector<Mem_, DT_> LUx(num_cols_per_row * stride);
        DenseVector<Mem_, Index> LUj(num_cols_per_row * stride);
        DenseVector<Mem_, Index> LUrl(n);
        Index * pluj(LUj.elements());
        Index * plurl(LUrl.elements());

        Index used_elements(nnz);
        Index k1(0);

        for (Index i(0); i < n; ++i)
        {
          k1 = 0;
          for (it = ll[i].begin(); it != ll[i].end(); ++it, ++k1)
          {
            pluj[i + k1 * stride] = it->second;
          }
          plurl[i] = Index(ll[i].size());
        }

        SparseMatrixELL<Mem_, DT_> tLU(n, n, stride, num_cols_per_row,
                                       used_elements, LUx, LUj , LUrl);
        LU = tLU;

        delete[] ll;
        delete[] pldiag;
      } // symbolic_lu_factorisation

      void copy_entries(bool check = true)
      {
        if (check == false)
        {
          DT_ * plu(LU.Ax());
          const DT_ * pa(A.Ax());
          const Index data_length(LU.stride() * LU.num_cols_per_row());

          // initialize LU array to A
          for (Index i(0); i < data_length; ++i)
          {
            plu[i] = pa[i];
          }
        }
        else
        {
          DT_ * plu(LU.Ax());
          const Index * pluj(LU.Aj());
          const Index * plurl(LU.Arl());

          const DT_ * pa(A.Ax());
          const Index * paj(A.Aj());
          const Index * parl(A.Arl());

          const Index n(LU.rows());
          const Index stride(A.stride());

          Index k, ctr;

          // iteration over all rows
          for (Index row(0); row < n; ++row)
          {
            k = row;
            ctr = 0;
            for (Index j(0); j < plurl[row]; ++j)
            {
              plu[row + j * stride] = DT_(0.0);
              while (ctr <= parl[row] && paj[k] <= pluj[row + j * stride])
              {
                if (pluj[row + j * stride] == paj[k])
                {
                  plu[row + j * stride] = pa[k];
                  ++ctr;
                  k = k + stride;
                  break;
                }
                ++ctr;
                k = k + stride;
              }
            }
          }
        }
      } // function copy_entries
    };


    /**
     * \brief ILU(p)-Preconditioner for COO-matrices.
     *
     * This class specializes the ILU(p)-Preconditioner for COO-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class ILUPreconditioner<Algo::Generic, SparseMatrixCOO<Mem_, DT_>,
                            DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCOO<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      ILUPreconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                        DenseVector<Mem_, DT_> > _precond;
    public:
      virtual ~ILUPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] p level of fillin
       *           if p = 0, the layout of A is used for the ILU-decomposition
       *
       * Creates a ILU preconditioner to the given matrix and level of fillin
       */
      ILUPreconditioner(const SparseMatrixCOO<Mem_, DT_> & A, const int p) :
        _precond(SparseMatrixCSR<Mem_, DT_> (A), p)
      {
      }

      /**
       * \brief Constructor
       *
       * param[in] LU external LU-matrix
       *           the LU-decomposition is not calculated internally;
       *           so this preconditioner only solves \f$y \leftarrow L^{-1} U^{-1} x\f$
       *
       * Creates a ILU preconditioner to the given LU-decomposition
       */
      ILUPreconditioner(const SparseMatrixCOO<Mem_, DT_> & LU) :
        _precond(SparseMatrixCSR<Mem_, DT_> (LU))
      {
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "ILU_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        _precond.apply(out, in);
      }
    };


    /**
     * \brief SOR-Preconditioner.
     *
     * This class represents a dummy class for the SOR-Preconditioner.
     *
     * \author Christoph Lohmann
     */
    template <typename Algo_, typename MT_, typename VT_>
    class SORPreconditioner : public Preconditioner<Algo_, MT_, VT_>
    {
    };


    /**
     * \brief SOR-Preconditioner for CSR-matrices.
     *
     * This class specializes the SOR-Preconditioner for CSR-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class SORPreconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                            DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const SparseMatrixCSR<Mem_, DT_> & A;
      const DT_ _omega;

    public:
      virtual ~SORPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] omega A parameter of the preconditioner (default omega = 0.7)
       *
       * Creates a SOR preconditioner to the given matrix and parameter
       */
      SORPreconditioner(const SparseMatrixCSR<Mem_, DT_> & A,
                        DT_ omega = DT_(0.7)) :
        A(A),
        _omega(omega)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
      {
        return "SOR_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pa_val(A.val());
        const Index * pa_col_ind(A.col_ind());
        const Index * pa_row(A.row_ptr());

        Index row;

        // iteration over all rows
        for (Index i(0); i < A.rows(); ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          row = pa_row[i];

          while(pa_col_ind[row] < i)
          {
            pout[i] -= _omega * pa_val[row] * pout[pa_col_ind[row]];
            ++row;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pa_val[row];
        }

        out.template scale<Algo::Generic>(out, DT_(1.0) / _omega);
      }
    };


    /**
     * \brief SOR-Preconditioner for COO-matrices.
     *
     * This class specializes the SOR-Preconditioner for COO-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class SORPreconditioner<Algo::Generic, SparseMatrixCOO<Mem_, DT_>,
                            DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCOO<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const SparseMatrixCOO<Mem_, DT_> & A;
      const DT_ _omega;

    public:
      virtual ~SORPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] omega A parameter of the preconditioner (default omega = 0.7)
       *
       * Creates a SOR preconditioner to the given matrix and parameter
       */
      SORPreconditioner(const SparseMatrixCOO<Mem_, DT_> & A,
                        DT_ omega = DT_(0.7)) :
        A(A),
        _omega(omega)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
      {
        return "SOR_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pa_val(A.val());
        const Index * pa_column(A.column());
        const Index * pa_row(A.row());

        Index row(0);

        // Notice, that the elements of A are sorted by the column.
        // So we can solve the system with only one iteration over the data-array

        // iteration over all rows
        for (Index i(0); i < A.rows(); ++i)
        {
          // go to the i-th column of A
          while (pa_row[row] < i)
          {
            ++row;
          }

          // iteration over all elements on the left side of the main-diagonal
          while (pa_column[row] < i)
          {
            pout[i] -= _omega * pa_val[row] * out(pa_column[row]);
            ++row;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pa_val[row];
        }

        out.template scale<Algo::Generic>(out, DT_(1.0) / _omega);
      }
    };


    /**
     * \brief SOR-Preconditioner for ELL-matrices.
     *
     * This class specializes the SOR-Preconditioner for ELL-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class SORPreconditioner<Algo::Generic, SparseMatrixELL<Mem_, DT_>,
                            DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixELL<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const SparseMatrixELL<Mem_, DT_> & A;
      const DT_ _omega;

    public:
      virtual ~SORPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] omega A parameter of the preconditioner (default omega = 0.7)
       *
       * Creates a SOR preconditioner to the given matrix and parameter
       */
      SORPreconditioner(const SparseMatrixELL<Mem_, DT_> & A,
                        DT_ omega = DT_(0.7)) :
        A(A),
        _omega(omega)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
      {
        return "SOR_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pa_val(A.Ax());
        const Index * paj(A.Aj());

        Index row;
        const Index stride(A.stride());

        // iteration over all rows
        for (Index i(0); i < A.rows(); ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          row = i;

          while(paj[row] < i)
          {
            pout[i] -= _omega * pa_val[row] * pout[paj[row]];
            row += stride;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pa_val[row];
        }

        out.template scale<Algo::Generic>(out, DT_(1.0) / _omega);
      }
    };


    /**
     * \brief SPAI-Preconditioner.
     *
     * This class represents a dummy class for the SPAI-Preconditioner.
     *
     * \author Christoph Lohmann
     */
    template <typename Algo_, typename MT_, typename VT_>
    class SPAIPreconditioner : public Preconditioner<Algo_, MT_, VT_>
    {
    };


    /**
     * \brief SSOR-Preconditioner.
     *
     * This class represents a dummy class for the SSOR-Preconditioner.
     *
     * \author Christoph Lohmann
     */
    template <typename Algo_, typename MT_, typename VT_>
    class SSORPreconditioner : public Preconditioner<Algo_, MT_, VT_>
    {
    };


    /**
     * \brief SSOR-Preconditioner for CSR-matrices.
     *
     * This class specializes the SSOR-Preconditioner for CSR-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class SSORPreconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                             DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const SparseMatrixCSR<Mem_, DT_> & A;
      const DT_ _omega;
      DenseVector<Mem_, DT_> _diag;


    public:
      virtual ~SSORPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A The system-matrix
       * param[in] omega A parameter of the preconditioner (default omega = 1.3)
       *
       * Creates a SSOR preconditioner to the given matrix and parameter
       */
      SSORPreconditioner(const SparseMatrixCSR<Mem_, DT_> & A, const DT_ omega = DT_(1.3)) :
        A(A),
        _omega(omega),
        _diag(A.rows())
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        if (Math::abs(_omega - DT_(2.0)) < 1e-10)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "omega too close to 2!");
        }

        const Index n(A.rows());

        for (Index i(0) ; i < n ; ++i)
        {
          _diag(i, A(i, i));
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
      {
        return "SSOR_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(A.val());
        const Index * pcol_ind(A.col_ind());
        const Index * prow(A.row_ptr());
        const Index n(A.rows());

        Index col;

        // __forward-insertion__
        // iteration over all rows
        for (Index i(0); i < n; ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          col = prow[i];

          while(pcol_ind[col] < i)
          {
            pout[i] -= _omega * pval[col] * pout[pcol_ind[col]];
            ++col;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        out.template component_product<Algo::Generic>(_diag, in);

        // __backward-insertion__
        // iteration over all rows
        for (Index i(n); i > 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          col = prow[i+1]-1;

          while(pcol_ind[col] > i)
          {
            pout[i] -= _omega * pval[col] * pout[pcol_ind[col]];
            --col;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        out.template scale<Algo::Generic>(out, DT_(1.0) / (_omega * (DT_(2.0) - _omega)));
      } // function apply
    };


    /**
     * \brief SSOR-Preconditioner for ELL-matrices.
     *
     * This class specializes the SSOR-Preconditioner for ELL-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class SSORPreconditioner<Algo::Generic, SparseMatrixELL<Mem_, DT_>,
                             DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixELL<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const SparseMatrixELL<Mem_, DT_> & A;
      const DT_ _omega;
      DenseVector<Mem_, DT_> _diag;

    public:
      virtual ~SSORPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A The system-matrix
       * param[in] omega A parameter of the preconditioner (default omega = 1.3)
       *
       * Creates a SSOR preconditioner to the given matrix and parameter
       */
      SSORPreconditioner(const SparseMatrixELL<Mem_, DT_> & A, const DT_ omega = DT_(1.3)) :
        A(A),
        _omega(omega),
        _diag(A.rows())
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }

        if (Math::abs(_omega - DT_(2.0)) < 1e-10)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "omega too close to 2!");
        }

        const Index n(A.rows());

        for (Index i(0) ; i < n ; ++i)
        {
          _diag(i, A(i, i));
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
      {
        return "SSOR_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        // copy in-vector to out-vector
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pval(A.Ax());
        const Index * paj(A.Aj());
        const Index stride(A.stride());
        const Index * parl(A.Arl());
        const Index n(A.rows());


        Index col;

        // __forward-insertion__
        // iteration over all rows
        for (Index i(0); i < n; ++i)
        {
          // iteration over all elements on the left side of the main-diagonal
          col = i;

          while (paj[col] < i)
          {
            pout[i] -= _omega * pval[col] * pout[paj[col]];
            col += stride;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        out.template component_product<Algo::Generic>(_diag, in);

        // __backward-insertion__
        // iteration over all rows
        for (Index i(n); i > 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          col = i + stride * (parl[i] - 1);

          while (paj[col] > i)
          {
            pout[i] -= _omega * pval[col] * pout[paj[col]];
            col -= stride;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pval[col];
        }

        out.template scale<Algo::Generic>(out, DT_(1.0) / (_omega * (DT_(2.0) - _omega)));
      } // function apply
    };


    /**
     * \brief SSOR-Preconditioner for COO-matrices.
     *
     * This class specializes the SSOR-Preconditioner for COO-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class SSORPreconditioner<Algo::Generic, SparseMatrixCOO<Mem_, DT_>,
                             DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCOO<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      SSORPreconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                         DenseVector<Mem_, DT_> > _precond;
    public:
      virtual ~SSORPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A The system-matrix
       * param[in] omega A parameter of the preconditioner (default omega = 1.3)
       *
       * Creates a SSOR preconditioner to the given matrix and parameter
       */
      SSORPreconditioner(const SparseMatrixCOO<Mem_, DT_> & A, const DT_ omega = DT_(1.3)) :
        _precond(SparseMatrixCSR<Mem_, DT_> (A), omega)
      {
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
      {
        return "SSOR_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector to be preconditioned.
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        _precond.apply(out, in); // TODO: Warum geht das hier auf einmal nicht mehr? Segmentation fault
      }
    };


    /**
     * \brief SPAI-Preconditioner for CSR-matrices.
     *
     * This class specializes the SPAI-Preconditioner for CSR-matrices.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_>
    class SPAIPreconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                             DenseVector<Mem_, DT_> >
      : public Preconditioner<Algo::Generic, SparseMatrixCSR<Mem_, DT_>,
                              DenseVector<Mem_, DT_> >
    {
    private:
      const SparseMatrixCSR<Mem_, DT_> & A;
      const SparseLayout<SparseMatrixCSR<Mem::Main, bool> > &layout;
      SparseMatrixCSR<Mem_, DT_> M;

      const DT_ _eps_res;
      const Index _fill_in;
      const Index _max_iter;
      const DT_ _eps_res_comp;
      const DT_ _max_rho;
      const bool _column_structure;

    public:
      virtual ~SPAIPreconditioner()
      {
      }
      /**
       * \brief Constructor
       *
       * param[in] A system-matrix
       * param[in] layout the initial layout of the approximate inverse \f$M \approx A^{-1}\f$
       * param[in] eps_res stopping-criterion for new fill-in: norm of residuum
       * param[in] fill_in stopping-criterion for new fill-in: maximal number of fill-in per column
       * param[in] max_iter maximal number of iterations for creating new fill-in
       * param[in] eps_res_comp criterion for accepting a residual-component
       * param[in] max_rho criterion for acceptiong a rho-component
       *
       * Creates a SPAI preconditioner to the given matrix
       */
      SPAIPreconditioner(const SparseMatrixCSR<Mem_, DT_> & A,
                         const SparseLayout<SparseMatrixCSR<Mem::Main, bool> > & layout,
                         DT_ eps_res = 1e6, Index fill_in = 10, Index max_iter = 10,
                         DT_ eps_res_comp = 1e-6, DT_ max_rho = 1e-4,
                         bool column_structure = false) :
        A(A),
        layout(layout),
        _eps_res(eps_res),
        _fill_in(fill_in),
        _max_iter(max_iter),
        _eps_res_comp(eps_res_comp),
        _max_rho(max_rho),
        _column_structure(column_structure)
      {
        if (A.columns() != A.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
        if (layout._scalar_index[1] != layout._scalar_index[2])
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Precon-layout is not square!");
        }
        if (A.columns() != layout._scalar_index[1])
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Precon-layout and matrix do not match!");
        }

        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/1014
         *    section 5.6.1; Algo 5.21.; page 154
         */

        Index n(A.rows());
        const Index * pacol(A.col_ind());
        const Index * parow(A.row_ptr());
        const DT_ * pa(A.val());

        const Index * playoutcol(layout._indices[0]);
        const Index * playoutrow(layout._indices[1]);

        // type of list-entries
        typedef std::pair<DT_, Index> PAIR_;

        std::vector<std::list<PAIR_> > a_columnwise(n);
        std::list<PAIR_> * m_columns = new std::list<PAIR_>[n];

        typename std::list<PAIR_>::iterator it_I, it_J;


        if (_column_structure == true)
        {
          // 0. __create column-structure of matrix A__
          for (Index i(0); i < n; ++i)
          {
            for (Index l(parow[i]); l < parow[i + 1]; ++l)
            {
              a_columnwise[pacol[l]].emplace_back(pa[l], i);
            }
          }
        }

        // Iteration over each row of \f$M \approx A^{-1}\f$
        for (Index k = 0; k < n; ++k)
        {
          // Allocate memory for indices I and J of the algorithms
          // J saves the entries of the matrix \f$M \approx A^{-1}\f$, too,
          // I saves the residual \f$r_k\f$, too.
          // TODO: Speichern des Residuums in I nicht zwingend noetig.
          std::list<PAIR_> & J (m_columns[k]);
          std::list<PAIR_> I;

          // 1. __get initial structure of the j-th column of M__
          // Iteration over each row of A to get the initial layout of the k-th column of M
          for (Index row = 0; row < n; ++row)
          {
            Index i = playoutrow[row];

            while (playoutcol[i] < k && i < playoutrow[row + 1] - 1)
            {
              ++i;
            }

            if (playoutcol[i] == k)
            {
              J.emplace_back(DT_(0.0), row);
            }
          }

          // save size of J
          Index nn = Index(J.size());

          // 2. __get row-indices I of matching matrix-entries__
          if (_column_structure == true)
          {
            for (it_J = J.begin(); it_J != J.end(); ++it_J)
            {
              Index col = it_J->second;
              it_I = I.begin();
              for(auto it_col = a_columnwise[col].begin();
                  it_col != a_columnwise[col].end(); ++it_col)
              {
                Index row = it_col->second;
                while (it_I != I.end() && it_I->second < row)
                {
                  ++it_I;
                }
                if (it_I == I.end() || it_I->second != row)
                {
                  I.emplace(it_I, DT_(0.0), row);
                }
              }
            }
          }
          else
          {
            // Iteration over each row of A to get row-indices I matching to J
            for (Index i(0); i < n; ++i)
            {
              it_J = J.begin();
              for (Index l(parow[i]); l < parow[i + 1]; ++l)
              {
                while (it_J != J.end() && it_J->second < pacol[l])
                {
                  ++it_J;
                }

                if (it_J != J.end() && it_J->second == pacol[l])
                {
                  I.emplace_back(DT_(0.0), i);
                  break;
                }

                if (it_J == J.end())
                {
                  break;
                }
              }
            }
          }

          // save size of I
          Index mm = Index(I.size());

          // allocate dynamic matrix for saving the transped matrices \f$QR(I,J)^\top\f$ and \f$A(I,J)^\top\f$
          // TODO: Matrix local koennte noch entfernt werden; Mat-Vec-Mult kann mit QR durchgefuehrt werden
          std::vector<std::vector<DT_> > local(nn, std::vector<DT_>(mm));
          std::vector<std::vector<DT_> > qr(nn, std::vector<DT_>(mm));

          // fill temporary matrix local with entries of A(I,J)
          if (_column_structure == true)
          {
            it_J = J.begin();
            for (Index j(0); j < nn; ++j, ++it_J)
            {
              Index col = it_J->second;
              it_I = I.begin();
              Index i = 0;
              for(auto it_col = a_columnwise[col].begin();
                  it_col != a_columnwise[col].end(); ++it_col)
              {
                while (it_I->second < it_col->second)
                {
                  qr[j][i] = DT_(0.0);
                  local[j][i] = DT_(0.0);
                  ++it_I;
                  ++i;
                }
                qr[j][i] = it_col->first;
                local[j][i] = it_col->first;
                ++it_I;
                ++i;
              }
              while (i < mm)
              {
                qr[j][i] = DT_(0.0);
                local[j][i] = DT_(0.0);
                ++it_I;
                ++i;
              }
            }
          }
          else
          {
            it_I = I.begin();
            for (Index i(0); i < mm; ++i, ++it_I)
            {
              Index l = parow[it_I->second];
              it_J = J.begin();
              for (Index j(0); j < nn; ++j, ++it_J)
              {
                while (pacol[l] < it_J->second && l < parow[it_I->second + 1] - 1)
                {
                  ++l;
                }
                if (pacol[l] == it_J->second)
                {
                  qr[j][i] = pa[l];
                  local[j][i] = pa[l];
                  ++l;
                }
                else
                {
                  qr[j][i] = DT_(0.0);
                  local[j][i] = DT_(0.0);
                }
              }
            }
          }

          // 3. __calculate the QR-decomposition of \f$A(I,J)\f$
          /**
           * This algorithm has been adopted by
           *    Walter Gander - Algorithms for the QR-Decomposition; April 1980 (paper)
           *    page 12
           */
          // allocate dynamic vector for saving the diagonal entries of R
          std::vector<DT_> d(nn);

          for (Index j(0); j < nn; ++j)
          {
            DT_ s = DT_(0.0);
            for (Index i(j); i < mm; ++i)
            {
              s += Math::sqr(qr[j][i]);
            }
            s = Math::sqrt(s);

            if (qr[j][j] > 0)
            {
              d[j] = -s;
            }
            else
            {
              d[j] = s;
            }

            DT_ fak = Math::sqrt(s * (s + Math::abs(qr[j][j])));
            qr[j][j] -= d[j];

            for (Index l(j); l < mm; ++l)
            {
              qr[j][l] /= fak;
            }

            for (Index i(j + 1); i < nn; ++i)
            {
              DT_ s = DT_(0.0);
              for (Index l = j; l < mm; ++l)
              {
                s += qr[j][l] * qr[i][l];
              }
              for (Index l(j); l < mm; ++l)
              {
                qr[i][l] -= qr[j][l] * s;
              }
            }
          }

          // 4.a __calculate m_k as the solution of the least-squares-problem \f$A(I,J)^\top \cdot A(I,J) \cdot m_k = A^T(I,J) \cdot e_k\f$__
          // calculate \f$e_k\f$
          // TODO: Wir koennten auf den Vektor e verzichten und direkt in J->first rechnen; dafuer sollte J ein Vektor sein
          std::vector<DT_> e(mm);
          it_I = I.begin();
          for (Index i(0); i < mm; ++i, ++it_I)
          {
            if (it_I->second == k)
            {
              e[i] = DT_(1.0);
              it_I->first = DT_(-1.0);
            }
            else
            {
              e[i] = DT_(0.0);
              it_I->first = DT_(0.0);
            }
          }

          // calculate \f$e_k = Q^T * e_k\f$
          for (Index j(0); j < nn; ++j)
          {
            DT_ s = DT_(0.0);
            for (Index l(j); l < mm; ++l)
            {
              s += qr[j][l] * e[l];
            }
            for (Index l(j); l < mm; ++l)
            {
              e[l] -= qr[j][l] * s;
            }
          }

          // solve \f$R^{-1} e_k k= e_k\f$
          for (Index i(nn); i > 0;)
          {
            --i;
            for (Index j(i+1); j < nn; ++j)
            {
              e[i] -= qr[j][i] * e[j];
            }
            e[i] /= d[i];
          }

          // save values of e_k in J->first
          it_J = J.begin();
          for (Index j(0); j < nn; ++j, ++it_J)
          {
            it_J->first = e[j];
          }

          // 4.b __calculate the residual \f$r_k = A(I,J) \cdot J\to first - e_k\f$__
          it_I = I.begin();
          for (Index i(0); i < mm; ++i, ++it_I)
          {
            it_J = J.begin();
            for (Index j(0); j < nn; ++j, ++it_J)
            {
              it_I->first += local[j][i] * it_J->first;
            }
          }

          // 4.c __calculate the norm of the residual \f$res = \|r_k\|\f$
          DT_ res(DT_(0.0));
          for (it_I = I.begin(); it_I != I.end(); ++it_I)
          {
            res += Math::sqr(it_I->first);
          }
          res = Math::sqrt(res);

          // sorted list of I
          std::list<PAIR_> I_sorted(I);
          I_sorted.emplace_back(DT_(0.0), n);

          Index iter;
          // 5. __While \f$res = \|r_k\| > eps\f$, the maximal number of iterations and fill-in is not reached, repeat...__
          for (iter = 0; iter < _max_iter; ++iter)
          {
            // termination condition
            if (res < _eps_res)
            {
              break;
            }

            // allocate memory for \f$ \rho \f$
            std::vector<DT_> rho(mm, DT_(0.0));

            // 5.a __Search indices \f$i \in I\f$ with \f$i \not\in J\f$ and \f$r_k(i) \not=0\f$ and calculate \f$\rho_i\f$__
            // notice: if iter > 1 J is not sorted
            it_I = I.begin();
            for (Index i(0); i < mm; ++i, ++it_I)
            {
              it_J = J.begin();
              while (it_J->second != it_I->second && std::next(it_J) != J.end())
              {
                ++it_J;
              }

              if (it_J->second != it_I->second && Math::abs(it_I->first) > _eps_res_comp)
              {
                if (_column_structure == true)
                {
                  DT_ s(DT_(0.0));
                  for(auto it_col = a_columnwise[i].begin();
                      it_col != a_columnwise[i].end(); ++it_col)
                  {
                    s += Math::sqr(it_col->first);
                    auto it = I.begin();
                    while(it->second != it_col->second && std::next(it) != I.end())
                    {
                      ++it;
                    }
                    if (it->second == it_col->second)
                    {
                      rho[i] += it->first * it_col->first;
                    }
                  }
                  rho[i] = Math::sqr(res) - Math::sqr(rho[i]) / s;
                }
                else
                {
                  Index j;
                  DT_ s(DT_(0.0));
                  for (Index l(0); l < n; ++l)
                  {
                    j = parow[l];
                    while (pacol[j] < i && j < parow[l+1] - 1)
                    {
                      ++j;
                    }
                    if (pacol[j] == i)
                    {
                      s += Math::sqr(pa[j]);
                      auto it = I.begin();
                      while (it->second != l && std::next(it) != I.end())
                      {
                        ++it;
                      }
                      if (it->second == l)
                      {
                        rho[i] += it->first * pa[j];
                      }
                    }
                  }
                  rho[i] = Math::sqr(res) - Math::sqr(rho[i]) / s;
                }
              }
            }

            // save the iterators to the last entries of I and J
            auto it_I_end = std::prev(I.end());
            auto it_J_end = std::prev(J.end());

            // 5.b __add promising entries of \f$ \tilde J \f$ to \f$ J \f$
            // TODO: Auswahlkriterium ueberarbeiten
            it_I = I.begin();
            bool first = true;
            for (Index i(0); i < mm; ++i, ++it_I)
            {
              if (rho[i] > _max_rho && J.size() < _fill_in)
              {
                if (first == true)
                {
                  J.emplace_back(DT_(0.0), it_I->second);
                  first = false;
                  it_J = std::next(it_J_end);
                }
                else
                {
                  if (it_J->second > it_I->second)
                  {
                    it_J = std::next(it_J_end);
                  }
                  while (it_J != J.end() && it_J->second < it_I->second)
                  {
                    ++it_J;
                  }
                  it_J = J.emplace(it_J, DT_(0.0), it_I->second);
                }
              }
            }

            // save new size of J
            Index nn_new(Index(J.size()));

            // we can stop if no new entries have been added
            if (nn_new == nn)
            {
              break;
            }

            // calculate new indices \f$ \tilde I \f$
            if (_column_structure == true)
            {
              for (it_J = std::next(it_J_end); it_J != J.end(); ++it_J)
              {
                Index col = it_J->second;
                it_I = I_sorted.begin();
                for(auto it_col = a_columnwise[col].begin();
                    it_col != a_columnwise[col].end(); ++it_col)
                {
                  Index row = it_col->second;
                  while (it_I != I_sorted.end() && it_I->second < row)
                  {
                    ++it_I;
                  }
                  if (it_I == I_sorted.end() || it_I->second != row)
                  {
                    I_sorted.emplace(it_I, DT_(0.0), row);
                    auto it = std::next(it_I_end);
                    while (it != I.end() && it->second < row)
                    {
                      ++it;
                    }
                    I.emplace(it, DT_(0.0), row);
                  }
                }
              }
            }
            else
            {
              it_I = I_sorted.begin();
              Index i(0);
              for (it_I = I_sorted.begin(); it_I != I_sorted.end(); ++it_I)
              {
                for (; i < it_I->second; ++i)
                {
                  it_J = std::next(it_J_end);
                  for (Index l(parow[i]); l < parow[i + 1]; ++l)
                  {
                    while (it_J != J.end() && it_J->second < pacol[l])
                    {
                      ++it_J;
                    }

                    if (it_J != J.end() && it_J->second == pacol[l])
                    {
                      I.emplace_back(DT_(0.0), i);
                      I_sorted.emplace(it_I, DT_(0.0), i);
                      break;
                    }

                    if (it_J == J.end())
                    {
                      break;
                    }
                  }
                }
                ++i;
              }
            }

            // save new size of I
            Index mm_new(Index(I.size()));

            // resize matrices qr and local
            // now the storage locations of qr and local are similar to an upper triangular matrix
            for (Index j(nn); j < nn_new; ++j)
            {
              qr.push_back(std::vector<DT_>(mm_new));
              local.push_back(std::vector<DT_>(mm_new));
            }

            // fill temporary matrices local and qr with entries of A(I,J)
            it_I = I.begin();
            for (Index i(0); i < mm_new; ++i, ++it_I)
            {
              Index l = parow[it_I->second];
              it_J = std::next(it_J_end);
              for (Index j(nn); j < nn_new; ++j, ++it_J)
              {
                while (pacol[l] < it_J->second && l < parow[it_I->second + 1] - 1)
                {
                  ++l;
                }
                if (pacol[l] == it_J->second)
                {
                  qr[j][i] = pa[l];
                  local[j][i] = pa[l];
                }
                else
                {
                  qr[j][i] = DT_(0.0);
                  local[j][i] = DT_(0.0);
                }
              }
            }

            // mulitply last rows of matrix qr with \f$Q^\top\f$
            for (Index k(nn); k < nn_new; ++k)
            {
              // calculate \f$e_k = Q^\top \cdot e_k\f$
              for (Index j(0); j < nn; ++j)
              {
                DT_ s = DT_(0.0);
                for (Index l(j); l < qr[j].size(); ++l)
                {
                  s += qr[j][l] * qr[k][l];
                }
                for (Index l(j); l < qr[j].size(); ++l)
                {
                  qr[k][l] -= qr[j][l] * s;
                }
              }
            }

            // 5.c __calculate the qr-decomposition of \f$A(\tilde I, \tilde J\f$__
            // notice: we only have to calculate the qr-decomposition for the new columns \f$\tilde J \setminus J\f$
            // resizethe dynamic vector for saving the diagonal entries of R
            for (Index j(nn); j < nn_new; ++j)
            {
              d.push_back(DT_(0.0));
            }

            for (Index j(nn); j < nn_new; ++j)
            {
              DT_ s = DT_(0.0);
              for (Index i(j); i < mm_new; ++i)
              {
                s += Math::sqr(qr[j][i]);
              }
              s = Math::sqrt(s);

              if (qr[j][j] > 0)
              {
                d[j] = -s;
              }
              else
              {
                d[j] = s;
              }

              DT_ fak = Math::sqrt(s * (s + Math::abs(qr[j][j])));
              qr[j][j] -= d[j];

              for (Index l(j); l < mm_new; ++l)
              {
                qr[j][l] /= fak;
              }

              for (Index i(j + 1); i < nn_new; ++i)
              {
                DT_ s = DT_(0.0);
                for (Index l = j; l < mm_new; ++l)
                {
                  s += qr[j][l] * qr[i][l];
                }
                for (Index l(j); l < mm_new; ++l)
                {
                  qr[i][l] -= qr[j][l] * s;
                }
              }
            }

            // 5.d __calculate m_k as the solution of the least-squares-problem \f$A(I,J)^\top \cdot A(I,J) \cdot m_k = A^T(I,J) \cdot e_k\f$__
            // TODO: Wir koennten auf den Vektor e verzichten und direkt in J->first rechnen; dafuer sollte J ein Vektor sein
            // calculate \f$e_k\f$
            for (Index j(mm); j < mm_new; ++j)
            {
              e.push_back(DT_(0.0));
            }
            it_I = I.begin();
            for (Index i(0); i < mm_new; ++i, ++it_I)
            {
              if (it_I->second == k)
              {
                e[i] = DT_(1.0);
                it_I->first = DT_(-1.0);
              }
              else
              {
                e[i] = DT_(0.0);
                it_I->first = DT_(0.0);
              }
            }

            // calculate \f$e_k = Q^\top \cdot e_k\f$
            for (Index j(0); j < nn_new; ++j)
            {
              DT_ s = DT_(0.0);
              for (Index l(j); l < qr[j].size(); ++l)
              {
                s += qr[j][l] * e[l];
              }
              for (Index l(j); l < qr[j].size(); ++l)
              {
                e[l] -= qr[j][l] * s;
              }
            }

            // solve \f$R^{-1} e_k = e_k\f$
            for (Index i(nn_new); i > 0;)
            {
              --i;
              for (Index j(i+1); j < nn_new; ++j)
              {
                e[i] -= qr[j][i] * e[j];
              }
              e[i] /= d[i];
            }

            // save values of e_k in J->first
            it_J = J.begin();
            for (Index j(0); j < nn_new; ++j, ++it_J)
            {
              it_J->first = e[j];
            }

            // 5.e __calculate the residual \f$r_k = A(I,J) \cdot J\to first - e_k\f$__
            it_J = J.begin();
            for (Index j(0); j < nn_new; ++j, ++it_J)
            {
              it_I = I.begin();
              for (Index i(0); i < qr[j].size(); ++i, ++it_I)
              {
                it_I->first += local[j][i] * it_J->first;
              }
            }

            // 5.f __calculate the norm of the residual \f$res = \|r_k\|\f$__
            res = DT_(0.0);
            for (it_I = I.begin(); it_I != I.end(); ++it_I)
            {
              res += Math::sqr(it_I->first);
            }
            res = Math::sqrt(res);

            // set old dimensions of I and J to the new one
            mm = mm_new;
            nn = nn_new;
          }
        } // end for-loop over each row of \f$M \approx A^{-1}\f$

        // Create matrix M with calculated entries
        // calculate number of used elements
        Index nnz = 0;
        for (Index k(0); k < n; ++k)
        {
          nnz += Index(m_columns[k].size());
        }

        DenseVector<Mem_, Index> col_ind(nnz);
        DenseVector<Mem_, DT_> val(nnz);
        DenseVector<Mem_, Index> row_ptr(n+1);
        Index nz(0);
        row_ptr(0,0);

        for (Index k(0); k < n; ++k)
        {
          for (it_J = m_columns[k].begin(); it_J != m_columns[k].end(); ++it_J)
          {
            col_ind(nz, it_J->second);
            val(nz, it_J->first);
            ++nz;
          }
          row_ptr(k+1, nz);
        }

        SparseMatrixCSR<Mem_, DT_> tM(n, n, col_ind, val, row_ptr);
        M = tM;

        delete[] m_columns;
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
      {
        return "SPAI_Preconditioner";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector, which is applied to the preconditioning
       */
      virtual void apply(DenseVector<Mem_, DT_> & out,
                         const DenseVector<Mem_, DT_> & in)
      {
        if (in.elements() == out.elements())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Input- and output-vectors must be different!");
        }

        const Index n(M.rows());
        const Index * pmcol(M.col_ind());
        const Index * pmrow(M.row_ptr());
        const DT_ * pm(M.val());
        const DT_ * pin(in.elements());
        DT_ * pout(out.elements());

        for (Index i(0); i < n; i++)
        {
          pout[i] = DT_(0.0);
        }

        for (Index i(0); i < n; i++)
        {
          for (Index c(pmrow[i]); c < pmrow[i + 1]; c++)
          {
            pout[pmcol[c]] += pm[c] * pin[i];
          }
        }
      }
    };
  }// namespace LAFEM
} // namespace FEAST




#endif // KERNEL_LAFEM_PRECONDITIONER_HPP
