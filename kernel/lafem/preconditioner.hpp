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
#include <kernel/lafem/algorithm.hpp>
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
      static String type_name()
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
      JacobiPreconditioner(const MT_ & A, typename VT_::DataType damping) :
        _jac(A.rows())
      {
        if (A.columns() != A.rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");

        for (Index i(0) ; i < _jac.size() ; ++i)
          _jac(i, damping / A(i, i));
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
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
      DT_ damping;
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
      static String type_name()
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
            pout[i] -= pa_val[row] * pout[pa_col_ind[row]];
            ++row;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pa_val[row];
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
      DT_ damping;
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
       * Creates a Gauss preconditioner to the given matrix and damping-parameter
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
      static String type_name()
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
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pa_val(A.val());
        const Index * pa_column(A.column());
        const Index * pa_row(A.row());

        Index row = 0;

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
            pout[i] -= pa_val[row] * out(pa_column[row]);
            ++row;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pa_val[row];
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
      DT_ damping;
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
       * Creates a Gauss preconditioner to the given matrix and damping-parameter
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
      static String type_name()
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
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pa_val(A.Ax());
        const Index * paj(A.Aj());

        Index col;
        const Index stride(A.stride());

        // iteration over all rows
        for (Index i(0); i < A.rows(); ++i)
        {
          // search entry on the main-diagonal
          col = i;

          while(paj[col] < i)
          {
            col += stride;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pa_val[col];

          // iteration over all elements on the left side of the main-diagonal
          for (Index j(i+1); j < A.rows(); ++j)
          {
            col = j;
            while(paj[col] < i)
            {
              col += stride;
            }

            if (paj[col] == i)
            {
              pout[j] -= pa_val[col] * pout[i];
            }
          }
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
      VT_ aux;              // auxilary-vector
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
      static String type_name()
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
          copy(*pauxs[0], in);

          for (Index i = 1; i <= m; ++i)
          {
            pauxs[i%2]->template axpy<Algo::Generic>(DT_(1.0), in, *pauxs[(i-1)%2]);
            pauxs[i%2]->template axpy<Algo::Generic>(DT_(-1.0), A, *pauxs[(i-1)%2], *pauxs[i%2]);
          }
        }
        else
        {
          pauxs[0]->template component_product<Algo::Generic>(*pscale, in);
          out.template scale<Algo::Generic>(out, DT_(-1.0));

          for (Index i = 1; i <= m; ++i)
          {
            pauxs[i%2]->template axpy<Algo::Generic>(DT_(1.0), in, *pauxs[(i-1)%2]);
            pauxs[i%2]->template axpy<Algo::Generic>(*pscale, A, *pauxs[(i-1)%2], *pauxs[i%2]);
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
                        const SparseLayout<SparseMatrixCSR<Mem::Main, bool> > & layout) :
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
      static String type_name()
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
        copy(out, in);

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
            it2 = pldiag[col];
            while (std::next(it2) != ll[col].end())
            {
              ++it2;

              col2 = it2->second;
              l2 = it2->first;

              neues_level = 2*(int) n - l - l2 + 1;

              // if new entries must be created, find the correct position in the list
              if (neues_level <= p)
              {
                while (it1->second < col2)
                {
                  if (it1 == ll[row].end())
                  {
                    ll[row].emplace_back(neues_level, col2);
                    break;
                  }
                  ++it1;
                }

                if (it1->second != col2)
                {
                  it1 = ll[row].emplace(it1, neues_level, col2);
                }
              }
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
                        const SparseLayout<SparseMatrixELL<Mem::Main, bool> > & layout) :
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
      static String type_name()
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
        copy(out, in);

        // create pointers
        DT_ * pout(out.elements());
        const DT_ * pA_val(LU.Ax());
        const Index * pAj(LU.Aj());
        const Index stride(LU.stride());
        const Index * pArl(LU.Arl());
        const Index n(A.rows());

        Index col, ncol;

        // __forward-insertion__
        // iteration over all rows
        for (Index i(0); i < n; ++i)
        {
          // search entry on the main-diagonal
          col = i;

          while(pAj[col] < i)
          {
            col += stride;
          }

          // iteration over all cells under the main-diagonal
          for (Index j(i+1); j < n; ++j)
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

        // __backward-insertion__
        // iteration over all rows
        for (Index i(LU.rows()); i> 0;)
        {
          --i;

          // iteration over all elements on the right side of the main-diagonal
          ncol = pArl[i];
          col = i + stride * (ncol - 1);

          while (pAj[col] > i)
          {
            pout[i] -= pA_val[col] * pout[pAj[col]];
            col = col - stride;
          }

          // divide by the element on the main-diagonal
          pout[i] /= pA_val[col];
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
            it2 = pldiag[col];
            while (std::next(it2) != ll[col].end())
            {
              ++it2;

              col2 = it2->second;
              l2 = it2->first;

              neues_level = 2*(int) n - l - l2 + 1;

              // if new entries must be created, find the correct position in the list
              // and update if necessary the iterator of the last row-entry
              if (neues_level <= p)
              {
                while (it1->second < col2)
                {
                  if (it1 == ll[row].end())
                  {
                    ll[row].emplace_back(neues_level, col2);
                    break;
                  }
                  ++it1;
                }

                if (it1->second != col2)
                {
                  it1 = ll[row].emplace(it1, neues_level, col2);
                }
              }
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
                        DenseVector<Mem_, DT_> > precond;
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
        precond(SparseMatrixCSR<Mem_, DT_> (A), p)
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
        precond(SparseMatrixCSR<Mem_, DT_> (LU))
      {
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String type_name()
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
        precond.apply(out, in);
      }
    };
  }// namespace LAFEM
} // namespace FEAST




#endif // KERNEL_LAFEM_PRECONDITIONER_HPP
