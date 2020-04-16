// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_LEGACY_PRECONDITIONERS_HPP
#define KERNEL_SOLVER_LEGACY_PRECONDITIONERS_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <vector>

namespace FEAT
{
  namespace Solver
  {
    /**
     * Supported sparse precon types.
     */
    enum class SparsePreconType
    {
      pt_none = 0,
      pt_file,
      pt_diagonal,
      pt_jacobi,
      pt_gauss_seidel,
      pt_polynomial,
      pt_ilu,
      pt_sor,
      pt_ssor,
      pt_spai
    };

    /**
     * \brief Preconditioner base class
     *
     */
    template <typename MT_, typename VT_>
    class Preconditioner
    {
    public:
      virtual ~Preconditioner()
      {
      }

      virtual void apply(VT_ & out, const VT_ & in) = 0;
    };

    /// \cond internal
    namespace Intern
    {
      template <typename MT_, typename VT_>
      class SPAIPreconditionerMTdepending : public Preconditioner<MT_, VT_>
      {
      };

      template <typename DT_, typename IT_>
      class SPAIPreconditionerMTdepending<LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>,
                                          LAFEM::DenseVector<Mem::Main, DT_, IT_> >
        : public Preconditioner<LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>,
                                LAFEM::DenseVector<Mem::Main, DT_, IT_> >
      {
      private:
        typedef std::pair<DT_, IT_> PAIR_;
        typedef LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_> MatrixType;
        typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorType;

      protected:
        const MatrixType & _A;
        const LAFEM::SparseLayout<Mem::Main, IT_, MatrixType::layout_id> _layout;
        const Index _m;
        MatrixType _M;
        std::list<PAIR_> * _m_columns;
        std::vector<std::list<PAIR_> > _a_columnwise;

      public:
        SPAIPreconditionerMTdepending(const MatrixType & A, const Index m) :
          _A(A),
          _layout(_A.layout()),
          _m(m),
          _a_columnwise(_A.rows())
        {
        }

        SPAIPreconditionerMTdepending(const MatrixType & A,
                                      LAFEM::SparseLayout<Mem::Main, IT_, MatrixType::layout_id> && layout) :
          _A(A),
          _layout(std::move(layout)),
          _m(std::numeric_limits<Index>::max()),
          _a_columnwise(_A.rows())
        {
        }

        virtual ~SPAIPreconditionerMTdepending()
        {
        }

        void create_initial_m_columns ()
        {
          const Index n(_A.rows());

          const IT_ * playoutcol(_layout.get_indices().at(0));
          const IT_ * playoutrow(_layout.get_indices().at(1));

          for (Index i(0); i < n; ++i)
          {
            for (IT_ l(playoutrow[i]); l < playoutrow[i + 1]; ++l)
            {
              _m_columns[playoutcol[l]].emplace_back(DT_(0.0), IT_(i));
            }
          }
        } // function create_initial_m_columns

        void create_a_columnwise ()
        {
          const DT_ * pval(_A.val());
          const IT_ * pcol_ind(_A.col_ind());
          const IT_ * prow_ptr(_A.row_ptr());
          const Index n(_A.columns());

          for (Index i(0); i < n; ++i)
          {
            for (IT_ l(prow_ptr[i]); l < prow_ptr[i + 1]; ++l)
            {
              _a_columnwise[pcol_ind[l]].emplace_back(pval[l], IT_(i));
            }
          }
        } // function create_a_columnwise

        void create_m_transpose ()
        {
          const Index n(_A.rows());

          // calculate number of used elements
          Index nnz = 0;
          for (Index k(0); k < n; ++k)
          {
            nnz += Index(_m_columns[k].size());
          }

          LAFEM::DenseVector<Mem::Main, DT_, IT_> val(nnz);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> col_ind(nnz);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> row_ptr(n+1);
          DT_ * pval = val.elements();
          IT_ * pcol_ind = col_ind.elements();
          IT_ * prow_ptr = row_ptr.elements();
          IT_ nz(0);
          prow_ptr[0] = 0;

          for (Index k(0); k < n; ++k)
          {
            for (auto it_J = _m_columns[k].begin(); it_J != _m_columns[k].end(); ++it_J)
            {
              pcol_ind[nz] = it_J->second;
              pval[nz] = it_J->first;
              ++nz;
            }
            prow_ptr[k+1] = nz;
          }

          _M = MatrixType(n, n, col_ind, val, row_ptr);
        } // function create_m_transpose

        void create_m()
        {
          const Index n(_A.rows());

          LAFEM::DenseVector<Mem::Main, IT_, IT_> trow_ptr(n + 1, IT_(0));
          IT_ * ptrow_ptr(trow_ptr.elements());

          ptrow_ptr[0] = 0;

          Index used_elements(0);
          for (Index i(0); i < n; ++i)
          {
            used_elements += Index(_m_columns[i].size());
            for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
            {
              ++ptrow_ptr[it_J->second + 1];
            }
          }

          for (Index i(1); i < n - 1; ++i)
          {
            ptrow_ptr[i + 1] += ptrow_ptr[i];
          }

          LAFEM::DenseVector<Mem::Main, IT_, IT_> tcol_ind(used_elements);
          LAFEM::DenseVector<Mem::Main, DT_, IT_> tval(used_elements);
          IT_ * ptcol_ind(tcol_ind.elements());
          DT_ * ptval(tval.elements());

          for (IT_ i(0); i < IT_(n); ++i)
          {
            for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
            {
              const IT_ l(it_J->second);
              const IT_ j(ptrow_ptr[l]);
              ptval[j] = it_J->first;
              ptcol_ind[j] = i;
              ++ptrow_ptr[l];
            }
          }

          for (Index i(n); i > 0; --i)
          {
            ptrow_ptr[i] = ptrow_ptr[i - 1];
          }
          ptrow_ptr[0] = 0;

          _M = MatrixType(n, n, tcol_ind, tval, trow_ptr);
        } // function create_m

        void create_m_without_new_entries ()
        {
          const Index n(_A.rows());

          if (_m != std::numeric_limits<Index>::max())
          {
            const Index used_elements(n * (1 + 2 * _m) - _m * (_m + 1));

            LAFEM::DenseVector<Mem::Main, DT_, IT_> val(used_elements);
            LAFEM::DenseVector<Mem::Main, IT_, IT_> col_ind(used_elements);
            LAFEM::DenseVector<Mem::Main, IT_, IT_> row_ptr(n+1);
            DT_ * pval(val.elements());
            IT_ * pcol_ind(col_ind.elements());
            IT_ * prow_ptr(row_ptr.elements());

            prow_ptr[0] = IT_(0);

            IT_ k(0);
            for (IT_ i(0); i < IT_(n); ++i)
            {
              for (IT_ l((_m > i) ? 0 : i - IT_(_m)); l < Math::min(IT_(n), IT_(_m) + i + IT_(1)); ++l)
              {
                pcol_ind[k] = l;
                ++k;
              }
              prow_ptr[i+1] = k;
            }

            for (Index i(0); i < n; ++i)
            {
              for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
              {
                const IT_ tmp(std::min(it_J->second, IT_(_m)));
                pval[prow_ptr[it_J->second] + i - it_J->second + tmp] = it_J->first;
              }
            }

            _M = MatrixType(n, n, col_ind, val, row_ptr);
          }
          else
          {
            _M = MatrixType(_layout);

            DT_ * pval(_M.val());
            const IT_ * pcol_ind(_M.col_ind());
            const IT_ * prow_ptr(_M.row_ptr());

            for (Index i(0); i < n; ++i)
            {
              for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
              {
                IT_ k = prow_ptr[it_J->second];

                while (pcol_ind[k] != i)
                {
                  ++k;
                }

                pval[k] = it_J->first;
              }
            }
          }
        } // function create_m_without_new_entries

        void apply_m_transpose (VectorType & out, const VectorType & in)
        {
          const Index n(_M.rows());
          const IT_ * pmcol(_M.col_ind());
          const IT_ * pmrow(_M.row_ptr());
          const DT_ * pm(_M.val());
          const DT_ * pin(in.elements());
          DT_ * pout(out.elements());

          for (Index i(0); i < n; i++)
          {
            pout[i] = DT_(0.0);
          }

          for (Index i(0); i < n; i++)
          {
            for (IT_ c(pmrow[i]); c < pmrow[i + 1]; c++)
            {
              pout[pmcol[c]] += pm[c] * pin[i];
            }
          }
        } // function apply_m_transpose
      };


      template <typename DT_, typename IT_>
      class SPAIPreconditionerMTdepending<LAFEM::SparseMatrixCOO<Mem::Main, DT_, IT_>,
                                          LAFEM::DenseVector<Mem::Main, DT_, IT_> >
        : public Preconditioner<LAFEM::SparseMatrixCOO<Mem::Main, DT_, IT_>,
                                LAFEM::DenseVector<Mem::Main, DT_, IT_> >
      {
      private:
        typedef std::pair<DT_, IT_> PAIR_;
        typedef LAFEM::SparseMatrixCOO<Mem::Main, DT_, IT_> MatrixType;
        typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorType;

      protected:
        const MatrixType & _A;
        const LAFEM::SparseLayout<Mem::Main, IT_, MatrixType::layout_id> _layout;
        const Index _m;
        MatrixType _M;
        std::list<PAIR_> * _m_columns;
        std::vector<std::list<PAIR_> > _a_columnwise;

      public:
        SPAIPreconditionerMTdepending(const MatrixType & A, const Index m) :
          _A(A),
          _layout(_A.layout()),
          _m(m),
          _a_columnwise(_A.rows())
        {
        }

        SPAIPreconditionerMTdepending(const MatrixType & A,
                                      LAFEM::SparseLayout<Mem::Main, IT_, MatrixType::layout_id> && layout) :
          _A(A),
          _layout(std::move(layout)),
          _m(std::numeric_limits<Index>::max()),
          _a_columnwise(_A.rows())
        {
        }

        virtual ~SPAIPreconditionerMTdepending()
        {
        }

        void create_initial_m_columns ()
        {
          const IT_ * playoutcol(_layout.get_indices().at(1));
          const IT_ * playoutrow(_layout.get_indices().at(0));
          const Index used_elements(_layout.get_scalar_index().at(3));

          for (Index i(0); i < used_elements; ++i)
          {
            _m_columns[playoutcol[i]].emplace_back(DT_(0.0), playoutrow[i]);
          }
        } // function create_initial_m_columns

        void create_a_columnwise ()
        {
          const DT_ * pa(_A.get_elements().at(0));
          const IT_ * pacol(_A.get_indices().at(1));
          const IT_ * parow(_A.get_indices().at(0));
          const Index used_elements(_A.get_scalar_index().at(3));

          for (Index i(0); i < used_elements; ++i)
          {
            _a_columnwise[pacol[i]].emplace_back(pa[i], parow[i]);
          }
        } // function create_a_aolumnwise

        void create_m_transpose ()
        {
          const Index n(_A.rows());

          // calculate number of used elements
          Index nnz = 0;
          for (Index k(0); k < n; ++k)
          {
            nnz += Index(_m_columns[k].size());
          }

          LAFEM::DenseVector<Mem::Main, DT_, IT_> val(nnz);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> col_ind(nnz);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> row_ind(nnz);
          DT_ * pval = val.elements();
          IT_ * pcol_ind = col_ind.elements();
          IT_ * prow_ind = row_ind.elements();
          Index nz(0);

          for (IT_ k(0); k < IT_(n); ++k)
          {
            for (auto it_J = _m_columns[k].begin(); it_J != _m_columns[k].end(); ++it_J)
            {
              pcol_ind[nz] = it_J->second;
              prow_ind[nz] = k;
              pval[nz] = it_J->first;
              ++nz;
            }
          }

          _M = MatrixType(n, n, row_ind, col_ind, val);
        } // function create_m_transpose

        void create_m()
        {
          const Index n(_A.rows());

          LAFEM::DenseVector<Mem::Main, IT_, IT_> trow_ptr(n + 1, IT_(0));
          IT_ * ptrow_ptr(trow_ptr.elements());

          Index used_elements(0);
          for (Index i(0); i < n; ++i)
          {
            used_elements += Index(_m_columns[i].size());
            for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
            {
              ++ptrow_ptr[it_J->second + 1];
            }
          }
          ptrow_ptr[0] = IT_(0);

          for (Index i(1); i < n - 1; ++i)
          {
            ptrow_ptr[i + 1] += ptrow_ptr[i];
          }

          LAFEM::DenseVector<Mem::Main, IT_, IT_> tcol_ind(used_elements);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> trow_ind(used_elements);
          LAFEM::DenseVector<Mem::Main, DT_, IT_> tval(used_elements);
          IT_ * ptcol_ind(tcol_ind.elements());
          IT_ * ptrow_ind(trow_ind.elements());
          DT_ * ptval(tval.elements());

          for (IT_ i(0); i < IT_(n); ++i)
          {
            for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
            {
              const IT_ l(it_J->second);
              const IT_ j(ptrow_ptr[l]);
              ptval[j] = it_J->first;
              ptcol_ind[j] = i;
              ptrow_ind[j] = it_J->second;
              ++ptrow_ptr[l];
            }
          }

          _M = MatrixType(n, n, trow_ind, tcol_ind, tval);
        } // function create_m

        void create_m_without_new_entries ()
        {
          const Index n(_A.rows());

          if (_m != std::numeric_limits<Index>::max())
          {
            const Index used_elements(n * (1 + 2 * _m) - _m * (_m + 1));

            LAFEM::DenseVector<Mem::Main, DT_, IT_> val(used_elements);
            LAFEM::DenseVector<Mem::Main, IT_, IT_> col_ind(used_elements);
            LAFEM::DenseVector<Mem::Main, IT_, IT_> row_ind(used_elements);
            LAFEM::DenseVector<Mem::Main, IT_, IT_> row_ptr(n+1);
            DT_ * pval(val.elements());
            IT_ * pcol_ind(col_ind.elements());
            IT_ * prow_ind(row_ind.elements());
            IT_ * prow_ptr(row_ptr.elements());

            prow_ptr[0] = IT_(0);

            IT_ k(0);
            for (IT_ i(0); i < IT_(n); ++i)
            {
              for (IT_ l((_m > i) ? 0 : i - IT_(_m)); l < Math::min(IT_(n), IT_(_m) + i + IT_(1)); ++l)
              {
                pcol_ind[k] = l;
                prow_ind[k] = i;
                ++k;
              }
              prow_ptr[i+1] = k;
            }

            for (Index i(0); i < n; ++i)
            {
              for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
              {
                const IT_ tmp(std::min(it_J->second, IT_(_m)));
                pval[prow_ptr[it_J->second] + i - it_J->second + tmp] = it_J->first;
              }
            }

            _M = MatrixType(n, n, row_ind, col_ind, val);
          }
          else
          {
            _M = MatrixType(_layout);

            for (Index i(0); i < n; ++i)
            {
              for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
              {
                _M(Index(it_J->second), i, it_J->first);
              }
            }
          }
        } // function create_m_without_new_entries

        void apply_m_transpose (VectorType & out, const VectorType & in)
        {
          const Index used_elements(_M.used_elements());
          const Index n(_M.rows());
          const DT_ * pval(_M.val());
          const IT_ * pcol(_M.column_indices());
          const IT_ * prow(_M.row_indices());
          const DT_ * pin(in.elements());
          DT_ * pout(out.elements());

          for (Index i(0); i < n; i++)
          {
            pout[i] = DT_(0.0);
          }

          for (Index i(0); i < used_elements; i++)
          {
            pout[pcol[i]] += pval[i] * pin[prow[i]];
          }
        } // function apply_m_transpose
      };


      template <typename DT_, typename IT_>
      class SPAIPreconditionerMTdepending<LAFEM::SparseMatrixELL<Mem::Main, DT_, IT_>,
                                          LAFEM::DenseVector<Mem::Main, DT_, IT_> >
        : public Preconditioner<LAFEM::SparseMatrixELL<Mem::Main, DT_, IT_>,
                                LAFEM::DenseVector<Mem::Main, DT_, IT_> >
      {
      private:
        typedef std::pair<DT_, IT_> PAIR_;
        typedef LAFEM::SparseMatrixELL<Mem::Main, DT_, IT_> MatrixType;
        typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorType;

      protected:
        const MatrixType & _A;
        const LAFEM::SparseLayout<Mem::Main, IT_, MatrixType::layout_id> _layout;
        const Index _m;
        MatrixType _M;
        std::list<PAIR_> * _m_columns;
        std::vector<std::list<PAIR_> > _a_columnwise;

      public:
        SPAIPreconditionerMTdepending(const MatrixType & A, const Index m) :
          _A(A),
          _layout(_A.layout()),
          _m(m),
          _a_columnwise(_A.rows())
        {
        }

        SPAIPreconditionerMTdepending(const MatrixType & A,
                                      LAFEM::SparseLayout<Mem::Main, IT_, MatrixType::layout_id> && layout) :
          _A(A),
          _layout(std::move(layout)),
          _m(std::numeric_limits<Index>::max()),
          _a_columnwise(_A.rows())
        {
        }

        virtual ~SPAIPreconditionerMTdepending()
        {
        }

        void create_initial_m_columns ()
        {
          const IT_ n((IT_(_A.rows())));

          const IT_ * playoutcol_ind(_layout.get_indices().at(0));
          const IT_ * playoutcs(_layout.get_indices().at(1));
          const IT_ * playoutrl(_layout.get_indices().at(3));
          const IT_ C((IT_(_layout.get_scalar_index().at(3))));

          for (IT_ i(0); i < n; ++i)
          {
            for (IT_ l(playoutcs[i/C] + i%C); l < playoutcs[i/C] + i%C + playoutrl[i] * C; l += C)
            {
              _m_columns[playoutcol_ind[l]].emplace_back(DT_(0.0), i);
            }
          }
        } // function create_initial_m_columns

        void create_a_columnwise ()
        {
          const IT_ n((IT_(_A.rows())));
          const DT_ * pval(_A.val());
          const IT_ * pacol_ind(_A.col_ind());
          const IT_ * pacs(_A.cs());
          const IT_ * parl(_A.rl());
          const IT_ C((IT_(_A.C())));

          for (IT_ i(0); i < n; ++i)
          {
            for (IT_ l(pacs[i/C] + i%C); l < pacs[i/C] + i%C + parl[i] * C; l += C)
            {
              _a_columnwise[pacol_ind[l]].emplace_back(pval[l], i);
            }
          }
        } // function create_a_aolumnwise

        void create_m_transpose ()
        {
          const IT_ n((IT_(_A.rows())));
          const IT_ C((IT_(_A.C())));

          // calculate cl-array and fill rl-array
          Index num_of_chunks(Index(ceil(n / float(C))));
          LAFEM::DenseVector<Mem::Main, IT_, IT_> mcl(num_of_chunks, IT_(0));
          LAFEM::DenseVector<Mem::Main, IT_, IT_> mcs(num_of_chunks + 1);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> mrl(_A.rows());
          IT_ * pmcl(mcl.elements());
          IT_ * pmcs(mcs.elements());
          IT_ * pmrl(mrl.elements());

          Index nnz(0);

          for (IT_ i(0); i < n; ++i)
          {
            pmrl[i] = IT_(_m_columns[i].size());
            pmcl[i/C] = Math::max(pmcl[i/C], pmrl[i]);
            nnz += Index(_m_columns[i].size());
          }

          // calculuate cs-array
          pmcs[0] = IT_(0);
          for (Index i(0); i < num_of_chunks; ++i)
          {
            pmcs[i+1] = pmcs[i] + IT_(C) * pmcl[i];
          }

          Index val_size = Index(pmcs[num_of_chunks]);

          LAFEM::DenseVector<Mem::Main, DT_, IT_> mval(val_size);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> mcol_ind(val_size);
          DT_ * pmval    (mval.elements());
          IT_ * pmcol_ind(mcol_ind.elements());

          for (IT_ i(0); i < n; ++i)
          {
            IT_ k(pmcs[i/C] + i%C);
            for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J, k+=C)
            {
              pmcol_ind[k] = it_J->second;
              pmval    [k] = it_J->first;
            }
            for (k+=C; k < pmcs[i/C + 1]; ++k)
            {
              pmcol_ind[k] = IT_(0);
              pmval    [k] = DT_(0);
            }
          }

          _M = MatrixType(_A.rows(), _A.columns(), nnz, mval, mcol_ind, mcs, mcl, mrl, _A.C());
        } // function create_m_transpose

        void create_m()
        {
          const IT_ n((IT_(_A.rows())));
          const IT_ C((IT_(_A.C())));

          // calculate cl-array and fill rl-array
          Index num_of_chunks(Index(ceil(n / float(C))));
          LAFEM::DenseVector<Mem::Main, IT_, IT_> mcl(num_of_chunks, IT_(0));
          LAFEM::DenseVector<Mem::Main, IT_, IT_> mcs(num_of_chunks + 1);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> mrl(_A.rows(), IT_(0));
          IT_ * pmcl(mcl.elements());
          IT_ * pmcs(mcs.elements());
          IT_ * pmrl(mrl.elements());

          Index nnz(0);

          for (IT_ i(0); i < n; ++i)
          {
            nnz += Index(_m_columns[i].size());
            for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
            {
              ++pmrl[it_J->second];
            }
          }

          for (Index i(0); i < n; ++i)
          {
            pmcl[i/C] = Math::max(pmcl[i/C], pmrl[i]);
            pmrl[i] = IT_(0);
          }

          // calculuate cs-array
          pmcs[0] = IT_(0);
          for (Index i(0); i < num_of_chunks; ++i)
          {
            pmcs[i+1] = pmcs[i] + IT_(C) * pmcl[i];
          }

          Index val_size = Index(pmcs[num_of_chunks]);

          LAFEM::DenseVector<Mem::Main, DT_, IT_> mval(val_size);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> mcol_ind(val_size);
          DT_ * pmval    (mval.elements());
          IT_ * pmcol_ind(mcol_ind.elements());


          for (IT_ i(0); i < n; ++i)
          {
            for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
            {
              const IT_ k(it_J->second);
              pmcol_ind[pmcs[k/C] + k%C + pmrl[k] * C] = i;
              pmval    [pmcs[k/C] + k%C + pmrl[k] * C] = it_J->first;
              ++pmrl[k];
            }
          }

          for (IT_ i(0); i < n; ++i)
          {
            for (IT_ k(pmcs[i/C] + i%C + pmrl[i] * C); k < pmcs[i/C + 1]; k+=C)
            {
              pmcol_ind[k] = IT_(0);
              pmval    [k] = DT_(0);
            }
          }

          _M = MatrixType(_A.rows(), _A.columns(), nnz, mval, mcol_ind, mcs, mcl, mrl, _A.C());
        } // function create_m

        void create_m_without_new_entries ()
        {
          const IT_ n((IT_(_A.rows())));
          if (_m != std::numeric_limits<Index>::max())
          {
            const Index used_elements(_A.rows() * (1 + 2 * _m) - _m * (_m + 1));
            const IT_ C((IT_(_A.C())));

            // calculate cl-array and fill rl-array
            Index num_of_chunks(Index(ceil(n / float(C))));
            LAFEM::DenseVector<Mem::Main, IT_, IT_> mcl(num_of_chunks, IT_(0));
            LAFEM::DenseVector<Mem::Main, IT_, IT_> mcs(num_of_chunks + 1);
            LAFEM::DenseVector<Mem::Main, IT_, IT_> mrl(_A.rows());
            IT_ * pmcl(mcl.elements());
            IT_ * pmcs(mcs.elements());
            IT_ * pmrl(mrl.elements());

            for (IT_ i(0); i < n; ++i)
            {
              pmrl[i] = IT_(_m + 1 + std::min(std::min(i, IT_(_m)), n - i - IT_(1)));
              pmcl[i/C] = Math::max(pmcl[i/C], pmrl[i]);
            }

            // calculuate cs-array
            pmcs[0] = IT_(0);
            for (Index i(0); i < num_of_chunks; ++i)
            {
              pmcs[i+1] = pmcs[i] + C * pmcl[i];
            }

            Index val_size = Index(pmcs[num_of_chunks]);

            LAFEM::DenseVector<Mem::Main, DT_, IT_> mval(val_size);
            LAFEM::DenseVector<Mem::Main, IT_, IT_> mcol_ind(val_size);
            DT_ * pmval    (mval.elements());
            IT_ * pmcol_ind(mcol_ind.elements());

            for (IT_ i(0); i < n; ++i)
            {
              for (IT_ l(IT_((_m > i) ? 0 : i - _m)), k(pmcs[i/C] + i%C); l < Math::min(n, IT_(_m) + i + 1); ++l, k+= C)
              {
                pmcol_ind[k] = IT_(l);
              }
              for (IT_ k(pmcs[i/C] + i%C + pmrl[i] * C); k < pmcs[i/C + 1]; k+=C)
              {
                pmcol_ind[k] = IT_(0);
                pmval    [k] = DT_(0);
              }
            }

            for (IT_ i(0); i < n; ++i)
            {
              for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
              {
                const IT_ tmp(std::min(it_J->second, IT_(_m)));
                pmval[pmcs[it_J->second/C] + it_J->second%C + C * (i - it_J->second + tmp)] = it_J->first;
              }
            }

            _M = MatrixType(_A.rows(), _A.columns(), used_elements, mval, mcol_ind, mcs, mcl, mrl, _A.C());
          }
          else
          {
            _M = MatrixType(_layout);

            DT_ * pmval(_M.val());
            const IT_ * pmcol_ind(_M.col_ind());
            const IT_ * pmcs(_M.cs());
            const IT_ C((IT_(_M.C())));

            for (IT_ i(0); i < n; ++i)
            {
              for (auto it_J = _m_columns[i].begin(); it_J != _m_columns[i].end(); ++it_J)
              {
                IT_ k(pmcs[it_J->second/C] + it_J->second%C);

                while (pmcol_ind[k] != i)
                {
                  k += C;
                }

                pmval[k] = it_J->first;
              }
            }
          }
        } // function create_m_without_new_entries

        void apply_m_transpose (VectorType & out, const VectorType & in)
        {
          const IT_ n((IT_(_M.rows())));
          const IT_ C((IT_(_M.C())));
          const DT_ * pval(_M.val());
          const IT_ * pcol_ind(_M.col_ind());
          const IT_ * pcs(_M.cs());
          const IT_ * prl(_M.rl());
          const DT_ * pin(in.elements());
          DT_ * pout(out.elements());

          for (IT_ i(0); i < n; i++)
          {
            pout[i] = DT_(0.0);
          }

          for (IT_ i(0); i < n; i++)
          {
            for (IT_ c(pcs[i/C] + i%C); c < pcs[i/C] + i%C + prl[i] * C; c += C)
            {
              pout[pcol_ind[c]] += pval[c] * pin[i];
            }
          }
        } // function apply_m_transpose
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief SPAI-Preconditioner.
     *
     * This class represents the SPAI-Preconditioner \f$M \approx A^{-1}\f$.
     *
     * \author Christoph Lohmann
     */
    template <typename MT_, typename VT_>
    class SPAIPreconditioner : public Intern::SPAIPreconditionerMTdepending<MT_, VT_>
    {
    private:
      typedef typename MT_::DataType DT_;
      typedef typename MT_::MemType Mem_;
      typedef typename MT_::IndexType IT_;
      typedef std::pair<DT_, IT_> PAIR_;

      using Intern::SPAIPreconditionerMTdepending<MT_, VT_>::_A;
      using Intern::SPAIPreconditionerMTdepending<MT_, VT_>::_m;
      using Intern::SPAIPreconditionerMTdepending<MT_, VT_>::_layout;
      using Intern::SPAIPreconditionerMTdepending<MT_, VT_>::_M;
      using Intern::SPAIPreconditionerMTdepending<MT_, VT_>::_m_columns;
      using Intern::SPAIPreconditionerMTdepending<MT_, VT_>::_a_columnwise;


      const DT_ _eps_res;
      const Index _fill_in;
      const Index _max_iter;
      const DT_ _eps_res_comp;
      const DT_ _max_rho;
      const bool _transpose;

    public:
      /// Our datatype
      typedef typename MT_::DataType DataType;
      /// Our indextype
      typedef typename MT_::IndexType IndexType;
      /// Our memory architecture type
      typedef typename MT_::MemType MemType;
      /// Our vectortype
      typedef VT_ VectorType;
      /// Our matrixtype
      typedef MT_ MatrixType;
      /// Our used precon type
      const static SparsePreconType PreconType = SparsePreconType::pt_spai;

      virtual ~SPAIPreconditioner()
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] A system-matrix
       *
       * \param[in] m full band-matrix with m diagonals on either side as initial layout (\f$2m + 1\f$ bands)
       *
       * \param[in] max_iter maximal number of iterations for creating new fill-in (default max_iter = 10)
       *
       * \param[in] eps_res stopping-criterion for new fill-in: norm of residuum (default eps_res = 1e-2)
       *
       * \param[in] fill_in stopping-criterion for new fill-in: maximal number of fill-in per column (default fill_in = 10)
       *
       * \param[in] eps_res_comp criterion for accepting a residual-component (default eps_res_comp = 1e-3)
       *
       * \param[in] max_rho criterion for acceptiong a rho-component (default max_rho = 1e-3)
       *
       * \param[in] transpose If you do only want to calculate _M^T, set transpose = true (default transpose = false)
       *
       * Creates a SPAI preconditioner to the given matrix and the initial layout defined by a band-matrix with \f$2m + 1\f$ bands
       */
      SPAIPreconditioner(const MT_ & A,
                         const Index m,
                         const Index max_iter = 10,
                         const DT_ eps_res = 1e-2,
                         const Index fill_in = 10,
                         const DT_ eps_res_comp = 1e-3,
                         const DT_ max_rho = 1e-3,
                         const bool transpose = false) :
        Intern::SPAIPreconditionerMTdepending<MT_, VT_>(A, m),
        _eps_res(eps_res),
        _fill_in(fill_in),
        _max_iter(max_iter),
        _eps_res_comp(eps_res_comp),
        _max_rho(max_rho),
        _transpose(transpose)
      {
        if (_A.columns() != _A.rows())
        {
          XABORTM("Matrix is not square!");
        }

        const Index n(_A.columns());

        _m_columns = new std::list<PAIR_>[n];

        // __get initial structure of the j-th column of M__
        for (Index i(0); i < n; ++i)
        {
          for (Index l((m > i) ? 0 : i - m); l < Math::min(n, m + i + 1); ++l)
          {
            _m_columns[i].emplace_back(DT_(0.0), l);
          }
        }

        _create_M();
      } // constructor

      /**
       * \brief Constructor
       *
       * \param[in] A system-matrix
       *
       * \param[in] layout the initial layout of the approximate inverse \f$M \approx A^{-1}\f$
       *
       * \param[in] max_iter maximal number of iterations for creating new fill-in (default max_iter = 10)
       *
       * \param[in] eps_res stopping-criterion for new fill-in: norm of residuum (default eps_res = 1.0)
       *
       * \param[in] fill_in stopping-criterion for new fill-in: maximal number of fill-in per column (default fill_in = 10)
       *
       * \param[in] eps_res_comp criterion for accepting a residual-component (default eps_res_comp = 1e-3)
       *
       * \param[in] max_rho criterion for acceptiong a rho-component (default max_rho = 1e-3)
       *
       * \param[in] transpose If you do only want to calculate _M^T, set transpose = true (default transpose = false)
       *
       * Creates a SPAI preconditioner to the given matrix and given initial layout
       */
      SPAIPreconditioner(const MT_ & A,
                         //LAFEM::SparseLayout<Mem_, typename MT_::IndexType, MT_::layout_id> && layout,
                         const Index max_iter = 10,
                         const DT_ eps_res = DT_(1.0),
                         const Index fill_in = 10,
                         const DT_ eps_res_comp = DT_(1e-3),
                         const DT_ max_rho = DT_(1e-3),
                         const bool transpose = false) :
        Intern::SPAIPreconditionerMTdepending<MT_, VT_>(A, std::move(A.layout())),
        _eps_res(eps_res),
        _fill_in(fill_in),
        _max_iter(max_iter),
        _eps_res_comp(eps_res_comp),
        _max_rho(max_rho),
        _transpose(transpose)
      {
        if (_A.columns() != _A.rows())
        {
          XABORTM("Matrix is not square!");
        }
        if (_layout.get_scalar_index().at(1) != _layout.get_scalar_index().at(2))
        {
          XABORTM("Precon-layout is not square!");
        }
        if (_A.columns() != _layout.get_scalar_index().at(1))
        {
          XABORTM("Precon-layout and matrix do not match!");
        }

        _m_columns = new std::list<PAIR_>[_A.rows()];

        // __get initial structure of the j-th column of M__
        this->create_initial_m_columns();

        _create_M();
      } // constructor

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SPAI";
      }

      /**
       * \brief apply the preconditioner
       *
       * \param[out] out The preconditioner result.
       * \param[in] in The vector, which is applied to the preconditioning
       */
      virtual void apply(LAFEM::DenseVector<Mem_, DT_, IT_> & out,
                         const LAFEM::DenseVector<Mem_, DT_, IT_> & in) override
      {
        if (_max_iter > 0 && _transpose == true)
        {
          if (in.elements() == out.elements())
          {
            XABORTM("Input- and output-vectors must be different!");
          }

          this->apply_m_transpose(out, in);
        }
        else
        {
          _M.apply(out, in);
        }
      }

      /// Returns the actual spai matrix.
      const MT_ & get_M() const
      {
        return _M;
      }

    private:
      void _create_M()
      {
        /**
         * This algorithm has been adopted by
         *    Dominik Goeddeke - Schnelle Loeser (Vorlesungsskript) 2013/2014
         *    section 5.6.1; Algo 5.21.; page 154
         */

        Index n(_A.rows());
        Index mm, nn, mm_new, nn_new;
        DT_ res;
        std::vector<DT_> d(0);

        typename std::list<PAIR_>::iterator it_I, it_J, it_I_end, it_J_end;
        typename std::list<std::pair<Index, IT_> >::iterator it_I_sorted;

        // __create column-structure of matrix A__
        this->create_a_columnwise();

        // Iteration over each row of \f$M \approx A^{-1}\f$
        for (Index k(0); k < n; ++k)
        {
          nn = 0;
          mm = 0;

          // Allocate memory for indices I and J of the algorithms
          // J saves the entries of the matrix \f$M \approx A^{-1}\f$, too,
          // I saves the residual \f$r_k\f$, too.
          std::list<PAIR_> & J (_m_columns[k]);
          std::list<PAIR_> I;

          it_I_end = I.begin();
          it_J_end = J.begin();

          // save size of J
          nn_new = Index(J.size());

          // __get row-indices I of matching matrix-entries__
          for (it_J = J.begin(); it_J != J.end(); ++it_J)
          {
            IT_ col = it_J->second;
            it_I = I.begin();
            for (auto it_col = _a_columnwise[col].begin();
                 it_col != _a_columnwise[col].end(); ++it_col)
            {
              IT_ row = it_col->second;
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

          // save size of I
          mm_new = Index(I.size());

          // save sorted list of I
          std::list<std::pair<Index, IT_> > I_sorted(I.size());
          it_I_sorted = I_sorted.begin();
          it_I = I.begin();
          for (Index i(0); i < mm_new; ++i, ++it_I_sorted, ++it_I)
          {
            it_I_sorted->second = it_I->second;
            it_I_sorted->first = i;
          }

          // allocate dynamic matrix for saving the transposed matrices \f$QR(I,J)^\top\f$ and \f$A(I,J)^\top\f$
          std::vector<std::vector<DT_> > qr;
          std::vector<std::vector<DT_> > local;

          // __While \f$res = \|r_k\| > eps\f$, the maximal number of iterations and fill-in is not reached, repeat...__
          Index iter(0);
          while (true)
          {
            // resize matrices qr and local
            qr.resize(nn_new, std::vector<DT_>(mm_new, DT_(0.0)));
            local.resize(nn_new, std::vector<DT_>(mm_new, DT_(0.0)));

            // fill temporary matrix qr and local with entries of A(I,J)
            it_J = ((it_J_end == J.begin()) ? J.begin() : std::next(it_J_end));
            for (Index j(nn); j < nn_new; ++j, ++it_J)
            {
              IT_ col = it_J->second;
              it_I_sorted = I_sorted.begin();
              for (auto it_col = _a_columnwise[col].begin();
                   it_col != _a_columnwise[col].end(); ++it_col)
              {
                while (it_I_sorted->second < it_col->second)
                {
                  ++it_I_sorted;
                }
                qr[j][it_I_sorted->first] = it_col->first;
                local[j][it_I_sorted->first] = it_col->first;
              }
            }

            // mulitply last/new rows of matrix qr with \f$Q^\top\f$
            // (the resctriction of nn > 1 could not be correct, but prevents the creation of NaNs)
            if (nn > 1)
            {
              for (Index k1(nn); k1 < nn_new; ++k1)
              {
                // calculate \f$e_k1 = Q^\top \cdot e_k1\f$
                for (Index j(0); j < nn; ++j)
                {
                  DT_ s = DT_(0.0);
                  for (Index l(j); l < qr[j].size(); ++l)
                  {
                    s += qr[j][l] * qr[k1][l];
                  }
                  for (Index l(j); l < qr[j].size(); ++l)
                  {
                    qr[k1][l] -= qr[j][l] * s;
                  }
                }
              }
            }

            // __calculate the qr-decomposition of \f$A(\tilde I, \tilde J\f$__
            // notice: we only have to calculate the qr-decomposition for the new columns \f$\tilde J \setminus J\f$
            // resize the dynamic vector for saving the diagonal entries of R
            d.resize(nn_new);

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
                s = DT_(0.0);
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

            // __calculate m_k as the solution of the least-squares-problem \f$A(I,J)^\top \cdot A(I,J) \cdot m_k = A^T(I,J) \cdot e_k\f$__
            // calculate \f$e_k\f$
            std::vector<DT_> e(mm_new);
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

            // termination condition
            if (iter >= _max_iter || nn_new >= _fill_in)
            {
              break;
            }
            ++iter;

            // __calculate the residual \f$r_k = A(I,J) \cdot J\to first - e_k\f$__
            it_J = J.begin();
            for (Index j(0); j < nn_new; ++j, ++it_J)
            {
              it_I = I.begin();
              for (Index i(0); i < qr[j].size(); ++i, ++it_I)
              {
                it_I->first += local[j][i] * it_J->first;
              }
            }

            // __calculate the norm of the residual \f$res = \|r_k\|\f$
            res = DT_(0.0);
            for (it_I = I.begin(); it_I != I.end(); ++it_I)
            {
              res += Math::sqr(it_I->first);
            }
            res = Math::sqrt(res);

            // set old dimensions of I and J to the new one
            mm = mm_new;
            nn = nn_new;

            // termination condition
            if (res < _eps_res)
            {
              break;
            }

            // allocate memory for \f$ \rho \f$
            std::vector<std::pair<DT_, IT_> > rho(mm);

            // __Search indices \f$i \in I\f$ with \f$i \not\in J\f$ and \f$r_k(i) \not=0\f$ and calculate \f$\rho_i\f$__
            // notice: if iter > 1 J is not sorted
            it_I = I.begin();
            for (Index i(0); i < mm; ++i, ++it_I)
            {
              if (Math::abs(it_I->first) < _eps_res_comp)
              {
                rho[i].first = DT_(0.0);
                rho[i].second = it_I->second;
                continue;
              }

              it_J = J.begin();
              while (it_J->second != it_I->second && std::next(it_J) != J.end())
              {
                ++it_J;
              }
              if (it_J->second == it_I->second)
              {
                continue;
              }

              DT_ s(DT_(0.0));
              auto it = I.begin();
              for (auto it_col = _a_columnwise[it_I->second].begin();
                   it_col != _a_columnwise[it_I->second].end(); ++it_col)
              {
                s += Math::sqr(it_col->first);
                while (it->second != it_col->second && std::next(it) != I.end())
                {
                  ++it;
                }
                if (it->second == it_col->second)
                {
                  rho[i].first += it->first * it_col->first;
                }
              }
              rho[i].first = Math::sqr(res) - Math::sqr(rho[i].first) / s;
              rho[i].second = it_I->second;
            }

            // save the iterators to the last entries of I and J
            it_I_end = std::prev(I.end());
            it_J_end = std::prev(J.end());

            // __add promising entries of \f$ \tilde J \f$ to \f$ J \f$
            bool first = true;
            while (J.size() < _fill_in)
            {
              // search maximal value in rho
              DT_ max_val = DT_(0.0);
              Index max_ind(0);
              IT_ max_sec(0);

              for (Index i(0); i < mm; ++i)
              {
                if (rho[i].first > max_val)
                {
                  max_val = rho[i].first;
                  max_sec = rho[i].second;
                  max_ind = i;
                }
              }

              if (max_val > _max_rho)
              {
                rho[max_ind].first = DT_(0.0);
                if (first == true)
                {
                  J.emplace_back(DT_(0.0), max_sec);
                  first = false;
                  it_J = std::next(it_J_end);
                }
                else
                {
                  if (it_J->second > max_sec)
                  {
                    it_J = std::next(it_J_end);
                  }
                  while (it_J != J.end() && it_J->second < max_sec)
                  {
                    ++it_J;
                  }
                  it_J = J.emplace(it_J, DT_(0.0), max_sec);
                }
              }
              else
              {
                break;
              }
            }

            // save new size of J
            nn_new = Index(J.size());

            // we can stop if no new entries have been added
            if (nn_new == nn)
            {
              break;
            }

            // calculate new indices \f$ \tilde I \f$
            for (it_J = std::next(it_J_end); it_J != J.end(); ++it_J)
            {
              IT_ col = it_J->second;
              it_I_sorted = I_sorted.begin();
              for (auto it_col = _a_columnwise[col].begin();
                   it_col != _a_columnwise[col].end(); ++it_col)
              {
                IT_ row = it_col->second;
                while (it_I_sorted != I_sorted.end() && it_I_sorted->second < row)
                {
                  ++it_I_sorted;
                }
                if (it_I_sorted == I_sorted.end() || it_I_sorted->second != row)
                {
                  I_sorted.emplace(it_I_sorted, std::numeric_limits<Index>::max(), row);
                  it_I = std::next(it_I_end);
                  while (it_I != I.end() && it_I->second < row)
                  {
                    ++it_I;
                  }
                  I.emplace(it_I, DT_(0.0), row);
                }
              }
            }

            // save new size of I
            mm_new = Index(I.size());

            // fill sorted vector sorted_I with new entries of I
            for (it_I_sorted = I_sorted.begin(); it_I_sorted != I_sorted.end(); ++it_I_sorted)
            {
              if (it_I_sorted->first != std::numeric_limits<Index>::max())
              {
                continue;
              }
              it_I = std::next(it_I_end);
              Index i(mm);
              while (it_I->second != it_I_sorted->second)
              {
                ++it_I;
                ++i;
              }
              it_I_sorted->first = i;
            }
          }
        } // end for-loop over each row of \f$M \approx A^{-1}\f$

        // Create matrix M with calculated entries
        if (_max_iter > 0)
        {
          if (_transpose == true)
          {
            this->create_m_transpose();
          }
          else
          {
            this->create_m();
          }
        }
        else
        {
          this->create_m_without_new_entries();
        }

        delete[] _m_columns;
      } // function _create_m
    };


    /**
     * \brief Polynomial-Preconditioner.
     *
     * This class represents the Neumann-Polynomial-Preconditioner \f$M^{-1} = \sum_{k=0}^m (I - \tilde M^{-1}A)^k \tilde M^{-1}\f$
     *
     * \author Christoph Lohmann
     */
    template <typename MT_, typename VT_>
    class PolynomialPreconditioner : public Preconditioner<MT_, VT_>
    {
    private:
      typedef typename MT_::MemType Mem_;
      typedef typename MT_::DataType DT_;

      const MT_ & _A;                              // system-matrix
      const Index _m;                              // order m of preconditioner
      const Index _num_of_auxs;                    // number of auxilary-vectors
                                                   // 1 for NonePreconditioner
                                                   // 2 if out and in can be the same vectors for _precon.apply
                                                   // 3 else
      VT_ _aux1, _aux2, _aux3;                     // auxilary-vector
      Preconditioner<MT_, VT_> * _precond;

    public:
      /// Our datatype
      typedef typename MT_::DataType DataType;
      /// Our memory architecture type
      typedef typename MT_::MemType MemType;
      /// Our vectortype
      typedef VT_ VectorType;
      /// Our matrixtype
      typedef MT_ MatrixType;
      /// Our used precon type
      const static SparsePreconType PreconType = SparsePreconType::pt_polynomial;

      /**
       * \brief Destructor
       *
       * deletes the preconditionier _precond
       */
      virtual ~PolynomialPreconditioner()
      {
        delete _precond;
      }

      /**
       * \brief Constructor of Neumann-Polynomial-Preconditioner.
       *
       * \param[in] A system-matrix
       * \param[in] m order of the polynom
       * \param[in] precond A preconditioner used for the Polynomial preconditioner
       *
       * Creates a Polynomial preconditioner to the given matrix and the given order
       */
      PolynomialPreconditioner(const MT_ & A, Index m, Preconditioner<MT_, VT_> * precond) :
        _A(A),
        _m(m),
        _num_of_auxs(3),
        _aux1(_A.rows()),
        _aux2(_A.rows()),
        _aux3(_A.rows()),
        _precond(precond)
      {
        if (_A.columns() != _A.rows())
        {
          XABORTM("Matrix is not square!");
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
      virtual void apply(VT_ & out, const VT_ & in) override
      {
        /*
         * preconditioner is given by
         *   \f$ M^-1 = \left(I + (I - \tilde M^{-1}A) + ... + (I - \tilde M^{-1 }A)^m\right) \tilde M^{-1} \f$
         *
         * the preconditioner only works, if
         *   ||I - \tilde M^{-1} A||_2 < 1.
         */

        _precond->apply(out, in);
        _aux3.copy(out);

        for (Index i = 1; i <= _m; ++i)
        {
          _A.apply(_aux1, out);
          _precond->apply(_aux2, _aux1);
          out.axpy(out, _aux3);
          out.axpy(_aux2, out, DataType(-1.0));
        }
      } // function apply
    };
  }// namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_LEGACY_PRECONDITIONERS_HPP
