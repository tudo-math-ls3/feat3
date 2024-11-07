// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/mean_filter.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/cudss.hpp>
#include <kernel/solver/mkl_dss.hpp>
#include <kernel/solver/umfpack.hpp>

// includes, system
#include <deque>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Core class for conversion of Stokes systems to CSR format for direct solvers
     *
     * See the specialization of this class for details.
     *
     * \author Peter Zajac
     */
    template<
      typename SolverDT_, typename SolverIT_,
      typename MatrixA_, typename MatrixB_, typename MatrixD_,
      typename MatrixS_ = LAFEM::SparseMatrixCSR<typename MatrixA_::DataType, typename MatrixA_::IndexType>>
    class DirectStokesCore;

    /**
     * \brief Core class for conversion of BCSR Stokes systems to CSR format for direct solvers
     *
     * This class can be used to convert a Stokes saddle-point system of the form
     * \f[ \begin{bmatrix} A & B\\D & S \end{bmatrix} \cdot \begin{bmatrix}v\\p\end{bmatrix} = \begin{bmatrix}f\\g\end{bmatrix}\f]
     * consisting of BCSR-matrices A, B and D as well as a scalar matrix S, which is usually an empty matrix,
     * into a scalar CSR system of the form
     * \f[ \begin{bmatrix} \tilde{A} & \tilde{B} & Q & 0\\\tilde{D} & \tilde{S} & 0 & P\\Q^\top & 0 & 0 & 0\\0 & P^\top & 0 & 0 \end{bmatrix}
       \cdot \begin{bmatrix}v\\p\\q\\\lambda\end{bmatrix} = \begin{bmatrix}f\\g\\0\\0\end{bmatrix}\f]
       where
     * - \f$\tilde{A}\f$ corresponds to the matrix block \e A, which has been filtered by all velocity unit filters,
     *   i.e. each row that corresponds to a velocity unit filter DOF has been replaced by a unit row.
     * - \f$\tilde{B}\f$ corresponds to the matrix block \e B, which has been off-diagonally filtered by all velocity
     *   unit filters, i.e. each row that corresponds to a velocity unit filter DOF has been replaced by an empty row.
     * - \f$\tilde{D}\f$ corresponds to the matrix block \e D, which has been off-diagonally filtered by all pressure
     *   unit filters, i.e. each row that corresponds to a pressure unit filter DOF has been replaced by an empty row.
     * - \f$\tilde{S}\f$ corresponds to the matrix block \e S, which has been filtered by all pressure unit filters,
     *   i.e. each row that corresponds to a pressure unit filter DOF has been replaced by a unit row.
     * - \e Q corresponds to the <i>nv</i>-by-<i>nq</i> Lagrange multiplier matrices arising from the velocity slip filters,
     *   where \e nv denotes the total number of (scalar) velocity DOFs and \e nq denotes the total number
     *   of (blocked) velocity slip filter DOFs.
     * - \e P corresponds to the <i>np</i>-by-1 Lagrange multiplier vector arising from the pressure mean filter,
     *   where \e np denoted the total number of pressure DOFs
     * - \e q is the (vector-valued) Lagrange multiplier for the velocity slip boundary conditions.
     *   Note that the equation for this unknown only exists if there exists at least one non-empty velocity slip filter.
     * - \f$\lambda\f$ is the scalar Lagrange multiplier for the pressure mean filter.
     *   Note that the equation for this unknown only exists if there exists a non-empty pressure mean filter.
     *
     * This class is capable of handling and incorporating the following filters into the system:
     * - velocity unit filters for the incorporation of inflow and no-flow (aka no-slip) boundary conditions
     * - velocity slip filters for the incorporation of slip boundary conditions
     * - pressure unit filters for the incorporation of pressure Dirichlet conditions (rarely used)
     * - pressure mean filters for the incorporation of vanishing pressure integral mean condition
     *
     * \note The existence of pressure unit filters and a pressure mean filter are mutually exclusive, because they
     * both aim at resolving the issue of singular pressure in the case of a "no-outflow" system. However, this class
     * does not explicitly check for this mutual exclusivity.
     *
     * The primary use case of this class is as an internal class for the DirectStokesSolver implementations, however,
     * this class may also be used stand-alone, if one desires to do so. If you intend to use this class directly,
     * then the following points are important to remember:
     * - Before you call the init_symbolic() function, you have to ensure that (at least) the layouts of all matrices
     *   are initialized and that the velocity/pressure filters have also been initialized and set by calling the
     *   #set_filters() function.
     * - If the values of the velocity and/or pressure unit/slip filters change after the symbolic initialization has
     *   been performed, you have to update the filter contents by calling the #update_filters() function before
     *   performing the numerical initialization by calling the #init_numeric() function. Note that there is no harm
     *   in always updating the filters prior to a #init_numeric() call, even if it is redundant.
     *
     * \attention
     * As explained above, this class incorporates all boundary conditions and constraints, which are given by filters,
     * directly into the solver system, and as a direct consequence of this fact, the filter layouts are \b not allowed
     * to change between symbolic and numeric initialization! This means: if the indices of the velocity/pressure DOFs,
     * which are affected by unit or slip filters, change after symbolic initialization, then you \b must re-perform the
     * the symbolic initialization -- it is \b not sufficient to only re-perform the numeric initialization in this case!
     *
     * \tparam SolverDT_
     * The desired data type of the direct solver; this is usually double.
     *
     * \tparam SolverIT_
     * The desired index type of the direct solver
     *
     * \author Peter Zajac
     */
    template<typename SolverDT_, typename SolverIT_, typename DT_, typename IT_, int dim_>
    class DirectStokesCore<
      SolverDT_, SolverIT_,
      LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, dim_>,
      LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, 1>,
      LAFEM::SparseMatrixBCSR<DT_, IT_, 1, dim_>,
      LAFEM::SparseMatrixCSR<DT_, IT_>>
    {
    public:
      /// the type of the matrix block A
      typedef LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, dim_> MatrixTypeA;
      /// the type of the matrix block B
      typedef LAFEM::SparseMatrixBCSR<DT_, IT_, dim_, 1> MatrixTypeB;
      /// the type of the matrix block D
      typedef LAFEM::SparseMatrixBCSR<DT_, IT_, 1, dim_> MatrixTypeD;
      /// the type of the matrix block S
      typedef LAFEM::SparseMatrixCSR<DT_, IT_> MatrixTypeS;

      /// the type of the velocity vector
      typedef LAFEM::DenseVectorBlocked<DT_, IT_, dim_> VectorTypeV;
      /// the type of the pressure vector
      typedef LAFEM::DenseVector<DT_, IT_> VectorTypeP;

      /// the unit filter type for the velocity
      typedef LAFEM::UnitFilterBlocked<DT_, IT_, dim_> UnitFilterTypeV;
      /// the slip filter type for the velocity
      typedef LAFEM::SlipFilter<DT_, IT_, dim_> SlipFilterTypeV;
      /// the none filter type for the velocity
      typedef LAFEM::NoneFilterBlocked<DT_, IT_, dim_> NoneFilterTypeV;

      /// the unit filter type for the pressure
      typedef LAFEM::UnitFilter<DT_, IT_> UnitFilterTypeP;
      //typedef LAFEM::MeanFilter<DT_, IT_> MeanFilterTypeP;
      //typedef Global::MeanFilter<DT_, IT_> MeanFilterTypePG;
      /// the none filter type for the pressure
      typedef LAFEM::NoneFilter<DT_, IT_> NoneFilterTypeP;

      /// the data type to be used by the solver
      typedef SolverDT_ SolverDataType;
      /// the index type to be used by the solver
      typedef SolverIT_ SolverIndexType;

      /// the solver matrix type
      typedef LAFEM::SparseMatrixCSR<SolverDataType, SolverIndexType> SolverMatrixType;
      /// the solver vector type
      typedef LAFEM::DenseVector<SolverDataType, SolverIndexType> SolverVectorType;
      /// the solver filter type
      typedef LAFEM::NoneFilter<SolverDataType, SolverIndexType> SolverFilterType;

    protected:
      /// reference to matrix block A
      const MatrixTypeA& _matrix_a;
      /// reference to matrix block B
      const MatrixTypeB& _matrix_b;
      /// reference to matrix block D
      const MatrixTypeD& _matrix_d;
      /// pointer to matrix block S (may be nullptr or an empty matrix)
      const MatrixTypeS* _matrix_s;

      /// the solver matrix
      SolverMatrixType _solver_matrix;
      /// the solver filter (a NoneFilter)
      SolverFilterType _solver_filter;

      /// deque of velocity slip filters
      std::deque<SlipFilterTypeV> _slip_filters_v;
      /// deque of velocity unit filters
      std::deque<UnitFilterTypeV> _unit_filters_v;
      /// deque of pressure unit filters
      std::deque<UnitFilterTypeP> _unit_filters_p;
      /// pressure mean filter (may be nullptr or an empty vector if no pressure mean filter is to be applied)
      VectorTypeP _mean_filter_p_dual;

      /// have the filters been set yet?
      bool _filters_already_set = false;

      /// do we have to use the pressure mean filter?
      bool _have_mean_filter_p = false;

      /// do we need a main diagonal in S?
      bool _need_main_diagonal = false;

      /// unit filter masks for velocity and pressure
      std::vector<int> _unit_mask_v, _unit_mask_p;

      /// slip filter matrix Q entries
      std::vector<IT_> _slip_row_ptr, _slip_col_idx, _slip_fil_dqu, _slip_fil_idx;

      /// type for filter counting
      typedef std::array<std::size_t, 3u> FilterUpdateCounter;

    public:
      /**
       * \brief Constructor
       *
       * This constructor sets the references to the sub-matrices. Note that the matrices do not need to be initialized
       * when this constructor is called.
       *
       * \param[in] matrix_a, matrix_b, matrix_d
       * The \resident references to the matrix blocks A, B and D.
       *
       * \param[in] matrix_s
       * A \resident pointer to the matrix block S.
       * This may be set to \c nullptr if the S block is empty, which is the usual case.
       */
      explicit DirectStokesCore(const MatrixTypeA& matrix_a, const MatrixTypeB& matrix_b,
        const MatrixTypeD& matrix_d, const MatrixTypeS* matrix_s = nullptr) :
        _matrix_a(matrix_a),
        _matrix_b(matrix_b),
        _matrix_d(matrix_d),
        _matrix_s(matrix_s)
      {
      }

      /// virtual destructor
      virtual ~DirectStokesCore()
      {
      }

      /// no copy allowed
      DirectStokesCore(const DirectStokesCore&) = delete;
      /// no copy allowed
      DirectStokesCore& operator=(const DirectStokesCore&) = delete;

      /**
       * \brief Specifies whether the solver matrix needs to contain all main diagonal entries
       *
       * \param[in] need_main_diagonal
       * Specifies whether the solver matrix structure should always contain all main diagonal entries.
       * This is required by some solvers, e.g. by the MKL DSS solver.
       */
      void set_need_main_diagonal(bool need_it)
      {
        this->_need_main_diagonal = need_it;
      }

      /**
       * \brief Returns a const reference to the solver matrix.
       */
      const SolverMatrixType& get_solver_matrix() const
      {
        return this->_solver_matrix;
      }

      /**
       * \brief Returns a const reference to the solver filter, which is always an empty NoneFilter
       */
      const SolverFilterType& get_solver_filter() const
      {
        return this->_solver_filter;
      }

      /**
       * \brief Creates a new solver vector, whose contents are left uninitialized.
       */
      SolverVectorType create_solver_vector() const
      {
        return this->_solver_matrix.create_vector_r();
      }

      /**
       * \brief Uploads a velocity/pressure vector pair to a solver vector.
       *
       * \param[inout] vector
       * A \transient reference to the solver vector that receives the uploaded data.
       *
       * \param[in] vector_v, vector_p
       * The \transient references to the velocity/pressure vectors that are to be uploaded.
       */
      virtual void upload_vector(SolverVectorType& vector, const VectorTypeV& vector_v, const VectorTypeP& vector_p) const
      {
        const IT_ nv = IT_(vector_v.template size<LAFEM::Perspective::pod>());
        const IT_ np = IT_(vector_p.size());
        const IT_ nx = IT_(vector.size());

        XASSERTM(nx == IT_(this->_solver_matrix.rows()), "invalid solver vector size");
        XASSERTM(nv + np <= nx, "invalid velocity/pressure vector size");

        SolverDataType* vec = vector.elements();
        const DT_* vec_v = vector_v.template elements<LAFEM::Perspective::pod>();
        const DT_* vec_p = vector_p.elements();
        IT_ j(0u);
        // copy velocity DOFs
        for(IT_ i(0); i < nv; ++i, ++j)
          vec[j] = SolverDataType(vec_v[i]);
        // copy pressure DOFs
        for(IT_ i(0); i < np; ++i, ++j)
          vec[j] = SolverDataType(vec_p[i]);
        // format Lagrange multipliers
        for(; j < nx; ++j)
          vec[j] = SolverDataType(0);
      }

      /**
       * \brief Downloads a velocity/pressure vector pair from a solver vector.
       *
       * \param[in] vector
       * A \transient reference to the solver vector that provides the download.
       *
       * \param[inout] vector_v, vector_p
       * The \transient references to the velocity/pressure vectors that receive the downloaded data.
       */
      virtual void download_vector(const SolverVectorType& vector, VectorTypeV& vector_v, VectorTypeP& vector_p) const
      {
        const IT_ nv = IT_(vector_v.template size<LAFEM::Perspective::pod>());
        const IT_ np = IT_(vector_p.size());
        const IT_ nx = IT_(vector.size());

        XASSERTM(nx == IT_(this->_solver_matrix.rows()), "invalid solver vector size");
        XASSERTM(nv + np <= nx, "invalid velocity/pressure vector size");

        const SolverDataType* vec = vector.elements();
        DT_* vec_v = vector_v.template elements<LAFEM::Perspective::pod>();
        DT_* vec_p = vector_p.elements();
        IT_ j(0u);
        // copy velocity DOFs
        for(IT_ i(0); i < nv; ++i, ++j)
          vec_v[i] = DT_(vec[j]);
        // copy pressure DOFs
        for(IT_ i(0); i < np; ++i, ++j)
          vec_p[i] = DT_(vec[j]);
        // Lagrange multipliers are ignored
      }

      /**
       * \brief Sets the velocity and pressure filters
       *
       * This function must be called before init_symbolic() and all filter layouts must already have been assembled.
       *
       * \param[in] filter
       * A \resident reference to the filter pair.
       */
      template<typename FilterV_, typename FilterP_>
      void set_filters(const LAFEM::TupleFilter<FilterV_, FilterP_>& filter)
      {
        XASSERTM(!this->_filters_already_set, "Filter have already been set");
        this->_add_filter_velo(filter.template at<0>(), nullptr);
        this->_add_filter_pres(filter.template at<1>(), nullptr);
        this->_filters_already_set = true;
      }

      /**
       * \brief Updates the velocity and pressure filters
       *
       * This function can be called between init_symoblic() and init_numeric() if the filters have been reassembled.
       * Note that only the values of the filters are allowed to change, but their layouts must remain unchanged!
       *
       * \param[in] filter
       * A \resident reference to the filter pair.
       */
      template<typename FilterV_, typename FilterP_>
      void update_filters(const LAFEM::TupleFilter<FilterV_, FilterP_>& filter)
      {
        XASSERTM(this->_filters_already_set, "Filters have not been set yet");
        FilterUpdateCounter fuc{0u, 0u, 0u};
        this->_add_filter_velo(filter.template at<0>(), &fuc);
        this->_add_filter_pres(filter.template at<1>(), &fuc);
      }

      /**
       * \brief Performs the symbolic initialization
       *
       * This function computes the sparsity pattern of the solver matrix.
       */
      virtual void init_symbolic()
      {
        // get the number of DOFs
        const IT_ idim = IT_(dim_);
        const IT_ num_dofs_v = IT_(_matrix_a.rows());
        const IT_ num_dofs_p = IT_(_matrix_d.rows());

        XASSERT(_matrix_a.columns() == num_dofs_v);
        XASSERT(_matrix_b.columns() == num_dofs_p);
        XASSERT(_matrix_d.columns() == num_dofs_v);
        XASSERT(_matrix_b.rows() == num_dofs_v);
        XASSERT(_matrix_d.rows() == num_dofs_p);

        if(_matrix_s != nullptr)
        {
          XASSERT(_matrix_s->rows() == num_dofs_p);
          XASSERT(_matrix_s->columns() == num_dofs_p);
        }

        // is the mean filter dual vector non-empty?
        this->_have_mean_filter_p = !this->_mean_filter_p_dual.empty();

        // compute number of slip DOFs
        const IT_ num_slip_dofs = this->_compute_slip_matrix_layout();

        // mask velocity and pressure unit dofs
        this->_compute_unit_mask_v();
        this->_compute_unit_mask_p();

        // compute total number of scalar DOFs
        const IT_ num_total_dofs = idim*num_dofs_v + num_dofs_p + num_slip_dofs + IT_(_have_mean_filter_p ? 1 : 0);

        // get number of matrix NZEs
        const IT_ num_nze_a = IT_(_matrix_a.used_elements());
        const IT_ num_nze_b = IT_(_matrix_b.used_elements());
        const IT_ num_nze_d = IT_(_matrix_d.used_elements());
        const IT_ num_nze_s = IT_(_matrix_s != nullptr ? _matrix_s->used_elements() : 0u);

        // compute upper bound for non-zero entry count
        IT_ max_nze = idim*idim*num_nze_a + idim*(num_nze_b + num_nze_d) + num_nze_s;

        // Lagrange multiplier for velocity slip filters
        max_nze += IT_(2)*num_slip_dofs*idim;

        // Lagrange multiplier for pressure mean filter (when necessary)
        if(this->_have_mean_filter_p)
          max_nze += IT_(2)*num_dofs_p;

        // need a main diagonal?
        if(this->_need_main_diagonal)
          max_nze += num_dofs_p + IT_(2)*num_slip_dofs*idim + IT_(1);

        // allocate row-pointer and column index arrays
        std::vector<IT_> row_ptr, col_idx;
        row_ptr.resize(num_total_dofs + IT_(1u), ~IT_(0));
        col_idx.resize(max_nze, ~IT_(0));
        row_ptr[0u] = IT_(0u);

        // the system has the following structure:
        //   | A   B   Q   0 |
        //   | D   S   0   P |
        //   | Q^T 0   0   0 |
        //   | 0   P^T 0   0 |
        //
        // where:
        // - A, B, D and S are the original matrices from the Stokes system
        // - Q is the Lagrange multiplier for the velocity slip DOFs
        // - P is the Lagrange multiplier for the pressure mean filter

        // get all the matrix arrays
        const IT_* row_ptr_a = _matrix_a.row_ptr();
        const IT_* col_idx_a = _matrix_a.col_ind();
        const IT_* row_ptr_b = _matrix_b.row_ptr();
        const IT_* col_idx_b = _matrix_b.col_ind();
        const IT_* row_ptr_d = _matrix_d.row_ptr();
        const IT_* col_idx_d = _matrix_d.col_ind();
        const IT_* row_ptr_s = (_matrix_s != nullptr ? _matrix_s->row_ptr() : nullptr);
        const IT_* col_idx_s = (_matrix_s != nullptr ? _matrix_s->col_ind() : nullptr);

        IT_ row(0u), op(0u);

        // okay, loop over all rows of A/B/Q
        for(IT_ i(0); i < num_dofs_v; ++i)
        {
          // is this a unit row?
          if(_unit_mask_v[i] != 0)
          {
            // set the entire block to an identity block
            for(IT_ k(0); k < idim; ++k)
            {
              col_idx[op] = i*idim + k;
              row_ptr[++row] = ++op;
            }
            continue;
          }

          // not a unit row, so copy rows of A, B and Q
          for(IT_ j(0); j < idim; ++j)
          {
            // copy row of A
            for(IT_ ip(row_ptr_a[i]); ip < row_ptr_a[i+1]; ++ip)
              for(IT_ k(0); k < idim; ++k, ++op)
                col_idx[op] = idim* col_idx_a[ip] + k;
            // copy row of B
            for(IT_ ip(row_ptr_b[i]); ip < row_ptr_b[i+1]; ++ip, ++op)
              col_idx[op] = idim*num_dofs_v + col_idx_b[ip];
            // copy row of Q
            for(IT_ ip(_slip_row_ptr[i]); ip < _slip_row_ptr[i+1]; ++ip, ++op)
              col_idx[op] = idim*num_dofs_v + num_dofs_p + _slip_col_idx[ip];
            // set row pointer
            row_ptr[++row] = op;
          }
        }

        XASSERT(row == idim*num_dofs_v);

        // next, loop over all rows of D/S/P
        for(IT_ i(0); i < num_dofs_p; ++i)
        {
          // is this a unit row?
          if(_unit_mask_p[i] != 0)
          {
            col_idx[op] = idim*num_dofs_v + i;
            row_ptr[++row] = ++op;
            continue;
          }

          // not a unit row, so copy rows of D, S and P
          // copy row of D
          for(IT_ ip(row_ptr_d[i]); ip < row_ptr_d[i+1]; ++ip)
            for(IT_ k(0); k < idim; ++k, ++op)
              col_idx[op] = idim*col_idx_d[ip] + k;
          // copy row of S (if it exists)
          if(row_ptr_s != nullptr)
          {
            if(this->_need_main_diagonal)
            {
              // copy row i of S and insert a main diagonal entry if it is not present in S
              bool have_diag = false;
              for(IT_ ip(row_ptr_s[i]); ip < row_ptr_s[i+1]; ++ip, ++op)
              {
                if(!have_diag && (col_idx_s[ip] > i))
                {
                  col_idx[op++] = row;
                  have_diag = true;
                }
                col_idx[op] = idim*num_dofs_v + col_idx_s[ip];
                have_diag = have_diag || (col_idx_s[ip] == i);
              }
            }
            else
            {
              for(IT_ ip(row_ptr_s[i]); ip < row_ptr_s[i+1]; ++ip, ++op)
                col_idx[op] = idim*num_dofs_v + col_idx_s[ip];
            }
          }
          else if(this->_need_main_diagonal)
          {
            // insert a main diagonal entry
            col_idx[op++] = row;
          }
          // copy row of P (if we have a pressure mean filter)
          if(this->_have_mean_filter_p)
            col_idx[op++] = idim*num_dofs_v + num_slip_dofs + num_dofs_p;
          // set row pointer
          row_ptr[++row] = op;
        }

        // next, add all slip filter entries Q^T
        for(auto it(_slip_filters_v.begin()); it != _slip_filters_v.end(); ++it)
        {
          const IT_ nsx = IT_(it->used_elements());
          const IT_* idx = it->get_indices();
          for(IT_ j(0); j < nsx; ++j)
          {
            for(IT_ k(0); k < idim; ++k, ++op)
              col_idx[op] = idim*idx[j] + k;
            if(this->_need_main_diagonal)
              col_idx[op++] = row;
            // set row pointer
            row_ptr[++row] = op;
          }
        }

        // finally, add the pressure mean filter vector P^T (if it exists)
        if(this->_have_mean_filter_p)
        {
          for(IT_ j(0); j < num_dofs_p; ++j, ++op)
            col_idx[op] = idim*num_dofs_v + j;
          if(this->_need_main_diagonal)
            col_idx[op++] = row;
          row_ptr[++row] = op;
        }

        // the number of rows must match, but we may have less non-zeros due to unit filtering
        XASSERT(row == num_total_dofs);
        XASSERT(op <= max_nze);

        const Index num_nze = op;

        // ok, allocate the matrix
        this->_solver_matrix = SolverMatrixType(num_total_dofs, num_total_dofs, num_nze);

        // copy row pointer
        SolverIndexType* row_ptr_x = this->_solver_matrix.row_ptr();
        for(IT_ i(0); i <= num_total_dofs; ++i)
          row_ptr_x[i] = SolverIndexType(row_ptr[i]);

        // copy column indices
        SolverIndexType* col_idx_x = this->_solver_matrix.col_ind();
        for(IT_ i(0); i < num_nze; ++i)
          col_idx_x[i] = SolverIndexType(col_idx[i]);
      }

      /**
       * \brief Performs the numeric initialization
       *
       * This function computes the numerical entries of the solver matrix.
       */
      virtual void init_numeric()
      {
        // get the number of DOFs
        const IT_ num_dofs_v = IT_(_matrix_a.rows());
        const IT_ num_dofs_p = IT_(_matrix_d.rows());

        // the system has the following structure:
        //   | A   B   Q   0 |
        //   | D   S   0   P |
        //   | Q^T 0   0   0 |
        //   | 0   P^T 0   0 |
        //
        // where:
        // - A, B, D and S are the original matrices from the Stokes system
        // - Q is the Lagrange multiplier for the velocity slip DOFs
        // - P is the Lagrange multiplier for the pressure mean filter

        // get all the matrix arrays
        const IT_* row_ptr_a = _matrix_a.row_ptr();
        const auto* val_a    = _matrix_a.val();
        const IT_* row_ptr_b = _matrix_b.row_ptr();
        const auto* val_b    = _matrix_b.val();
        const IT_* row_ptr_d = _matrix_d.row_ptr();
        const auto* val_d    = _matrix_d.val();
        const IT_* row_ptr_s = (_matrix_s != nullptr ? _matrix_s->row_ptr() : nullptr);
        const IT_* col_idx_s = (_matrix_s != nullptr ? _matrix_s->col_ind() : nullptr);
        const auto* val_s    = (_matrix_s != nullptr ? _matrix_s->val() : nullptr);
        SolverDataType* val_x = _solver_matrix.val();

        IT_ op(0u);

        // okay, loop over all rows of A/B/Q
        for(IT_ i(0); i < num_dofs_v; ++i)
        {
          // is this a unit row?
          if(_unit_mask_v[i] != 0)
          {
            // set the entire block to an identity block
            for(int k(0); k < dim_; ++k, ++op)
              val_x[op] = SolverDataType(1);
            continue;
          }

          // not a unit row, so copy rows of A, B and Q
          for(int j(0); j < dim_; ++j)
          {
            // copy row of A
            for(IT_ ip(row_ptr_a[i]); ip < row_ptr_a[i+1]; ++ip)
              for(int k(0); k < dim_; ++k, ++op)
                val_x[op] = SolverDataType(val_a[ip][j][k]);
            // copy row of B
            for(IT_ ip(row_ptr_b[i]); ip < row_ptr_b[i+1]; ++ip, ++op)
              val_x[op] = SolverDataType(val_b[ip][j][0]);
            // copy row of Q
            for(IT_ ip(_slip_row_ptr[i]); ip < _slip_row_ptr[i+1]; ++ip, ++op)
            {
              // get the normal vector of this slip filter entry
              const auto& nu = this->_slip_filters_v.at(_slip_fil_dqu[ip]).get_values()[_slip_fil_idx[ip]];
              // set the normal vector component
              val_x[op] = SolverDataType(nu[j]);
            }
          }
        }

        // next, loop over all rows of D/S/P
        for(IT_ i(0); i < num_dofs_p; ++i)
        {
          // is this a unit row?
          if(_unit_mask_p[i] != 0)
          {
            val_x[op++] = SolverDataType(1);
            continue;
          }

          // not a unit row, so copy rows of D, S and P
          // copy row of D
          for(IT_ ip(row_ptr_d[i]); ip < row_ptr_d[i+1]; ++ip)
            for(int k(0); k < dim_; ++k, ++op)
              val_x[op] = SolverDataType(val_d[ip][0][k]);
          // copy row of S (if it exists)
          if(row_ptr_s != nullptr)
          {
            if(this->_need_main_diagonal)
            {
              bool have_diag = false;
              for(IT_ ip(row_ptr_s[i]); ip < row_ptr_s[i+1]; ++ip, ++op)
              {
                if(!have_diag && (col_idx_s[ip] > i))
                {
                  val_x[op++] = SolverDataType(0);
                  have_diag = true;
                }
                val_x[op] = SolverDataType(val_s[ip]);
                have_diag = have_diag || (col_idx_s[ip] == i);
              }
            }
            else
            {
              for(IT_ ip(row_ptr_s[i]); ip < row_ptr_s[i+1]; ++ip, ++op)
                val_x[op] = SolverDataType(val_s[ip]);
            }
          }
          else if(this->_need_main_diagonal)
          {
            // insert a main diagonal entry
            val_x[op++] = SolverDataType(0);
          }

          // copy row of P (if we have a pressure mean filter)
          if(this->_have_mean_filter_p)
          {
            const auto* v = _mean_filter_p_dual.elements();
            val_x[op++] = SolverDataType(v[i]);
          }
        }

        // next, add all slip filter entries Q^T
        for(auto it(_slip_filters_v.begin()); it != _slip_filters_v.end(); ++it)
        {
          const IT_ nsx = IT_(it->used_elements());
          const auto* nu = it->get_values();
          for(IT_ j(0); j < nsx; ++j)
          {
            for(int k(0); k < dim_; ++k, ++op)
              val_x[op] = SolverDataType(nu[j][k]);
            if(this->_need_main_diagonal)
              val_x[op++] = SolverDataType(0);
          }
        }

        // finally, add the pressure mean filter vector P^T (if it exists)
        if(this->_have_mean_filter_p)
        {
          const auto* v = _mean_filter_p_dual.elements();
          for(IT_ j(0); j < num_dofs_p; ++j, ++op)
            val_x[op] = SolverDataType(v[j]);
          if(this->_need_main_diagonal)
            val_x[op++] = SolverDataType(0);
        }

        // sanity check
        XASSERT(op == _solver_matrix.used_elements());
      }

      /// \brief Releases all data allocated during numerical initialization
      virtual void done_numeric()
      {
        // nothing to do here
      }

      /**
       * \brief Releases all data allocated during symbolic initialization
       *
       * This function also clears all internal filter structures, i.e. the filters have to be set
       * again manually by calling the set_filters() function before another call to init_symoblic()
       * can be performed.
       */
      virtual void done_symbolic()
      {
        _slip_row_ptr.clear();
        _slip_col_idx.clear();
        _slip_fil_dqu.clear();
        _slip_fil_idx.clear();
        _slip_filters_v.clear();
        _unit_filters_v.clear();
        _unit_filters_p.clear();
        _mean_filter_p_dual.clear();
        _unit_mask_v.clear();
        _unit_mask_p.clear();
        _have_mean_filter_p = false;
        _filters_already_set = false;
      }

    protected:
      /// adds or updates a velocity unit filter
      void _add_filter_velo(const UnitFilterTypeV& filter, FilterUpdateCounter* fuc)
      {
        if(fuc != nullptr)
        {
          // get the old unit filter
          UnitFilterTypeV& old_filter = this->_unit_filters_v.at(fuc->at(0u)++);

          // compare layout of old and new filter
          if(!old_filter.get_filter_vector().compare_layout(filter.get_filter_vector()))
            XABORTM("velocity unit filter layout has changed");

          // update the filter
          old_filter.clone(filter, LAFEM::CloneMode::Shallow);
        }
        else
        {
          // push new filter
          this->_unit_filters_v.push_back(filter.clone(LAFEM::CloneMode::Shallow));
        }
      }

      /// adds or updates a velocity slip filter
      void _add_filter_velo(const SlipFilterTypeV& filter, FilterUpdateCounter* fuc)
      {
        if(fuc != nullptr)
        {
          // get the old slip filter
          SlipFilterTypeV& old_filter = this->_slip_filters_v.at(fuc->at(1u)++);

          // compare layout of old and new filter
          if(!old_filter.get_filter_vector().compare_layout(filter.get_filter_vector()))
            XABORTM("velocity slip filter layout has changed");

          // update the filter
          old_filter.clone(filter, LAFEM::CloneMode::Shallow);
        }
        else
        {
          // push new filter
          this->_slip_filters_v.push_back(filter.clone(LAFEM::CloneMode::Shallow));
        }
      }

      /// adds or updates a velocity none filter (dummy function)
      void _add_filter_velo(const NoneFilterTypeV&, FilterUpdateCounter*)
      {
        // nothing to do here
      }

      /// adds or updates a velocity filter chain
      template<typename First_>
      void _add_filter_velo(const LAFEM::FilterChain<First_>& filter, FilterUpdateCounter* fuc)
      {
        this->_add_filter_velo(filter.first(), fuc);
      }

      /// adds or updates a velocity filter chain
      template<typename First_, typename Second_, typename... Rest_>
      void _add_filter_velo(const LAFEM::FilterChain<First_, Second_, Rest_...>& filter, FilterUpdateCounter* fuc)
      {
        this->_add_filter_velo(filter.first(), fuc);
        this->_add_filter_velo(filter.rest(), fuc);
      }

      /// adds or updates a velocity filter sequence
      template<typename SubFilter_>
      void _add_filter_velo(const LAFEM::FilterSequence<SubFilter_>& filter, FilterUpdateCounter* fuc)
      {
        for(auto it = filter.begin(); it != filter.end(); ++it)
        {
          this->_add_filter_velo(it->second, fuc);
        }
      }

      /// adds or updates a global velocity filter
      template<typename LocFilter_, typename Mirror_>
      void _add_filter_velo(const Global::Filter<LocFilter_, Mirror_>& filter, FilterUpdateCounter* fuc)
      {
        this->_add_filter_velo(filter.local(), fuc);
      }

      /// adds or updates a pressure unit filter
      void _add_filter_pres(const UnitFilterTypeP& filter, FilterUpdateCounter* fuc)
      {
        if(fuc != nullptr)
        {
          // get the old unit filter
          UnitFilterTypeP& old_filter = this->_unit_filters_p.at(fuc->at(2u)++);

          // compare layout of old and new filter
          if(!old_filter.get_filter_vector().compare_layout(filter.get_filter_vector()))
            XABORTM("pressure unit filter layout has changed");

          // update the filter
          old_filter.clone(filter, LAFEM::CloneMode::Shallow);
        }
        else
        {
          this->_unit_filters_p.push_back(filter.clone(LAFEM::CloneMode::Shallow));
        }
      }

      /// adds or updates a pressure mean filter
      void _add_filter_pres(const LAFEM::MeanFilter<DT_, IT_>& filter, FilterUpdateCounter* fuc)
      {
        if(fuc != nullptr)
        {
          XASSERTM(this->_mean_filter_p_dual.empty() == filter.get_vec_dual().empty(), "pressure mean filter has changed");
        }
        this->_mean_filter_p_dual.clone(filter.get_vec_dual(), LAFEM::CloneMode::Shallow);
      }

      /// adds or updates a pressure mean filter
      void _add_filter_pres(const Global::MeanFilter<DT_, IT_>& filter, FilterUpdateCounter* fuc)
      {
        if(fuc != nullptr)
        {
          XASSERTM(this->_mean_filter_p_dual.empty() == filter.get_vec_dual().empty(), "pressure mean filter has changed");
        }
        this->_mean_filter_p_dual.clone(filter.get_vec_dual(), LAFEM::CloneMode::Shallow);
      }

      /// adds or updates a pressure none filter (dummy function)
      void _add_filter_pres(const NoneFilterTypeP&, FilterUpdateCounter*)
      {
        // nothing to do here
      }

      /// adds or updates a pressure filter chain
      template<typename First_>
      void _add_filter_pres(const LAFEM::FilterChain<First_>& filter, FilterUpdateCounter* fuc)
      {
        this->_add_filter_pres(filter.first(), fuc);
      }

      /// adds or updates a pressure filter chain
      template<typename First_, typename Second_, typename... Rest_>
      void _add_filter_pres(const LAFEM::FilterChain<First_, Second_, Rest_...>& filter, FilterUpdateCounter* fuc)
      {
        this->_add_filter_pres(filter.first(), fuc);
        this->_add_filter_pres(filter.rest(), fuc);
      }

      /// adds or updates a pressure filter sequence
      template<typename SubFilter_>
      void _add_filter_pres(const LAFEM::FilterSequence<SubFilter_>& filter, FilterUpdateCounter* fuc)
      {
        for(auto it = filter.begin(); it != filter.end(); ++it)
        {
          this->_add_filter_pres(it->second, fuc);
        }
      }

      /// adds or updates a global pressure filter
      template<typename LocFilter_, typename Mirror_>
      void _add_filter_pres(const Global::Filter<LocFilter_, Mirror_>& filter, FilterUpdateCounter* fuc)
      {
        this->_add_filter_pres(filter.local(), fuc);
      }

      /// computes the slip matrix layout Q^T and returns the total number of slip DOFs
      IT_ _compute_slip_matrix_layout()
      {
        const IT_ nv = IT_(_matrix_a.rows());

        _slip_row_ptr.clear();
        _slip_col_idx.clear();
        _slip_fil_dqu.clear();
        _slip_fil_idx.clear();
        _slip_row_ptr.resize(nv+1u, IT_(0));

        if(this->_slip_filters_v.empty())
          return IT_(0);

        // count the total number of slip filter entries
        IT_ count(0u);
        for(std::size_t it(0u); it < this->_slip_filters_v.size(); ++it)
        {
          const IT_ ni = IT_(this->_slip_filters_v.at(it).used_elements());
          const IT_* idx = this->_slip_filters_v.at(it).get_indices();
          for(IT_ i(0); i < ni; ++i)
            ++_slip_row_ptr.at(idx[i]+1);
          count += ni;
        }

        // only empty filters?
        if(count == IT_(0))
          return count;

        // compute row pointer array from aux
        for(IT_ i(0); i < nv; ++i)
        {
          _slip_row_ptr[i+1u] += _slip_row_ptr[i];
        }

        // allocate vectors
        _slip_col_idx.resize(count);
        _slip_fil_dqu.resize(count);
        _slip_fil_idx.resize(count);

        // copy row pointer
        std::vector<IT_> aux(_slip_row_ptr);

        // loop over all slip filters
        IT_ offset(0u);
        for(std::size_t it(0u); it < this->_slip_filters_v.size(); ++it)
        {
          const IT_ ni = IT_(this->_slip_filters_v.at(it).used_elements());
          const IT_* idx = this->_slip_filters_v.at(it).get_indices();

          // and emplace them
          for(IT_ i(0); i < ni; ++i)
          {
            IT_& px = aux[idx[i]];
            _slip_col_idx[px] = offset + i;
            _slip_fil_dqu[px] = IT_(it);
            _slip_fil_idx[px] = i;
            ++px;
          }
          offset += ni;
        }

        XASSERT(count == offset);
        return count;
      }

      /// computes the velocity unit filter mask
      void _compute_unit_mask_v()
      {
        // reset mask to 0
        _unit_mask_v.clear();
        _unit_mask_v.resize(_matrix_a.rows(), 0);

        for(const auto& it : this->_unit_filters_v)
        {
          const IT_ n = IT_(it.used_elements());
          const IT_* idx = it.get_indices();
          for(IT_ k(0); k < n; ++k)
            _unit_mask_v.at(idx[k]) = 1;
        }
      }

      /// computes the pressure unit filter mask
      void _compute_unit_mask_p()
      {
        // reset mask to 0
        _unit_mask_p.clear();
        _unit_mask_p.resize(_matrix_d.rows(), 0);

        for(const auto& it : this->_unit_filters_p)
        {
          const IT_ n = IT_(it.used_elements());
          const IT_* idx = it.get_indices();
          for(IT_ k(0); k < n; ++k)
            _unit_mask_p.at(idx[k]) = 1;
        }
      }
    }; // class DirectStokesCore

    /**
     * \brief Direct Stokes solver class template
     *
     * \see See the documentation of the specializations for details.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class DirectStokesSolver;

    /**
     * \brief Direct Stokes solver class for SaddlePointMatrix systems
     *
     * This specialization implements a direct solver for (Navier-)Stokes systems, which are discretized into a (local)
     * SaddlePointMatrix consisting of matrix blocks A, B and D.
     *
     * \attention
     * Currently, this solver class can only be used if FEAT is compiled and linked with support for the UMFPACK library.
     *
     * \author Peter Zajac
     */
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename FilterV_, typename FilterP_>
    class DirectStokesSolver<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>, LAFEM::TupleFilter<FilterV_, FilterP_>> :
      public Solver::SolverBase<LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL>>
    {
    public:
      /// our system matrix type
      typedef LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_> SystemMatrixType;
      /// our system vector type
      typedef LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL> SystemVectorType;
      /// our system filter type
      typedef LAFEM::TupleFilter<FilterV_, FilterP_> SystemFilterType;

      /// our base class
      typedef Solver::SolverBase<SystemVectorType> BaseClass;

      /// the direct Stokes core class type
      typedef DirectStokesCore<double, Index, MatrixA_, MatrixB_, MatrixD_> StokesCoreType;

      /// the solver matrix type
      typedef typename StokesCoreType::SolverMatrixType SolverMatrixType;
      /// the solver vector type
      typedef typename StokesCoreType::SolverVectorType SolverVectorType;
      /// the solver filter type; this is always a NoneFilter
      typedef typename StokesCoreType::SolverFilterType SolverFilterType;

    protected:
      /// the system matrix
      const SystemMatrixType& _matrix_sys;
      /// the system filter
      const SystemFilterType& _filter_sys;

      /// the direct stokes core
      StokesCoreType _stokes_core;

      /// our RHS/SOL vectors
      SolverVectorType _vec_sol, _vec_rhs;

      /// the actual direct solver; opaque SolverBase pointer
      std::shared_ptr<Solver::SolverBase<SolverVectorType>> _direct_solver;

#if defined(FEAT_HAVE_CUDSS) || defined(DOXYGEN)
      /// our CUDSS solver object; used if preferred backend is PreferredBackend::cuda
      std::shared_ptr<Solver::CUDSS> _cudss_solver;
#endif // FEAT_HAVE_CUDSS
#if defined(FEAT_HAVE_MKL) || defined(DOXYGEN)
      /// our MKLDSS solver object; used if preferred backend is PreferredBackend::mkl
      std::shared_ptr<Solver::MKLDSS> _mkldss_solver;
#endif // FEAT_HAVE_MKL
#if defined(FEAT_HAVE_UMFPACK) || defined(DOXYGEN)
      /// our UMFPACK solver object; used if no other direct solver was chosen
      std::shared_ptr<Solver::Umfpack> _umfpack_solver;
#endif // FEAT_HAVE_UMFPACK

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix_sys
       * A \resident reference to the system matrix object. The matrix does not need to be initialized yet.
       *
       * \param[in] filter_sys
       * A \resident reference to the system filter object. The filter does not need to be initialized yet.
       */
      explicit DirectStokesSolver(const SystemMatrixType& matrix_sys, const SystemFilterType& filter_sys) :
        _matrix_sys(matrix_sys),
        _filter_sys(filter_sys),
        _stokes_core(_matrix_sys.block_a(), _matrix_sys.block_b(), _matrix_sys.block_d()),
        _direct_solver()
      {
#ifdef FEAT_HAVE_CUDSS
        if(Backend::get_preferred_backend() == PreferredBackend::cuda)
        {
          // create a CUDSS solver
          _cudss_solver = Solver::new_cudss(_stokes_core.get_solver_matrix());
          _direct_solver = _cudss_solver;
          // main diagonal is not absolutely necessary for CUDSS, but it seems to improve precision a bit
          _stokes_core.set_need_main_diagonal(true);
        }
        else
#endif // FEAT_HAVE_CUDSS
#ifdef FEAT_HAVE_MKL
        if(Backend::get_preferred_backend() == PreferredBackend::mkl)
        {
          // create a MKLDSS solver
          _mkldss_solver = Solver::new_mkl_dss(_stokes_core.get_solver_matrix());
          _direct_solver = _mkldss_solver;
          // MKL emphasized that main diagonal is absolutely mandatory
          _stokes_core.set_need_main_diagonal(true);
        }
        else
#endif // FEAT_HAVE_MKL
#ifdef FEAT_HAVE_UMFPACK
        // fallback or generic: try to create UMFPACK solver
        {
          // create an UMFPACK solver
          _umfpack_solver = Solver::new_umfpack(_stokes_core.get_solver_matrix());
          _direct_solver = _umfpack_solver;
        }
#else
        {
          XABORTM("DirectStokesSolver can only be used if UMFPACK is available");
        }
#endif
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
#ifdef FEAT_HAVE_CUDSS
        if(_cudss_solver)
          return "DirectStokesSolver<SaddlePointMatrix>[CUDSS]";
#endif
#ifdef FEAT_HAVE_MKL
        if(_mkldss_solver)
          return "DirectStokesSolver<SaddlePointMatrix>[MKLDSS]";
#endif
#ifdef FEAT_HAVE_UMFPACK
        if(_umfpack_solver)
          return "DirectStokesSolver<SaddlePointMatrix>[UMFPACK]";
#endif
        return "DirectStokesSolver<SaddlePointMatrix>[???]";
      }

      /// Performs the symbolic initialization
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // add the filter
        _stokes_core.set_filters(_filter_sys);

        // initialize
        _stokes_core.init_symbolic();
        if(_direct_solver)
          _direct_solver->init_symbolic();

        // create vectors
        _vec_sol = _stokes_core.create_solver_vector();
        _vec_rhs = _stokes_core.create_solver_vector();
      }

      /// Performs the numeric initialization
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // update filter
        _stokes_core.update_filters(_filter_sys);

        // initialize
        _stokes_core.init_numeric();
        if(_direct_solver)
          _direct_solver->init_numeric();
      }

      /// Releases all data allocated during numeric initialization
      virtual void done_numeric() override
      {
        if(_direct_solver)
          _direct_solver->done_numeric();
        _stokes_core.done_numeric();
        BaseClass::done_numeric();
      }

      /// Releases all data allocated during symbolic initialization
      virtual void done_symbolic() override
      {
        _vec_rhs.clear();
        _vec_sol.clear();
        if(_direct_solver)
          _direct_solver->done_symbolic();
        _stokes_core.done_symbolic();
        BaseClass::done_symbolic();
      }

      /**
       * \brief Applies the direct Stokes solver
       *
       * \param[out] vec_cor
       * The vector that shall receive the solution of the linear system. It is assumed to be allocated, but
       * its numerical contents may be undefined upon calling this method.
       *
       * \param[in] vec_def
       * The vector that represents the right-hand-side of the linear system to be solved.
       *
       * \attention vec_cor and vec_def must \b not refer to the same vector object!
       *
       * \returns
       * A Status code that represents the status of the solution step.
       */
      virtual Solver::Status apply(SystemVectorType& vec_cor, const SystemVectorType& vec_def) override
      {
        // upload RHS vector
        _stokes_core.upload_vector(this->_vec_rhs, vec_def.template at<0>(), vec_def.template at<1>());

        // apply the actual solver
        Solver::Status status = _direct_solver->apply(this->_vec_sol, this->_vec_rhs);

        // download solution vector
        _stokes_core.download_vector(this->_vec_sol, vec_cor.template at<0>(), vec_cor.template at<1>());

        // return status
        return status;
      }
    }; // class DirectStokesSolver<LAFEM::SaddlePointMatrix<...>, LAFEM::TupleFilter<...>>

    /**
     * \brief Direct Stokes solver class for SaddlePointMatrix systems
     *
     * \attention
     * Currently, this solver class is merely a wrapper around a local direct solver, i.e. this class can only be used
     * to solve systems on a single process, e.g. in the case of a single-process coarse grid solver.
     * An "actual" parallel direct solver may or may not be implemented in future.
     *
     * \attention
     * Currently, this solver class can only be used if FEAT is compiled and linked with support for the UMFPACK library.
     *
     * \author Peter Zajac
     */
    template<typename LocalMatrix_, typename LocalFilter_, typename Mirror_>
    class DirectStokesSolver<Global::Matrix<LocalMatrix_, Mirror_, Mirror_>, Global::Filter<LocalFilter_, Mirror_>> :
      public Solver::SolverBase<Global::Vector<typename LocalMatrix_::VectorTypeR, Mirror_>>
    {
    public:
      /// our system matrix type
      typedef Global::Matrix<LocalMatrix_, Mirror_, Mirror_> SystemMatrixType;
      /// our system vector type
      typedef Global::Vector<typename LocalMatrix_::VectorTypeR, Mirror_> SystemVectorType;
      /// our system filter type
      typedef Global::Filter<LocalFilter_, Mirror_> SystemFilterType;

      /// our base class
      typedef Solver::SolverBase<SystemVectorType> BaseClass;

    protected:
      /// the system matrix
      const SystemMatrixType& _matrix_sys;
      /// the system filter
      const SystemFilterType& _filter_sys;

      /// our local solver
      DirectStokesSolver<LocalMatrix_, LocalFilter_> _local_solver;

    public:
      explicit DirectStokesSolver(const SystemMatrixType& matrix_sys, const SystemFilterType& filter_sys) :
        _matrix_sys(matrix_sys),
        _filter_sys(filter_sys),
        _local_solver(_matrix_sys.local(), _filter_sys.local())
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "DirectStokesSolver<Global::Matrix>";
      }

      /// Performs the symbolic initialization
      virtual void init_symbolic() override
      {
        // ensure that we are on one process only...
        const Dist::Comm* comm = _matrix_sys.get_comm();
        if(comm != nullptr)
        {
          XASSERTM(comm->size() <= 1, "DirectStokesSolver can only work on one process");
        }

        BaseClass::init_symbolic();
        _local_solver.init_symbolic();
      }

      /// Performs the numeric initialization
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        _local_solver.init_numeric();
      }

      /// Releases all data allocated during numeric initialization
      virtual void done_numeric() override
      {
        _local_solver.done_numeric();
        BaseClass::done_numeric();
      }

      /// Releases all data allocated during symbolic initialization
      virtual void done_symbolic() override
      {
        _local_solver.done_symbolic();
        BaseClass::done_symbolic();
      }

      /**
       * \brief Applies the direct Stokes solver
       *
       * \param[out] vec_cor
       * The vector that shall receive the solution of the linear system. It is assumed to be allocated, but
       * its numerical contents may be undefined upon calling this method.
       *
       * \param[in] vec_def
       * The vector that represents the right-hand-side of the linear system to be solved.
       *
       * \attention vec_cor and vec_def must \b not refer to the same vector object!
       *
       * \returns
       * A Status code that represents the status of the solution step.
       */
      virtual Solver::Status apply(SystemVectorType& vec_sol, const SystemVectorType& vec_rhs) override
      {
        return _local_solver.apply(vec_sol.local(), vec_rhs.local());
      }
    }; // class DirectStokesSolver<Global::Matrix<...>, Global::Filter<...>>

    /**
     * \brief Creates a new DirectStokesSolver object
     *
     * \param[in] matrix
     * A \resident reference to the system matrix.
     *
     * \param[in] filter
     * A \resident reference to the system filter.
     *
     * \returns
     * A shared pointer to a new DirectStokesSolver object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<DirectStokesSolver<Matrix_, Filter_>> new_direct_stokes_solver(const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<DirectStokesSolver<Matrix_, Filter_>>(matrix, filter);
    }
  } // namespace Solver
} // namespace FEAT
