// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/global/alg_dof_parti_system.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/pmdcdsc_matrix.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Base-Class for solvers based on Algebraic-DOF-Partitioning
     *
     * This class template specifies a common interface for solvers based on algebraic DOF partitioning.
     * There are two specializations of this class, which depend on the type of the system matrix.
     *
     * \attention
     * This class is only defined for this doxygen documentation! There exists no generic implementation
     * of this class template, but only the three specializations for Global::Matrix, Global::PMDCDSCMatrix
     * and any purely local LAFEM container, however, all of these specializations implement the same protected
     * interface, which is documented in this generic template.
     *
     * \tparam Matrix_
     * The type of the system matrix that is to be used by the solver.
     *
     * \tparam Filter_
     * The type of the system filter.
     *
     * \tparam SolverBase_
     * Specifies the base-class that this class should derive from. Defaults to Solver::SolverBase, but
     * it may also be set to Solver::IterativeSolver if the derived class represents an iterative solver.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_, typename SolverBase_ = Solver::SolverBase<typename Matrix_::VectorTypeL>>
#ifndef DOXYGEN
    class ADPSolverBase;
#else // DOXYGEN
    class ADPSolverBase :
      public SolverBase_
    {
    protected:
      /// constructor
      explicit ADPSolverBase(const GlobalMatrixType& matrix, const GlobalFilterType& filter);

      /// \returns A pointer to the underlying communicator
      const Dist::Comm* _get_comm() const;

      /// \returns The size of the local ADP vector; equal to _get_num_owned_dofs()
      Index _get_adp_vector_size() const;

      /// \returns The number of rows of the local ADP matrix; equal to _get_num_owned_dofs()
      Index _get_adp_matrix_num_rows() const;

      /// \returns The number of columns of the local ADP matrix; equal to _get_num_global_dofs()
      Index _get_adp_matrix_num_cols() const;

      /// \returns The number of nonzero entries of the local ADP matrix
      Index _get_adp_matrix_num_nzes() const;

      /// \returns The total number of nonzero entries of the global ADP matrix
      Index _get_num_global_nonzeros() const;

      /// \returns The total number of global DOFs on all processes
      Index _get_num_global_dofs() const;

      /// \returns The number of global DOFs owned by this process
      Index _get_num_owned_dofs() const;

      /// \returns The index of the first global DOF owned by this process
      Index _get_global_dof_offset() const;

      /// \returns The block information of the algebraic dof partitioning as an XML string.
      String _get_adp_block_information() const;

      /**
       * \brief Uploads the ADP matrix structure to the given arrays
       *
       * \param[out] row_ptr
       * A \transient pointer to an array of length _get_adp_matrix_num_rows()+1 that receives the
       * the row-pointer array of the local ADP matrix.
       *
       * \param[out] col_idx
       * A \transient pointer to an array of length _get_adp_matrix_num_nzes() that receives the
       * column index array of the local ADP matrix.
       */
      template<typename RPT_, typename CIT_>
      void _upload_symbolic(RPT_* row_ptr, CIT_* col_idx);

      /**
       * \brief Uploads the (filtered) ADP matrix values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_matrix_num_nzes() that receives the
       * filtered matrix values of the local ADP matrix.
       *
       * \param[in] row_ptr
       * A \transient pointer to the row-pointer array of the ADP matrix as returned by the
       * _upload_symbolic() function.
       *
       * \param[in] col_idx
       * A \transient pointer to the column-index array of the ADP matrix as returned by the
       * _upload_symbolic() function.
       */
      template<typename DTV_, typename RPT_, typename CIT_>
      void _upload_numeric(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx);

      /**
       * \brief Uploads the ADP vector values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_vector_size() that receives the
       * local ADP values of the input vector given by \p vector.
       *
       * \param[in] vector
       * A \transient reference to the local vector that is to be uploaded to the ADP value array.
       */
      template<typename DTV_>
      void _upload_vector(DTV_* val, const LocalVectorType& vector);

      /**
       * \brief Downloads the ADP vector values from the given array
       *
       * \param[inout] vector
       * A \transient reference to the local vector that receives the downloaded ADP value array.
       * It is assumed to be allocated to the correct size, but its contents on entry are ignored.
       *
       * \param[in] val
       * A \transient pointer to an array of length _get_adp_vector_size() that contains the
       * ADP values that are to be downloaded to the output vector given by \p vector.
       */
      template<typename DTV_>
      void _download_vector(LocalVectorType& vector, const DTV_* val);
    }; // class ADPSolverBase
#endif // DOXYGEN

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * \brief Implementation of ADPSolverBase interface for Global::Matrix instances
     *
     * This class implements the ADPSolverBase interface for Global::Matrix instances, which derive
     * from standard finite element discretizations. Most of the dirty work is outsourced
     * to the Global::AlgDofParti and Global::AlgDofPartiSystem classes.
     *
     * \author Peter Zajac
     */
    template<typename LocalMatrix_, typename Mirror_, typename LocalFilter_, typename SolverBase_>
    class ADPSolverBase<Global::Matrix<LocalMatrix_, Mirror_, Mirror_>, Global::Filter<LocalFilter_, Mirror_>, SolverBase_> :
      public SolverBase_
    {
    public:
      /// the local matrix type
      typedef LocalMatrix_ LocalMatrixType;
      /// the local vector type
      typedef typename LocalMatrixType::VectorTypeL LocalVectorType;
      /// the local filter type
      typedef LocalFilter_ LocalFilterType;
      /// the mirror type
      typedef Mirror_ MirrorType;

      /// the global matrix type
      typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;
      /// the global vector type
      typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;
      /// the global filter type
      typedef Global::Filter<LocalFilterType, MirrorType> GlobalFilterType;

      /// our base class
      typedef SolverBase_ BaseClass;

      /// the system matrix type; coincides with the global matrix type
      typedef GlobalMatrixType MatrixType;
      /// the system vector type; coincides with the global vector type
      typedef GlobalVectorType VectorType;
      /// the system filter type; coincides with the global filter type
      typedef GlobalFilterType FilterType;

      /// the algebraic dof partitioning instance for this class
      typedef Global::AlgDofParti<LocalVectorType, MirrorType> AlgDofPartiType;
      /// the algebraic dof partitioning system instance for this class
      typedef Global::AlgDofPartiSystem<LocalMatrixType, LocalFilterType, MirrorType> AlgDofPartiSystemType;

      /// our internal data type
      typedef typename AlgDofPartiType::DataType DataType;
      /// our internal index type
      typedef typename AlgDofPartiType::IndexType IndexType;

      /// the ADP matrix type; this is always a (globally partitioned) CSR matrix
      typedef LAFEM::SparseMatrixCSR<DataType, IndexType> ADPMatrixType;
      /// the ADP vector type; this is always a (globally partitioned) dense vector
      typedef LAFEM::DenseVector<DataType, IndexType> ADPVectorType;

    protected:
      /// the global system matrix
      const GlobalMatrixType& _system_matrix;
      /// the global system filter
      const GlobalFilterType& _system_filter;
      /// the algebraic DOF partitioning
      std::shared_ptr<AlgDofPartiSystemType> _adp_system;

      /**
       * \brief protected constructor
       *
       * \param[in] matrix
       * A \resident reference to the system matrix
       *
       * \param[in] filter
       * A \resident reference to the system filter
       */
      explicit ADPSolverBase(const GlobalMatrixType& matrix, const GlobalFilterType& filter) :
        BaseClass(),
        _system_matrix(matrix),
        _system_filter(filter),
        _adp_system()
      {
      }

    public:
      /// \returns A pointer to the underlying communicator
      const Dist::Comm* _get_comm() const
      {
        return this->_system_matrix.get_comm();
      }

      /**
       * \brief Symbolic initialization
       *
       * This function performs the symbolic initialization of the algebraic DOF partitioning, which pre-computes
       * the actual algebraic DOF partitioning. After this function call, the derived class has to call the
       * corresponding upload functions to retrieve the ADP matrix arrays.
       */
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // assemble the algebraic dof partitioning
        this->_adp_system = std::make_shared<AlgDofPartiSystemType>(this->_system_matrix, this->_system_filter);
        this->_adp_system->init_symbolic();
      }

      /**
       * \brief Symbolic finalization
       *
       * This function releases all data allocated during the symbolic initialization phase.
       */
      virtual void done_symbolic() override
      {
        _adp_system.reset();

        BaseClass::done_symbolic();
      }

    protected:
      /// \returns The size of the local ADP vector; equal to _get_num_owned_dofs()
      Index _get_adp_vector_size() const
      {
        return this->_adp_system->get_adp_vector_size();
      }

      /// \returns The number of rows of the local ADP matrix; equal to _get_num_owned_dofs()
      Index _get_adp_matrix_num_rows() const
      {
        return this->_adp_system->get_adp_matrix_rows();
      }

      /// \returns The number of columns of the local ADP matrix; equal to _get_num_global_dofs()
      Index _get_adp_matrix_num_cols() const
      {
        return this->_adp_system->get_adp_matrix_cols();
      }

      /// \returns The number of nonzero entries of the local ADP matrix
      Index _get_adp_matrix_num_nzes() const
      {
        return this->_adp_system->get_adp_matrix_nzes();
      }

      /// \returns The total number of nonzero entries of the global ADP matrix
      Index _get_num_global_nonzeros() const
      {
        return this->_adp_system->get_num_global_nonzeros();
      }

      /// \returns The total number of global DOFs on all processes
      Index _get_num_global_dofs() const
      {
        return this->_adp_system->get_num_global_dofs();
      }

      /// \returns The number of global DOFs owned by this process
      Index _get_num_owned_dofs() const
      {
        return this->_adp_system->get_num_owned_dofs();
      }

      /// \returns The index of the first global DOF owned by this process
      Index _get_global_dof_offset() const
      {
        return this->_adp_system->get_global_dof_offset();
      }

      /// \returns The block information of the algebraic dof partitioning as an XML string.
      String _get_adp_block_information() const
      {
        return this->_adp_system->get_alg_dof_parti()->get_block_information();
      }

      /**
       * \brief Uploads the ADP matrix structure to the given arrays
       *
       * \param[out] row_ptr
       * A \transient pointer to an array of length _get_adp_matrix_num_rows()+1 that receives the
       * the row-pointer array of the ADP matrix.
       *
       * \param[out] col_idx
       * A \transient pointer to an array of length _get_adp_matrix_num_nzes() that receives the
       * column index array of the ADP matrix.
       */
      template<typename RPT_, typename CIT_>
      void _upload_symbolic(RPT_* row_ptr, CIT_* col_idx)
      {
        this->_adp_system->upload_matrix_symbolic(row_ptr, col_idx);
      }

      /**
       * \brief Uploads the (filtered) ADP matrix values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_matrix_num_nzes() that receives the
       * filtered matrix values of the local ADP matrix.
       *
       * \param[in] row_ptr
       * A \transient pointer to the row-pointer array of the ADP matrix as returned by the
       * _upload_symbolic() function.
       *
       * \param[in] col_idx
       * A \transient pointer to the column-index array of the ADP matrix as returned by the
       * _upload_symbolic() function.
       */
      template<typename DTV_, typename RPT_, typename CIT_>
      void _upload_numeric(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx)
      {
        this->_adp_system->upload_matrix_numeric(val, row_ptr, col_idx);
        this->_adp_system->upload_filter();
        this->_adp_system->filter_matrix(val, row_ptr, col_idx);
      }

      /**
       * \brief Uploads the ADP vector values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_vector_size() that receives the
       * local ADP values of the input vector given by \p vector.
       *
       * \param[in] vector
       * A \transient reference to the global vector that is to be uploaded to the ADP value array.
       */
      template<typename DTV_>
      void _upload_vector(DTV_* val, const GlobalVectorType& vector)
      {
        this->_adp_system->upload_vector(val, vector.local());
      }

      /**
       * \brief Uploads the ADP vector values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_vector_size() that receives the
       * local ADP values of the input vector given by \p vector.
       *
       * \param[in] vector
       * A \transient reference to the local vector that is to be uploaded to the ADP value array.
       */
      template<typename DTV_>
      void _upload_vector(DTV_* val, const LocalVectorType& vector)
      {
        this->_adp_system->upload_vector(val, vector);
      }

      /**
       * \brief Downloads the ADP vector values from the given array
       *
       * \param[inout] vector
       * A \transient reference to the global vector that receives the downloaded ADP value array.
       * It is assumed to be allocated to the correct size, but its contents on entry are ignored.
       *
       * \param[in] val
       * A \transient pointer to an array of length _get_adp_vector_size() that contains the
       * ADP values that are to be downloaded to the output vector given by \p vector.
       */
      template<typename DTV_>
      void _download_vector(GlobalVectorType& vector, const DTV_* val)
      {
        this->_adp_system->download_vector(vector.local(), val);
      }

      /**
       * \brief Downloads the ADP vector values from the given array
       *
       * \param[inout] vector
       * A \transient reference to the local vector that receives the downloaded ADP value array.
       * It is assumed to be allocated to the correct size, but its contents on entry are ignored.
       *
       * \param[in] val
       * A \transient pointer to an array of length _get_adp_vector_size() that contains the
       * ADP values that are to be downloaded to the output vector given by \p vector.
       */
      template<typename DTV_>
      void _download_vector(LocalVectorType& vector, const DTV_* val)
      {
        this->_adp_system->download_vector(vector, val);
      }
    }; // class ADPSolverBase<Global::Matrix<...>, ...>

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * \brief Implementation of ADPSolverBase interface for Global::PMDCDSCMatrix instances
     *
     * This class implements the ADPSolverBase interface for Global::PMDCDSCMatrix instances,
     * which are used to implement the global Schur-complement matrices appearing in Navier-Stokes
     * equations with discontinuous pressure spaces. The most prominent use of this specialization
     * is the Pressure-Poisson-Problem solver of the infamous "PP" application.
     *
     * Most of the dirty work is outsourced to the Global::PMDCDSCMatrix class, most notably to
     * its three member functions Global::PMDCDSCMatrix::adp_compute_counts(),
     * Global::PMDCDSCMatrix::adp_upload_symbolic() and Global::PMDCDSCMatrix::adp_upload_numeric().
     *
     * \author Peter Zajac
     */
    template<typename MatrixB_, typename MatrixD_, typename GlobalFilter_, typename SolverBase_>
    class ADPSolverBase<Global::PMDCDSCMatrix<MatrixB_, MatrixD_>, GlobalFilter_, SolverBase_> :
      public SolverBase_
    {
    public:
      /// the global matrix type
      typedef Global::PMDCDSCMatrix<MatrixB_, MatrixD_> GlobalMatrixType;
      /// the global vector type
      typedef typename GlobalMatrixType::VectorTypeL GlobalVectorType;
      /// the global filter type
      typedef GlobalFilter_ GlobalFilterType;

      /// the local Schur-complement matrix type
      typedef typename GlobalMatrixType::LocalMatrixTypeS LocalMatrixTypeS;
      /// our internal local vector type; always a DenseVector<DataType, IndexType>
      typedef typename LocalMatrixTypeS::VectorTypeL LocalVectorType;

      /// our base class
      typedef SolverBase_ BaseClass;

      /// the system vector type
      typedef GlobalVectorType VectorType;

      /// our internal data type
      typedef typename LocalMatrixTypeS::DataType DataType;
      /// our internal index type
      typedef typename LocalMatrixTypeS::IndexType IndexType;

    protected:
      /// the global system matrix
      const GlobalMatrixType& _system_matrix;
      /// the global system filter
      const GlobalFilterType& _system_filter;

    private:
      /// the total global DOF count
      Index _global_dof_count;
      /// our global DOF offset
      Index _global_dof_offset;
      /// the owned DOF count
      Index _owned_dof_count;
      /// the number of owned non-zero entries
      Index _owned_num_nzes;
      /// the number of global non-zero entries
      Index _global_num_nzes;

    protected:
      /**
       * \brief protected constructor
       *
       * \param[in] matrix
       * A \resident reference to the system matrix
       *
       * \param[in] filter
       * A \resident reference to the system filter
       */
      explicit ADPSolverBase(const GlobalMatrixType& matrix, const GlobalFilterType& filter) :
        BaseClass(),
        _system_matrix(matrix),
        _system_filter(filter),
        _global_dof_count(0),
        _global_dof_offset(0),
        _owned_dof_count(0),
        _owned_num_nzes(0),
        _global_num_nzes(0)
      {
      }

    public:
      /// \returns A pointer to the underlying communicator
      const Dist::Comm* _get_comm() const
      {
        return this->_system_matrix.get_comm();
      }

      /**
       * \brief Symbolic initialization
       *
       * This function performs the symbolic initialization of the algebraic DOF partitioning, which boils down to
       * computing the ADP matrix dimensions for this specialization.
       */
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // compute the ADP dimensions
        _system_matrix.adp_compute_counts(_global_dof_offset, _global_dof_count, _owned_dof_count, _owned_num_nzes, _global_num_nzes);
      }

    protected:
      /// \returns The size of the local ADP vector; equal to _get_num_owned_dofs()
      Index _get_adp_vector_size() const
      {
        return _get_num_owned_dofs();
      }

      /// \returns The number of rows of the local ADP matrix; equal to _get_num_owned_dofs()
      Index _get_adp_matrix_num_rows() const
      {
        return _get_num_owned_dofs();
      }

      /// \returns The number of columns of the local ADP matrix; equal to _get_num_global_dofs()
      Index _get_adp_matrix_num_cols() const
      {
        return _get_num_global_dofs();
      }

      /// \returns The number of nonzero entries of the local ADP matrix
      Index _get_adp_matrix_num_nzes() const
      {
        return _owned_num_nzes;
      }

      /// \returns The total number of nonzero entries of the global ADP matrix
      Index _get_num_global_nonzeros() const
      {
        return _global_num_nzes;
      }

      /// \returns The total number of global DOFs on all processes
      Index _get_num_global_dofs() const
      {
        return _global_dof_count;
      }

      /// \returns The number of global DOFs owned by this process
      Index _get_num_owned_dofs() const
      {
        return _owned_dof_count;
      }

      /// \returns The index of the first global DOF owned by this process
      Index _get_global_dof_offset() const
      {
        return _global_dof_offset;
      }

      /// \returns The block information of the algebraic dof partitioning as an XML string.
      String _get_adp_block_information() const
      {
        XABORTM("Block information for ADPSolverBase<PMDCDSCMatrix> not available yet!");
        return "-N/A-";
      }

      /**
       * \brief Uploads the ADP matrix structure to the given arrays
       *
       * \param[out] row_ptr
       * A \transient pointer to an array of length _get_adp_matrix_num_rows()+1 that receives the
       * the row-pointer array of the ADP matrix.
       *
       * \param[out] col_idx
       * A \transient pointer to an array of length _get_adp_matrix_num_nzes() that receives the
       * column index array of the ADP matrix.
       */
      template<typename RPT_, typename CIT_>
      void _upload_symbolic(RPT_* row_ptr, CIT_* col_idx)
      {
        _system_matrix.adp_upload_symbolic(row_ptr, col_idx, _global_dof_offset);
      }

      /**
       * \brief Uploads the (filtered) ADP matrix values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_matrix_num_nzes() that receives the
       * filtered matrix values of the local ADP matrix.
       *
       * \param[in] row_ptr
       * A \transient pointer to the row-pointer array of the ADP matrix as returned by the
       * _upload_symbolic() function.
       *
       * \param[in] col_idx
       * A \transient pointer to the column-index array of the ADP matrix as returned by the
       * _upload_symbolic() function.
       */
      template<typename DTV_, typename RPT_, typename CIT_>
      void _upload_numeric(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx)
      {
        _system_matrix.adp_upload_numeric(val, row_ptr, col_idx);
        this->_filter_mat(val, row_ptr, col_idx, this->_system_filter.local());
      }

      /**
       * \brief Uploads the ADP vector values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_vector_size() that receives the
       * local ADP values of the input vector given by \p vector.
       *
       * \param[in] vector
       * A \transient reference to the local vector that is to be uploaded to the ADP value array.
       */
      template<typename DTV_>
      void _upload_vector(DTV_* val, const LAFEM::DenseVector<DataType, IndexType>& vector)
      {
        const Index n = _get_adp_vector_size();
        const DataType* vx = vector.elements();

        XASSERT(vector.size() == n);

        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < n; ++i)
          val[i] = DTV_(vx[i]);
      }

      /**
       * \brief Downloads the ADP vector values from the given array
       *
       * \param[inout] vector
       * A \transient reference to the local vector that receives the downloaded ADP value array.
       * It is assumed to be allocated to the correct size, but its contents on entry are ignored.
       *
       * \param[in] val
       * A \transient pointer to an array of length _get_adp_vector_size() that contains the
       * ADP values that are to be downloaded to the output vector given by \p vector.
       */
      template<typename DTV_>
      void _download_vector(LAFEM::DenseVector<DataType, IndexType>& vector, const DTV_* val)
      {
        const Index n = _get_adp_vector_size();
        DataType* vx = vector.elements();

        XASSERT(vector.size() == n);

        FEAT_PRAGMA_OMP(parallel for)
          for(Index i = 0; i < n; ++i)
            vx[i] = DataType(val[i]);
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_>
      void _filter_mat(DTV_*, const RPT_*, const CIT_*,
        const LAFEM::NoneFilter<DataType, IndexType>&)
      {
        // nothing to do here
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_>
      void _filter_mat(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx,
        const LAFEM::UnitFilter<DataType, IndexType>& filter)
      {
        const IndexType n = filter.used_elements();
        const IndexType* fil_idx = filter.get_indices();
        for(IndexType i = 0; i < n; ++i)
        {
          const Index row = fil_idx[i];
          for(RPT_ j = row_ptr[row]; j < row_ptr[row+1]; ++j)
            val[j] = DTV_(Index(col_idx[j]) == row ? 1 : 0);
        }
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_>
      void _filter_mat(DTV_*, const RPT_*, const CIT_*,
        const LAFEM::MeanFilter<DataType, IndexType>& filter)
      {
        XASSERTM(filter.empty(), "LAFEM::MeanFilter is not supported by ADPSolverBase yet!");
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_>
      void _filter_mat(DTV_*, const RPT_*, const CIT_*,
        const Global::MeanFilter<DataType, IndexType>& filter)
      {
        XASSERTM(filter.empty(), "LAFEM::MeanFilter is not supported by ADPSolverBase yet!");
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_, typename SubFilter_>
      void _filter_mat(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx,
        const LAFEM::FilterSequence<SubFilter_>& filter)
      {
        for(const auto& sf : filter)
          this->_filter_mat(val, row_ptr, col_idx, sf.second);
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_, typename First_, typename... Rest_>
      void _filter_mat(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx,
        const LAFEM::FilterChain<First_, Rest_...>& filter)
      {
        this->_filter_mat(val, row_ptr, col_idx, filter.first());
        this->_filter_mat(val, row_ptr, col_idx, filter.rest());
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_, typename First_>
      void _filter_mat(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx,
        const LAFEM::FilterChain<First_>& filter)
      {
        this->_filter_mat(val, row_ptr, col_idx, filter.first());
      }
    }; // class ADPSolverBase<Global::PMDCDSCMatrix<...>, ...>

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /**
     * \brief Implementation of ADPSolverBase interface (purely local) LAFEM matrices
     *
     * This class implements the ADPSolverBase interface for all possible instances of purely local LAFEM matrices,
     * including meta matrix containers such as e.g. LAFEM::TupleMatrix or LAFEM::SaddlePointMatrix.
     * In the purely local case, the matrices are not partitioned and in the current implementation, this class does
     * not perform any partitioning of the matrix and the derived solver class operates only on a single MPI process,
     * which is represented by the self-communicator.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_, typename SolverBase_>
    class ADPSolverBase :
      public SolverBase_
    {
      /// this specialization expects a local matrix
      static_assert(Matrix_::is_local, "invalid instantiation of ADPSolverBase for non-local matrix type!");

    public:
      /// the (local) matrix type
      typedef Matrix_ MatrixType;
      /// the (local) vector type
      typedef typename MatrixType::VectorTypeL VectorType;
      /// the (local) filter type
      typedef Filter_ FilterType;

      /// our base class
      typedef SolverBase_ BaseClass;

      /// our internal data type
      typedef typename MatrixType::DataType DataType;
      /// our internal index type
      typedef typename MatrixType::IndexType IndexType;

    protected:
      /**
       * \brief self-communicator object
       *
       * Many third-party libraries always expect a MPI_Comm object even if we operate only on a single
       * MPI process, therefore this class creates a custom self-communicator that can be passed on to
       * the third-party libraries.
       */
      Dist::Comm _comm_self;
      /// the system matrix
      const MatrixType& _system_matrix;
      /// the system filter
      const FilterType& _system_filter;

    protected:
      /**
       * \brief protected constructor
       *
       * \param[in] matrix
       * A \resident reference to the system matrix
       *
       * \param[in] filter
       * A \resident reference to the system filter
       */
      explicit ADPSolverBase(const MatrixType& matrix, const FilterType& filter) :
        _comm_self(Dist::Comm::self()),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

      /// \returns A pointer to the underlying (self) communicator
      const Dist::Comm* _get_comm() const
      {
        return &this->_comm_self;
      }

      /// \returns The size of the local ADP vector; equal to _get_num_owned_dofs()
      Index _get_adp_vector_size() const
      {
        return _system_matrix.rows();
      }

      /// \returns The number of rows of the local ADP matrix; equal to _get_num_owned_dofs()
      Index _get_adp_matrix_num_rows() const
      {
        return _system_matrix.rows();
      }

      /// \returns The number of columns of the local ADP matrix; equal to _get_num_global_dofs()
      Index _get_adp_matrix_num_cols() const
      {
        return _system_matrix.columns();
      }

      /// \returns The number of nonzero entries of the local ADP matrix
      Index _get_adp_matrix_num_nzes() const
      {
        return _system_matrix.template used_elements<LAFEM::Perspective::pod>();
      }

      /// \returns The total number of nonzero entries of the global ADP matrix
      Index _get_num_global_nonzeros() const
      {
        return _system_matrix.template used_elements<LAFEM::Perspective::pod>();
      }

      /// \returns The total number of global DOFs on all processes
      Index _get_num_global_dofs() const
      {
        return _system_matrix.rows();
      }

      /// \returns The number of global DOFs owned by this process
      Index _get_num_owned_dofs() const
      {
        return _system_matrix.rows();
      }

      /// \returns The index of the first global DOF owned by this process
      Index _get_global_dof_offset() const
      {
        return Index(0);
      }

      /// \returns The block information of the algebraic dof partitioning as an XML string.
      String _get_adp_block_information() const
      {
        XABORTM("Block information for ADPSolverBase<...> not available yet!");
        return "-N/A-";
      }

      /**
       * \brief Uploads the ADP matrix structure to the given arrays
       *
       * \param[out] row_ptr
       * A \transient pointer to an array of length _get_adp_matrix_num_rows()+1 that receives the
       * the row-pointer array of the ADP matrix.
       *
       * \param[out] col_idx
       * A \transient pointer to an array of length _get_adp_matrix_num_nzes() that receives the
       * column index array of the ADP matrix.
       */
      template<typename RPT_, typename CIT_>
      void _upload_symbolic(RPT_* row_ptr, CIT_* col_idx)
      {
        const Index num_rows = _system_matrix.rows();
        row_ptr[0] = RPT_(0);
        for(Index i(0); i < num_rows; ++i)
        {
          row_ptr[i+1u] = row_ptr[i] + RPT_(_system_matrix.get_row_col_indices(i, &col_idx[row_ptr[i]], CIT_(0)));
        }
      }

      /**
       * \brief Uploads the (filtered) ADP matrix values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_matrix_num_nzes() that receives the
       * filtered matrix values of the local ADP matrix.
       *
       * \param[in] row_ptr
       * A \transient pointer to the row-pointer array of the ADP matrix as returned by the
       * _upload_symbolic() function.
       *
       * \param[in] col_idx
       * A \transient pointer to the column-index array of the ADP matrix as returned by the
       * _upload_symbolic() function.
       */
      template<typename DTV_, typename RPT_, typename CIT_>
      void _upload_numeric(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx)
      {
        const Index num_rows = _system_matrix.rows();
        for(Index i(0); i < num_rows; ++i)
        {
          _system_matrix.get_row_values(i, &val[row_ptr[i]]);
        }

        this->_filter_mat(val, row_ptr, col_idx, this->_system_filter);
      }

      /**
       * \brief Uploads the ADP vector values to the given array
       *
       * \param[out] val
       * A \transient pointer to an array of length _get_adp_vector_size() that receives the
       * local ADP values of the input vector given by \p vector.
       *
       * \param[in] vector
       * A \transient reference to the local vector that is to be uploaded to the ADP value array.
       */
      template<typename DTV_>
      void _upload_vector(DTV_* val, const VectorType& vector)
      {
        vector.set_vec(val);
      }

      /**
       * \brief Downloads the ADP vector values from the given array
       *
       * \param[inout] vector
       * A \transient reference to the local vector that receives the downloaded ADP value array.
       * It is assumed to be allocated to the correct size, but its contents on entry are ignored.
       *
       * \param[in] val
       * A \transient pointer to an array of length _get_adp_vector_size() that contains the
       * ADP values that are to be downloaded to the output vector given by \p vector.
       */
      template<typename DTV_>
      void _download_vector(VectorType& vector, const DTV_* val)
      {
        vector.set_vec_inv(val);
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_>
      void _filter_mat(DTV_*, const RPT_*, const CIT_*,
        const LAFEM::NoneFilter<DataType, IndexType>&)
      {
        // nothing to do here
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_>
      void _filter_mat(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx,
        const LAFEM::UnitFilter<DataType, IndexType>& filter)
      {
        const IndexType n = filter.used_elements();
        const IndexType* fil_idx = filter.get_indices();
        for(IndexType i = 0; i < n; ++i)
        {
          const Index row = fil_idx[i];
          for(RPT_ j = row_ptr[row]; j < row_ptr[row+1]; ++j)
            val[j] = DTV_(Index(col_idx[j]) == row ? 1 : 0);
        }
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_>
      void _filter_mat(DTV_*, const RPT_*, const CIT_*,
        const LAFEM::MeanFilter<DataType, IndexType>& filter)
      {
        XASSERTM(filter.empty(), "LAFEM::MeanFilter is not supported by ADPSolverBase yet!");
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_>
      void _filter_mat(DTV_*, const RPT_*, const CIT_*,
        const Global::MeanFilter<DataType, IndexType>& filter)
      {
        XASSERTM(filter.empty(), "LAFEM::MeanFilter is not supported by ADPSolverBase yet!");
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_, typename SubFilter_>
      void _filter_mat(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx,
        const LAFEM::FilterSequence<SubFilter_>& filter)
      {
        for(const auto& sf : filter)
          this->_filter_mat(val, row_ptr, col_idx, sf.second);
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_, typename First_, typename... Rest_>
      void _filter_mat(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx,
        const LAFEM::FilterChain<First_, Rest_...>& filter)
      {
        this->_filter_mat(val, row_ptr, col_idx, filter.first());
        this->_filter_mat(val, row_ptr, col_idx, filter.rest());
      }

      /// auxiliary function: filters the local ADP matrix
      template<typename DTV_, typename RPT_, typename CIT_, typename First_>
      void _filter_mat(DTV_* val, const RPT_* row_ptr, const CIT_* col_idx,
        const LAFEM::FilterChain<First_>& filter)
      {
        this->_filter_mat(val, row_ptr, col_idx, filter.first());
      }
    }; // class ADPSolverBase<...>
  } // namespace Solver
} // namespace FEAT
