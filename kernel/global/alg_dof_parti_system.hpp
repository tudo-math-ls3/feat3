// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/adjacency/dynamic_graph.hpp>
#include <kernel/global/alg_dof_parti.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/mean_filter_blocked.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/global/mean_filter.hpp>
#include <kernel/lafem/matrix_mirror.hpp>

#include <algorithm>

namespace FEAT
{
  namespace Global
  {
    /// \cond internal
    namespace Intern
    {
      /// auxiliary class: alg-dof-parti data owner matrix mirror
      template<typename Matrix_>
      class ADPDOMM;

      /// auxiliary class for ADP matrix-related functions
      template<typename Matrix_>
      class ADPMatAux;

      /// auxiliary class for ADP filter-related functions
      template<typename Filter_>
      class ADPFilAux;
    }
    /// \endcond

    /**
     * \brief Algebraic DOF partitioning linear system conversion class
     *
     * This class implements the functionality that is required to convert a linear system
     * represented by a global system matrix, a global system filter and a global RHS vector from
     * the domain decomposition approach to a row-wise partitioned global CSR (RPG-CSR) system
     * based on the algebraic dof-partitioning, which can then be used as input for various
     * MPI-parallel third-party library solvers.
     *
     * This class is used by the Solver::ADPSolverBase class to do the dirty work of actually
     * converting the linear system. Unless you intend to work with the ADP system conversion on
     * a lower level yourself, it is recommended to work with the Solver::ADPSolverBase class
     * rather than this class directly.
     *
     * This class is implemented in a way that the user has to allocate and manage the necessary
     * arrays of the RPG-CSR matrix as well as the vectors manually, whereas this class only offers
     * functions that determine the dimensions of the corresponding arrays as well as functions
     * that fill the arrays with the proper entries.
     *
     * In the context of this class, the term \b upload always refers to the conversion of data
     * from the Global::Matrix, Global::Filter and/or Global::Vector objects to their corresponding
     * ADP array counterparts, whereas \b download refers to the conversion of an ADP vector array
     * to its corresponding Global::Vector object counterpart.
     *
     * This class is meant to be used as follows:
     * -# Create an object of this class for a given system matrix and system filter object.
     *    Optionally, you can also specify an already assembled Global::AlgDofParti object to the
     *    constructor, which obviously has to match the system matrix structure.
     * -# Once the matrix structure is assembled, call the init_symbolic() function.
     *    If you did not pass an AlgDofParti object to the constructor of this class, then a new
     *    AlgDofParti object will be created internally based on the row-gate of the system matrix
     *    during this phase automatically.
     * -# Allocate the following arrays for the linear system:
     *    - An integer row-pointer array of length get_adp_matrix_rows()+1 entries.
     *    - An integer column-index array of length get_adp_matrix_nzes() entries.
     *    - A floating point matrix value array of length get_adp_matrix_nzes() entries.
     *    - One or two ADP vector value arrays of length get_num_owned_dofs() entries.
     * -# Fill the row-pointer and column-index arrays by calling the upload_matrix_symbolic()
     *    function.
     * -# Once the system matrix values are assembled, fill the matrix values array by calling
     *    the upload_matrix_numeric() function.
     * -# Once the system filter values are assembled, upload the filter values to the ADP system
     *    by calling the upload_filter() function.
     * -# Once both the matrix and filter values have been uploaded (note that both of these steps
     *    are independent of each other and can be performed in any order), apply the uploaded
     *    filter to the matrix value array by calling the filter_matrix() function.
     * -# Once the system defect vector is assembled, upload its values to the ADP defect vector
     *    value array by calling the upload_vector() function.
     * -# Apply the defect filter by calling the filter_vec_def() function.
     * -# Now both the matrix and the corresponding defect vector have been converted fully, so you
     *    are now free to apply whatever solver algorithm you want onto them.
     * -# Once you have obtained a correction vector, you may want to apply the correction filter
     *    onto it by calling the filter_vec_cor() function if your solver algorithm has a tendency
     *    of not solving the system exactly, and then download the ADP correction vector value
     *    array to a system correction vector by calling the download_vector() function.
     * -# If the numerical values of the system matrix or the system filter change, e.g. after
     *    reassembling them during a non-linear and/or unsteady simulation, you can simply
          re-upload the matrix and/or filter values by using the above upload functions.
     * -# Once you are done, you can simply destroy the AlgDofPartiSystem object without any prior
     *    cleanup.
     *
     * \note
     * If the matrix structure changes for whatever reason, then you have to destroy both the
     * AlgDofPartiSystem object as well as its underlying AlgDofParti object and re-create them
     * for the new matrix structure from scratch.
     *
     * \attention
     * In its current state, only unit-filters (both scalar and blocked) and filter chains and
     * sequences thereof are supported by this class, i.e. neither mean filters nor slip filters
     * are supported yet. However, note that mean filters and slips filters are allowed to exist
     * in a filter chain passed to this class, but they \b must be empty and therefore inactive,
     * otherwise this class will abort with an assertion failure!
     *
     * \todo add support for LAFEM::PowerXXX containers
     * \todo add support for LAFEM/Global::MeanFilter
     * \todo add support for LAFEM::SlipFilter
     *
     * \author Peter Zajac
     */
    template<typename LocalMatrix_, typename LocalFilter_, typename Mirror_>
    class AlgDofPartiSystem
    {
    public:
      /// our data type
      typedef typename LocalMatrix_::DataType DataType;
      /// our index type
      typedef typename LocalMatrix_::IndexType IndexType;

      /// the local matrix type
      typedef LocalMatrix_ LocalMatrixType;
      /// the vector mirror type
      typedef Mirror_ MirrorType;
      /// the global matrix type
      typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;

      /// the local vector type
      typedef typename LocalMatrixType::VectorTypeR LocalVectorType;

      /// the local filter type
      typedef LocalFilter_ LocalFilterType;

      /// the global filter type
      typedef Global::Filter<LocalFilterType, MirrorType> GlobalFilterType;

      /// the buffer vector type used for communication
      typedef LAFEM::DenseVector<DataType, IndexType> BufferVectorType;
      /// the buffer matrix type used for communication
      typedef LAFEM::MatrixMirrorBuffer<DataType, IndexType> BufferMatrixType;

      /// the matrix type used for the internal storage of our owned matrix rows
      typedef LAFEM::SparseMatrixCSR<DataType, IndexType> OwnedMatrixType;

      /// the algebraic DOF partitioning
      typedef Global::AlgDofParti<LocalVectorType, MirrorType> AlgDofPartiType;

      /// auxiliary index vector type
      typedef typename AlgDofPartiType::IndexVectorType IndexVectorType;

      /// donee neighbor mirror type
      typedef LAFEM::VectorMirror<DataType, IndexType> DoneeMirrorType;
      /// owner neighbor mirror type
      typedef Intern::ADPDOMM<LocalMatrixType> OwnerMirrorType;

      typedef Intern::ADPMatAux<LocalMatrixType> ADPMatAuxType;
      typedef Intern::ADPFilAux<LocalFilterType> ADPFilAuxType;

    //protected:
      /// the global matrix that we want to convert
      const GlobalMatrixType& _global_matrix;
      /// the global filter that we want to convert
      const GlobalFilterType& _global_filter;

      /// the algebraic dof-partitioning
      std::shared_ptr<AlgDofPartiType> _alg_dof_parti;

      /// adjacency graph containing the owned matrix structure
      std::shared_ptr<Adjacency::Graph> _owned_graph;

      /// number of non-zero entries owned by this process
      Index _num_owned_non_zeros;

      /// number of non-zero entries in total
      Index _num_global_non_zeros;

      /// the matrix buffers for our donee and owner neighbors; both buffers use global DOF indices
      std::vector<BufferMatrixType> _donee_bufs;
      std::vector<BufferMatrixType> _owner_bufs;

      /// the data array mirrors for our donee and owner neighbors
      std::vector<DoneeMirrorType> _donee_data_mirs;
      std::vector<OwnerMirrorType> _owner_data_mirs;

      /// the data array mirror for this process
      OwnerMirrorType _owned_data_mir;

      /// indices of all unit-filtered matrix rows; we store this as a mirror
      LAFEM::UnitFilter<DataType, IndexType> _unit_filter_rows;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A \resident reference to the global system matrix
       *
       * \param[in] filter
       * A \resident reference to the global system filter
       */
      explicit AlgDofPartiSystem(const GlobalMatrixType& matrix, const GlobalFilterType& filter) :
        _global_matrix(matrix),
        _global_filter(filter),
        _alg_dof_parti(),
        _owned_graph(),
        _num_owned_non_zeros(0u),
        _num_global_non_zeros(0u)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A \resident reference to the global system matrix
       *
       * \param[in] filter
       * A \resident reference to the global system filter
       *
       * \param[in] adp
       * A shared pointer containing the ADP object to be used
       */
      explicit AlgDofPartiSystem(const GlobalMatrixType& matrix, const GlobalFilterType& filter,
        std::shared_ptr<AlgDofPartiType> adp) :
        _global_matrix(matrix),
        _global_filter(filter),
        _alg_dof_parti(adp),
        _owned_graph(),
        _num_owned_non_zeros(0u),
        _num_global_non_zeros(0u)
      {
      }

      /// no copy, no problems
      AlgDofPartiSystem(const AlgDofPartiSystem&) = delete;
      /// no copy, no problems
      AlgDofPartiSystem& operator=(const AlgDofPartiSystem&) = delete;

      /// \returns The size of this object in bytes
      std::size_t bytes() const
      {
        std::size_t r = 0u;
        if(this->_alg_dof_parti)
          r += this->_alg_dof_parti->bytes();
        if(this->_owned_graph)
          r += this->_owned_graph->bytes();
        for(const auto& x : this->_donee_bufs)
          r += x.bytes();
        for(const auto& x : this->_owner_bufs)
          r += x.bytes();
        for(const auto& x : this->_donee_data_mirs)
          r += x.bytes();
        for(const auto& x : this->_owner_data_mirs)
          r += x.bytes();
        r += _owned_data_mir.bytes();
        r += _unit_filter_rows.bytes();
        return r;
      }

      /**
       * \brief Returns the algebraic dof partitioning that is used internally.
       *
       * \note
       * The ADP object is (in general) allocated and initialized by the init_symbolic() function,
       * so this function will (in general) return a null pointer if called before init_symbolic().
       */
      std::shared_ptr<AlgDofPartiType> get_alg_dof_parti()
      {
        return this->_alg_dof_parti;
      }

      /// \returns The number of local DOFs owned by this process.
      Index get_num_owned_dofs() const
      {
        return this->_alg_dof_parti->get_num_owned_dofs();
      }

      /// \returns The total number of global DOFs.
      Index get_num_global_dofs() const
      {
        return this->_alg_dof_parti->get_num_global_dofs();
      }

      /// \returns The index of the first global DOF owned by this process.
      Index get_global_dof_offset() const
      {
        return this->_alg_dof_parti->get_global_dof_offset();
      }

      /// \returns The size of a local ADP vector (corresponds to the number of owned DOFs).
      Index get_adp_vector_size() const
      {
        return this->_alg_dof_parti->get_num_owned_dofs();
      }

      /// \returns The number of rows of the local ADP matrix (corresponds to the number of owned DOFs).
      Index get_adp_matrix_rows() const
      {
        return this->_alg_dof_parti->get_num_owned_dofs();
      }

      /// \returns The number of columns of the local ADP matrix (corresponds to the total number of global DOFs).
      Index get_adp_matrix_cols() const
      {
        return this->_alg_dof_parti->get_num_global_dofs();
      }

      /// \returns The number of nonzero entries of the local ADP matrix
      Index get_adp_matrix_nzes() const
      {
        return this->_num_owned_non_zeros;
      }

      /// \returns The total number of nonzero entries of the global ADP matrix
      Index get_num_global_nonzeros() const
      {
        return this->_num_global_non_zeros;
      }

      /// \returns A pointer to the underlying communicator.
      const Dist::Comm* get_comm() const
      {
        XASSERT(this->_alg_dof_parti != nullptr);
        return this->_alg_dof_parti->get_comm();
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
      void upload_vector(DTV_* val, const LocalVectorType& vector)
      {
        XASSERT(this->_alg_dof_parti != nullptr);
        this->_alg_dof_parti->upload_vector(val, vector);
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
      void download_vector(LocalVectorType& vector, const DTV_* val)
      {
        XASSERT(this->_alg_dof_parti != nullptr);
        this->_alg_dof_parti->download_vector(val, vector);
      }

      /**
       * \brief Performs a global matrix-vector product with this matrix.
       *
       * \attention
       * As convenient and tempting as this function may look, please be aware that
       * using this function will completely destroy the scalability of our code when
       * running on a large number of processes. In addition to the scalability issues,
       * please be aware that this function may cause a fatal out-of-memory event, as
       * the underlying (naive) implementation gathers the full multiplicand vector
       * on all processes using the AlgDofPartiVector::allgather() function.\n
       * Do not use this function except for <b>small-scale debugging</b> purposes. You have been warned.
       *
       * \param[out] vec_r
       * A reference to the algebraic-dof-partitioned vector that receives
       * the result of the matrix-vector product.
       *
       * \param[in] vec_x
       * A const reference to the algebraic-dof-partitioned vector that
       * acts as a multiplicand vector.
       */
      template<typename DTX_, typename RPT_, typename CIT_>
      void apply(DTX_* vec_r, const DTX_* vec_x, const DTX_* val_a, const RPT_* row_ptr, const CIT_* col_idx) const
      {
        XASSERTM(!this->_alg_dof_parti->get_all_global_dof_counts().empty(),
          "You did not ask to assemble the required AlgDofParti allgather data");

        // create and gather full multiplicand vector
        std::vector<DTX_> vec_full(this->_alg_dof_parti->get_num_global_dofs());
        this->get_comm()->allgatherv(vec_x,
          this->_alg_dof_parti->get_num_owned_dofs(),
          vec_full.data(),
          this->_alg_dof_parti->get_all_global_dof_counts().data(),
          this->_alg_dof_parti->get_all_global_dof_offsets().data());

        // apply our owned matrix
        LAFEM::Arch::Apply::csr_generic(vec_r, DTX_(1), vec_full.data(), DTX_(0), vec_r, val_a, col_idx, row_ptr,
          this->get_adp_matrix_rows(), this->get_adp_matrix_cols(), this->get_adp_matrix_nzes(), false);
      }

      /**
       * \brief Performs the symbolic initialization
       *
       * This function assembles the algebraic dof partitioning based on the gate of the system matrix
       * (unless an already assembled ADP object was passed to the constructor) and assembles the
       * required buffers and mirrors for the ADP system conversion.
       */
      void init_symbolic()
      {
        if(this->_alg_dof_parti == nullptr)
        {
          this->_alg_dof_parti = std::make_shared<AlgDofPartiType>();
          this->_alg_dof_parti->assemble_by_gate(*this->_global_matrix.get_row_gate());
        }

        this->_assemble_buffers();
        this->_assemble_structure();
        this->_assemble_data_mirrors();
      }

      /**
       * \brief Uploads the row-pointer and column-index arrays of the ADP CSR system matrix
       *
       * \param[inout] row_ptr
       * A \transient pointer to an array containing at least get_adp_matrix_rows()+1 entries that
       * receives the CSR row-pointer array of the local slice of the row-wise partitioned global
       * CSR matrix. The array must be allocated to the correct size by the caller, but the contents
       * of the array upon entry are ignored.
       *
       * \param[in] col_idx
       * A \transient pointer to an array containing at least get_adp_matrix_nzes() entries that
       * receives the CSR column-index array of the local slice of the row-wise partitioned global
       * CSR matrix. The array must be allocated to the correct size by the caller, but the contents
       * of the array upon entry are ignored.
       *
       * \param[in] keep_graph
       * If set to \c false, the internal adjacency graph of the matrix is discarded at the end of
       * this function call, as it serves no further internal purpose and therefore there is no
       * need to keep the graph in memory unless the user wants to access it after the symbolic
       * upload for some debugging purposes.
       *
       * \note
       * The column indices in each row are always sorted in ascending order.
       */
      template<typename RPT_, typename CIT_>
      void upload_matrix_symbolic(RPT_* row_ptr, CIT_* col_idx, bool keep_graph = false)
      {
        XASSERT(row_ptr != nullptr);
        XASSERT(col_idx != nullptr);
        XASSERT(this->_owned_graph != nullptr);

        // maximum allowed row-pointer/column index values assuming signed int types
        static constexpr std::uint64_t max_rpt = 1ull << (8*sizeof(RPT_) - 1);
        static constexpr std::uint64_t max_cit = 1ull << (8*sizeof(CIT_) - 1);

        // prevent "unused variable" warnings in non-debug builds
        (void)max_rpt;
        (void)max_cit;

        const Index* dom_ptr = this->_owned_graph->get_domain_ptr();
        const Index* img_idx = this->_owned_graph->get_image_idx();

        const Index n = this->_owned_graph->get_num_nodes_domain();
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0u; i <= n; ++i)
        {
          ASSERTM(std::uint64_t(dom_ptr[i]) < max_rpt, "row-pointer exceeds RPT_ type range!");
          row_ptr[i] = RPT_(dom_ptr[i]);
        }

        const Index m = this->_owned_graph->get_num_indices();
        FEAT_PRAGMA_OMP(parallel for)
        for(Index j = 0; j < m; ++j)
        {
          ASSERTM(std::uint64_t(img_idx[j]) < max_cit, "column-index exceeds CIT_ type range!");
          col_idx[j] = CIT_(img_idx[j]);
        }

        // release owned graph
        if(!keep_graph)
          this->_owned_graph.reset();
      }

      /**
       * \brief Uploads the data value array of the ADP CSR system matrix
       *
       * \param[inout] val
       * A \transient pointer to an array containing at least get_adp_matrix_nzes() entries that
       * receives the CSR data value array of the local slice of the row-wise partitioned global
       * CSR matrix. The array must be allocated to the correct size by the caller, but the contents
       * of the array upon entry are ignored.
       *
       * \param[in] row_ptr
       * A \transient pointer to the row-pointer array that has been initialized by the
       * upload_matrix_symbolic() function.
       *
       * \param[in] col_idx
       * A \transient pointer to the column-index array that has been initialized by the
       * upload_matrix_symbolic() function.
       */
      template<typename DTA_, typename RPT_, typename CIT_>
      void upload_matrix_numeric(DTA_* val, const RPT_* row_ptr, const CIT_* col_idx)
      {
        this->upload_matrix_numeric(val, row_ptr, col_idx, this->_global_matrix.local());
      }

      /**
       * \brief Copies the matrix entry values from the input matrix into this algebraic-dof-partitioned matrix.
       *
       * \param[inout] val
       * A \transient pointer to an array containing at least get_adp_matrix_nzes() entries that
       * receives the CSR data value array of the local slice of the row-wise partitioned global
       * CSR matrix. The array must be allocated to the correct size by the caller, but the contents
       * of the array upon entry are ignored.
       *
       * \param[in] row_ptr
       * A \transient pointer to the row-pointer array that has been initialized by the
       * upload_matrix_symbolic() function.
       *
       * \param[in] col_idx
       * A \transient pointer to the column-index array that has been initialized by the
       * upload_matrix_symbolic() function.
       *
       * \param[in] matrix
       * The input matrix whose entry values are to be copied.
       *
       * \attention
       * The input matrix is assumed to be an unfiltered type-0 matrix.
       */
      template<typename DTA_, typename RPT_, typename CIT_>
      void upload_matrix_numeric(DTA_* val, const RPT_* DOXY(row_ptr), const CIT_* DOXY(col_idx), const LocalMatrixType& matrix)
      {
        XASSERT(this->_alg_dof_parti != nullptr);

        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // get our communicator
        const Dist::Comm& comm = *this->get_comm();

        // format our internal matrix
        memset(val, 0, sizeof(DTA_)*this->_num_owned_non_zeros);

        const Index num_neigh_owner = adp.get_num_owner_neighbors();
        const Index num_neigh_donee = adp.get_num_donee_neighbors();

        Dist::RequestVector recv_reqs(num_neigh_donee);
        Dist::RequestVector send_reqs(num_neigh_owner);

        std::vector<BufferVectorType> donee_vbufs(num_neigh_donee);
        std::vector<BufferVectorType> owner_vbufs(num_neigh_owner);

        // create receive buffers and post receives
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          // create buffer
          BufferVectorType& buf = donee_vbufs.at(i);
          buf = BufferVectorType(this->_donee_data_mirs.at(i).num_indices());

          // post receive
          recv_reqs[i] = comm.irecv(buf.elements(), buf.size(), adp.get_donee_rank(i));
        }

        // post sends
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          // create buffer
          BufferVectorType& buf = owner_vbufs.at(i);
          buf = BufferVectorType(this->_owner_data_mirs.at(i).num_indices(), DataType(0));

          // gather from mirror
          this->_owner_data_mirs.at(i).gather_owner_data(buf, matrix);

          // post send
          send_reqs[i] = comm.isend(buf.elements(), buf.size(), adp.get_owner_rank(i));
        }

        // upload our own data
        this->_owned_data_mir.upload_owned_data(val, matrix);

        // process all pending receives
        for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
        {
          // scatter from buffer
          this->_scatter_donee_data(val, donee_vbufs.at(idx), this->_donee_data_mirs.at(idx));
        }

        // wait for all sends to finish
        send_reqs.wait_all();
      }

      /// Uploads the values of the system filter to the ADP system
      void upload_filter()
      {
        upload_filter(this->_global_filter.local());
      }

      /**
       * \brief Uploads the values of the given filter to the ADP system
       *
       * \param[in] filter
       * A \transient reference to the filter whose values are to be uploaded
       */
      void upload_filter(const LocalFilterType& filter)
      {
        _unit_filter_rows = LAFEM::UnitFilter<DataType, IndexType>(this->get_num_owned_dofs());
        ADPFilAuxType::upload_filter(_unit_filter_rows, filter, _alg_dof_parti->get_owned_mirror() ,0u);
      }

      /**
       * \brief Applies the ADP filter onto an ADP defect vector
       *
       * \param[inout] vec_def
       * A \transient pointer to the defect vector data array that is to be filtered.
       */
      template<typename DTV_>
      void filter_vec_def(DTV_* vec_def)
      {
        if(_unit_filter_rows.used_elements() <= Index(0))
          return;

        XASSERT(vec_def != nullptr);

        LAFEM::Arch::UnitFilter::filter_def(vec_def, _unit_filter_rows.get_indices(), _unit_filter_rows.used_elements());
      }

      /**
       * \brief Applies the ADP filter onto an ADP correction vector
       *
       * \param[inout] vec_def
       * A \transient pointer to the defect vector data array that is to be filtered.
       */
      template<typename DTV_>
      void filter_vec_cor(DTV_* vec_cor)
      {
        if(_unit_filter_rows.used_elements() <= Index(0))
          return;

        XASSERT(vec_cor != nullptr);

        // UnitFilter has no filter_cor, because it is identical to filter_def
        LAFEM::Arch::UnitFilter::filter_def(vec_cor, _unit_filter_rows.get_indices(), _unit_filter_rows.used_elements());
      }

      /**
       * \brief Applies the ADP filter onto an ADP matrix
       *
       * \param[inout] val
       * A \transient pointer to the value data array of the matrix that is to be filtered.
       *
       * \param[in] row_ptr
       * A \transient pointer to the row-pointer array that has been initialized by the
       * upload_matrix_symbolic() function.
       *
       * \param[in] col_idx
       * A \transient pointer to the column-index array that has been initialized by the
       * upload_matrix_symbolic() function.
       */
      template<typename DTA_, typename RPT_, typename CIT_>
      void filter_matrix(DTA_* val, const RPT_* row_ptr, const CIT_* col_idx)
      {
        if(_unit_filter_rows.used_elements() <= Index(0))
          return;

        XASSERT(val != nullptr);
        XASSERT(row_ptr != nullptr);
        XASSERT(col_idx != nullptr);

        const Index row_off = this->get_alg_dof_parti()->get_global_dof_offset();
        const IndexType* row_idx = this->_unit_filter_rows.get_indices();
        const Index n = this->_unit_filter_rows.used_elements();

        // process all filtered rows
        for(Index i = 0; i < n; ++i)
        {
          const Index row = row_off + row_idx[i];
          RPT_ j = row_ptr[row_idx[i]];
          const RPT_ j_end = row_ptr[row_idx[i]+1];
          for(; j < j_end; ++j)
          {
            val[j] = DTA_(row == Index(col_idx[j]) ? 1 : 0);
          }
        }
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    protected:
      /**
       * \brief Assembles a matrix buffer for a given row mirror
       */
      static BufferMatrixType _asm_mat_buf(
        const LocalMatrixType& local_matrix,
        const MirrorType& row_mirror,
        const IndexVectorType& global_dof_idx,
        const Index num_global_dofs)
      {
        // compute number of rows of matrix buffers
        Index num_buf_rows = ADPMatAuxType::calc_mat_buf_num_rows(local_matrix, row_mirror);

        // allocate a temporary array
        std::vector<IndexType> aux_row(std::size_t(num_buf_rows+1u), IndexType(0));

        // count number of non-zeros per row
        ADPMatAuxType::calc_mat_buf_row_nze(aux_row.data(), local_matrix, row_mirror);

        // build buffer row pointer by exclusive scan
        feat_omp_ex_scan(std::size_t(num_buf_rows+1u), aux_row.data(), aux_row.data());
        IndexType num_buf_nze = aux_row.back();

        // allocate buffer
        BufferMatrixType buffer(num_buf_rows, num_global_dofs, num_buf_nze, 1);
        IndexType* buf_row_ptr = buffer.row_ptr();
        IndexType* buf_col_idx = buffer.col_ind();

        // copy buffer row pointer because we need a disposable copy for the next operation
        for(Index i = 0u; i <= num_buf_rows; ++i)
          buf_row_ptr[i] = aux_row[i];

        // gather the column indices
        // we have to use the auxiliary row pointer here because it is overwritten by the following call
        ADPMatAuxType::gather_mat_buf_col_idx(aux_row.data(), buf_col_idx, local_matrix, row_mirror, global_dof_idx);

        // sort the buffer column indices; this will make things a lot easier and more efficient later on
        for(Index i = 0u; i < num_buf_rows; ++i)
        {
          std::sort(&buf_col_idx[buf_row_ptr[i]], &buf_col_idx[buf_row_ptr[i+1]]);
        }

        // okay
        return buffer;
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      /**
       * \brief Assembles the owner and donee matrix buffers.
       */
      void _assemble_buffers()
      {
        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // get out communicator
        const Dist::Comm& comm = *this->get_comm();

        // get number of neighbors and allocate buffer vectors
        const Index num_neigh_owner = adp.get_num_owner_neighbors();
        const Index num_neigh_donee = adp.get_num_donee_neighbors();
        _owner_bufs.resize(num_neigh_owner);
        _donee_bufs.resize(num_neigh_donee);

        // create the owner buffers by using the auxiliary helper function.
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          // allocate matrix buffer
          this->_owner_bufs.at(i) = _asm_mat_buf(this->_global_matrix.local(),
            adp.get_owner_mirror(i), adp.get_global_dof_indices(), adp.get_num_global_dofs());
        }

        // We have assembled our owner-neighbor matrix buffers, however, we cannot
        // assemble the donee-neighbor matrix buffers on our own. Instead, each owner
        // neighbor has to receive its donee-buffers from its donee-neighbors, because
        // only the donee knows the matrix layout of those buffers.
        // Yes, this is brain-twisting, but believe me, it cannot be realized otherwise.

        // Note:
        // The following code is effectively just copy-&-paste from the SynchMatrix::init()
        // function, as the same buffer problem also arises when converting global type-0
        // matrices to local type-1 matrices.

        Dist::RequestVector send_reqs(num_neigh_owner);
        Dist::RequestVector recv_reqs(num_neigh_donee);

        // receive buffer dimensions vector
        std::vector<std::array<Index,4>> send_dims(num_neigh_owner);
        std::vector<std::array<Index,4>> recv_dims(num_neigh_donee);

        // post send-buffer dimension receives
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          recv_reqs[i] = comm.irecv(recv_dims.at(i).data(), std::size_t(4), adp.get_donee_rank(i));
        }

        // send owner-buffer dimensions
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          const BufferMatrixType& sbuf = this->_owner_bufs.at(i);
          send_dims.at(i)[0] = Index(sbuf.rows());
          send_dims.at(i)[1] = Index(sbuf.columns());
          send_dims.at(i)[2] = Index(sbuf.entries_per_nonzero());
          send_dims.at(i)[3] = Index(sbuf.used_elements());
          send_reqs[i] = comm.isend(send_dims.at(i).data(), std::size_t(4), adp.get_owner_rank(i));
        }

        // wait for all receives to finish
        recv_reqs.wait_all();

        // allocate donee-buffers
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          // get the receive buffer dimensions
          Index nrows = recv_dims.at(i)[0];
          Index ncols = recv_dims.at(i)[1];
          Index nepnz = recv_dims.at(i)[2];
          Index nnze  = recv_dims.at(i)[3];

          // allocate receive buffer
          this->_donee_bufs.at(i) = BufferMatrixType(nrows, ncols, nnze, nepnz);
        }

        // post donee-buffer row-pointer array receives
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          recv_reqs[i] = comm.irecv(this->_donee_bufs.at(i).row_ptr(),
            this->_donee_bufs.at(i).rows()+std::size_t(1), adp.get_donee_rank(i));
        }

        // wait for all previous sends to finish
        send_reqs.wait_all();

        // post owner-buffer row-pointer array sends
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          send_reqs[i] = comm.isend(this->_owner_bufs.at(i).row_ptr(),
            this->_owner_bufs.at(i).rows()+std::size_t(1), adp.get_owner_rank(i));
        }

        // wait for all previous receives to finish
        recv_reqs.wait_all();

        // post donee-buffer column-index array receives
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          recv_reqs[i] = comm.irecv(this->_donee_bufs.at(i).col_ind(),
            this->_donee_bufs.at(i).used_elements(), adp.get_donee_rank(i));
        }

        // wait for all previous sends to finish
        send_reqs.wait_all();

        // post owner-buffer column-index array sends
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          send_reqs[i] = comm.isend(this->_owner_bufs.at(i).col_ind(),
            this->_owner_bufs.at(i).used_elements(), adp.get_owner_rank(i));
        }

        // wait for all receives and sends to finish
        recv_reqs.wait_all();
        send_reqs.wait_all();
      }

      /**
       * \brief Assembles the matrix structure of the algebraic-DOF-partitioned matrix.
       */
      void _assemble_structure()
      {
        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // get our matrix dimensions:
        // * each row of our matrix corresponds to one owned DOF
        // * each column corresponds to one global DOF
        const Index num_rows = adp.get_num_owned_dofs();
        const Index num_cols = adp.get_num_global_dofs();

        // Unfortunately, assembling the matrix structure is not that easy,
        // because it is a union of our owned matrix structure and all of
        // our donee matrix buffers. We don't have a chance to determine how
        // many duplicate (i,j) pairs we will encounter here, so we use the
        // DynamicGraph class here to keep things simple.

        // create a dynamic graph for our matrix
        Adjacency::DynamicGraph dynamic_graph(num_rows, num_cols);

        // gather the owned rows for to obtain the initial structure
        ADPMatAuxType::gather_owned_struct(dynamic_graph, this->_global_matrix.local(),
          adp.get_owned_mirror(), adp.get_global_dof_indices(), Index(0));

        // process all donee-neighbor buffers
        const Index num_neigh_donee = adp.get_num_donee_neighbors();
        for(Index ineigh(0); ineigh < num_neigh_donee; ++ineigh)
        {
          // get the mirror and the matrix buffer of that donee neighbor
          // note: the donee mirror indices are owned DOF indices
          const auto& mir = adp.get_donee_mirror(ineigh);
          const BufferMatrixType& buf = this->_donee_bufs.at(ineigh);
          XASSERT(mir.num_indices() == buf.rows());
          const Index num_idx = mir.num_indices();
          const IndexType* mir_idx = mir.indices();
          const IndexType* buf_ptr = buf.row_ptr();
          const IndexType* buf_idx = buf.col_ind();

          // loop over all matrix buffer rows
          for(Index i(0); i < num_idx; ++i)
          {
            // the owned DOF index of our the buffer matrix row
            const IndexType row = mir_idx[i];

            // loop over all buffer indices
            for(IndexType j(buf_ptr[i]); j < buf_ptr[i + 1]; ++j)
            {
              // insert entry into our local matrix
              dynamic_graph.insert(row, buf_idx[j]);
            }
          }
        }

        // store number of non-zero entries
        this->_num_owned_non_zeros = dynamic_graph.get_num_indices();
        this->get_comm()->allreduce(&this->_num_owned_non_zeros, &this->_num_global_non_zeros, std::size_t(1), Dist::op_sum);

        // render into a standard graph
        this->_owned_graph = std::make_shared<Adjacency::Graph>(Adjacency::RenderType::as_is, dynamic_graph);
        dynamic_graph.clear();
      }

      /**
       * \brief Auxiliary function: assembles a data-mirror for a donee-neighbor
       */
      static void _asm_neighbor_donee_data_mir(
        DoneeMirrorType& data_mirror,
        const Adjacency::Graph& owned_graph,
        const DoneeMirrorType& mirror,
        const BufferMatrixType& buffer)
      {
        const Index buf_rows = buffer.rows();
        const Index num_idx = buffer.used_elements();
        const IndexType* row_idx = mirror.indices();
        const Index* row_ptr = owned_graph.get_domain_ptr();
        const Index* col_idx = owned_graph.get_image_idx();
        const IndexType* buf_ptr = buffer.row_ptr();
        const IndexType* buf_idx = buffer.col_ind();

        // allocate data mirror
        data_mirror = DoneeMirrorType(owned_graph.get_num_indices(), num_idx);
        IndexType* dat_idx = data_mirror.indices();

        // loop over all buffer matrix rows
        for(Index i(0); i < buf_rows; ++i)
        {
          IndexType       j     = buf_ptr[i];
          const IndexType j_end = buf_ptr[i + 1];
          Index       k     = row_ptr[row_idx[i]];
          const Index k_end = row_ptr[row_idx[i] + 1];

          // loop over local matrix and buffer matrix row in a merge-style fashion
          while((j < j_end) && (k < k_end))
          {
            if(buf_idx[j] < col_idx[k])
              ++j;
            else if(buf_idx[j] > col_idx[k])
              ++k;
            else
            {
              dat_idx[j] = IndexType(k);
              ++j;
              ++k;
            }
          }
        }
      }

      /**
       * \brief Assembles the data-mirrors for all owner- and donee-neighbors.
       *
       * This functions assembles the data-mirrors, which are used for the communication
       * of the numerical values of a matrix in the #upload_numeric() function.
       * This allows that the matrix values can be uploaded and synchronized without the
       * need to re-communicate the unchanged matrix buffer layouts every time.
       */
      void _assemble_data_mirrors()
      {
        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // assemble donee data mirrors
        this->_donee_data_mirs.resize(this->_donee_bufs.size());
        for(Index i(0); i < this->_donee_bufs.size(); ++i)
          _asm_neighbor_donee_data_mir(this->_donee_data_mirs.at(i), *this->_owned_graph,
            adp.get_donee_mirror(i), this->_donee_bufs.at(i));

        // assemble owner data mirrors
        this->_owner_data_mirs.resize(this->_owner_bufs.size());
        for(Index i(0); i < this->_owner_bufs.size(); ++i)
          this->_owner_data_mirs.at(i).asm_neighbor_owner_data_mir(this->_global_matrix.local(),
            adp.get_owner_mirror(i), this->_owner_bufs.at(i),
            adp.get_global_dof_indices(), Index(0));

        // assemble our owned data mirror
        this->_owned_data_mir.asm_owned_data_mir(this->_global_matrix.local(), *this->_owned_graph,
          adp.get_owned_mirror(), adp.get_global_dof_indices(), Index(0));
      }

      /**
       * \brief Scatters the values of a donee-neighbor buffer-vector into a local matrix
       */
      template<typename DTA_>
      static void _scatter_donee_data(DTA_* adp_val, const BufferVectorType& buffer,
        const DoneeMirrorType& data_mirror)
      {
        XASSERT(buffer.size() == data_mirror.num_indices());

        const DataType* buf_val = buffer.elements();
        const IndexType* mir_idx = data_mirror.indices();
        const Index n = buffer.size();
        for(Index i(0); i < n; ++i)
        {
          adp_val[mir_idx[i]] += DTA_(buf_val[i]);
        }
      }
    }; // class AlgDofPartiSystem<...>

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// \cond internal
    namespace Intern
    {
      template<typename DT_, typename IT_>
      class ADPDOMM<LAFEM::SparseMatrixCSR<DT_, IT_>>
      {
      public:
        LAFEM::VectorMirror<DT_, IT_> _buf_mir, _loc_mir;

        Index num_indices() const
        {
          return _buf_mir.num_indices();
        }

        Index buf_size() const
        {
          return _buf_mir.size();
        }

        Index loc_size() const
        {
          return _loc_mir.size();
        }

        String dump() const
        {
          String s;
          XASSERT(_buf_mir.num_indices() == _loc_mir.num_indices());
          const IT_* buf_idx = _buf_mir.indices();
          const IT_* loc_idx = _loc_mir.indices();
          for(Index i = 0; i < _buf_mir.num_indices(); ++i)
            s += "  " + stringify(buf_idx[i]) + ":" + stringify(loc_idx[i]);
          return "[" + s + "]*" + stringify(_buf_mir.num_indices());
        }

        /**
         * \brief Auxiliary function: assembles a data-mirror for an owner-neighbor
         *
         * \param[in] local_matrix
         * The local matrix that acts as a layout template.\n
         * Its layout must be initialized, but its numerical values are ignored.
         *
         * \param[in] mirror
         * The mirror of the owner-neighbor.
         *
         * \param[in] buffer
         * The matrix buffer for the owner-mirror.
         *
         * \param[in] glob_dof_idx
         * A vector containing the global DOF indices for all local DOFs.
         */
        Index asm_neighbor_owner_data_mir(
          const LAFEM::SparseMatrixCSR<DT_, IT_>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buffer,
          const LAFEM::DenseVector<IT_, IT_>& global_dof_idx,
          Index row_offset)
        {
          // get our mirror and matrix arrays
          const IT_* row_idx = row_mirror.indices();
          const IT_* row_ptr = local_matrix.row_ptr();
          const IT_* col_idx = local_matrix.col_ind();
          const IT_* buf_ptr = buffer.row_ptr();
          const IT_* buf_idx = buffer.col_ind();
          const IT_* glob_dof_idx = global_dof_idx.elements();
          const Index row_mir_num_idx = row_mirror.num_indices();

          // count the number of non-zeros that this matrix contributes to the buffer
          Index data_mir_num_idx = 0u;
          for(Index i = 0; i < row_mir_num_idx; ++i)
            data_mir_num_idx += row_ptr[row_idx[i]+1u] - row_ptr[row_idx[i]];

          XASSERT(data_mir_num_idx <= buffer.used_elements());

          // allocate data mirror
          std::vector<std::pair<IT_, IT_>> aux_data_mir;
          aux_data_mir.reserve(data_mir_num_idx);

          // vector for row pointers sorted by global indices
          std::vector<IT_> loc_it;
          loc_it.reserve(local_matrix.columns());

          // declare a lambda that sorts row pointers by global column indices
          auto sort_rel_loc = [&glob_dof_idx, &col_idx] (IT_ a, IT_ b)
          {
            return glob_dof_idx[col_idx[a]] < glob_dof_idx[col_idx[b]];
          };

          // loop over all buffer matrix rows
          for(Index i(0); i < row_mir_num_idx; ++i)
          {
            // get local matrix row index
            IT_       j     = buf_ptr[row_offset + i];
            const IT_ j_end = buf_ptr[row_offset + i + 1];
            const IT_ k_beg = row_ptr[row_idx[i]];
            const IT_ k_end = row_ptr[row_idx[i] + 1];
            loc_it.clear();
            for(IT_ k = k_beg; k < k_end; ++k)
              loc_it.push_back(k);
            std::sort(loc_it.begin(), loc_it.end(), sort_rel_loc);

            auto kit = loc_it.begin();

            // loop over local matrix and buffer matrix row in a merge-style fashion
            while((j < j_end) && (kit != loc_it.end()))
            {
              if(buf_idx[j] < glob_dof_idx[col_idx[*kit]])
                ++j;
              else if(buf_idx[j] > glob_dof_idx[col_idx[*kit]])
                ++kit;
              else
              {
                aux_data_mir.push_back(std::make_pair(j, *kit));
                ++j;
                ++kit;
              }
            }
          }

          // create data mirrors
          Index num_idx = Index(aux_data_mir.size());
          this->_buf_mir = LAFEM::VectorMirror<DT_, IT_>(buffer.used_elements(), num_idx);
          this->_loc_mir = LAFEM::VectorMirror<DT_, IT_>(
            local_matrix.template used_elements<LAFEM::Perspective::pod>(), num_idx);

          IT_* bm_idx = this->_buf_mir.indices();
          IT_* lm_idx = this->_loc_mir.indices();
          for(Index i = 0; i < num_idx; ++i)
          {
            bm_idx[i] = aux_data_mir[i].first;
            lm_idx[i] = aux_data_mir[i].second;
          }

          return row_mir_num_idx;
        }

        /**
         * \brief Auxiliary function: assembles a data-mirror for this process owned DOF's
         *
         * \param[in] local_matrix
         * The local matrix that acts as a layout template.\n
         * Its layout must be initialized, but its numerical values are ignored.
         *
         * \param[in] mirror
         * The mirror of the owner-neighbor.
         *
         * \param[in] glob_dof_idx
         * A vector containing the global DOF indices for all local DOFs.
         */
        Index asm_owned_data_mir(
          const LAFEM::SparseMatrixCSR<DT_, IT_>& local_matrix,
          const Adjacency::Graph& owned_graph,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVector<IT_, IT_>& global_dof_idx,
          Index row_offset)
        {
          // get matrix arrays
          const Index* own_row_ptr = owned_graph.get_domain_ptr();
          const Index* own_col_idx = owned_graph.get_image_idx();
          const IT_* loc_row_ptr = local_matrix.row_ptr();
          const IT_* loc_col_idx = local_matrix.col_ind();

          // allocate data mirror
          std::vector<std::pair<IT_, IT_>> aux_data_mir;
          aux_data_mir.reserve(owned_graph.get_num_indices());

          const Index num_owned_dofs = row_mirror.num_indices();
          const IT_* loc_dof_idx = row_mirror.indices();
          const IT_* glob_dof_idx = global_dof_idx.elements();

          //XASSERT(num_owned_dofs == adp_matrix.rows());

          // row iterator deque
          std::deque<IT_> row_it;

          // declare a lambda that sorts row pointers by global column indices
          auto sort_rel = [&glob_dof_idx, &loc_col_idx] (IT_ a, IT_ b)
          {
            return glob_dof_idx[loc_col_idx[a]] < glob_dof_idx[loc_col_idx[b]];
          };

          // loop over all our owned DOFS
          for(Index own_dof(0); own_dof < num_owned_dofs; ++own_dof)
          {
            const IT_ k_beg = loc_row_ptr[loc_dof_idx[own_dof]];
            const IT_ k_end = loc_row_ptr[loc_dof_idx[own_dof] + 1];
            row_it.clear();
            for(IT_ k = k_beg; k < k_end; ++k)
              row_it.push_back(k);
            std::sort(row_it.begin(), row_it.end(), sort_rel);

            // get the local DOF index for this owned DOF
            IT_       j     = IT_(own_row_ptr[row_offset + own_dof]);
            const IT_ j_end = IT_(own_row_ptr[row_offset + own_dof + 1]);
            auto kit = row_it.begin();

            // loop over local matrix and owned matrix row in a merge-style fashion
            while((j < j_end) && (kit != row_it.end()))
            {
              if(own_col_idx[j] < glob_dof_idx[loc_col_idx[*kit]])
                ++j;
              else if(own_col_idx[j] > glob_dof_idx[loc_col_idx[*kit]])
                ++kit;
              else
              {
                aux_data_mir.push_back(std::make_pair(j, *kit));
                ++j;
                ++kit;
              }
            }
          }

          // create data mirrors
          Index num_idx = Index(aux_data_mir.size());
          this->_buf_mir = LAFEM::VectorMirror<DT_, IT_>(local_matrix.used_elements(), num_idx);
          this->_loc_mir = LAFEM::VectorMirror<DT_, IT_>(local_matrix.used_elements(), num_idx);

          IT_* bm_idx = this->_buf_mir.indices();
          IT_* lm_idx = this->_loc_mir.indices();
          for(Index i = 0; i < num_idx; ++i)
          {
            bm_idx[i] = aux_data_mir[i].first;
            lm_idx[i] = aux_data_mir[i].second;
          }

          return num_owned_dofs;
        }

        template<typename DTA_>
        void upload_owned_data(DTA_* adp_val, const LAFEM::SparseMatrixCSR<DT_, IT_>& local_matrix) const
        {
          const DT_* loc_val = local_matrix.val();
          const IT_* adp_idx = _buf_mir.indices();
          const IT_* loc_idx = _loc_mir.indices();
          XASSERT(_buf_mir.num_indices() == _loc_mir.num_indices());
          Index n = _buf_mir.num_indices();
          for(Index i = 0; i < n; ++i)
          {
            ASSERT(loc_idx[i] < local_matrix.used_elements());
            adp_val[adp_idx[i]] = DTA_(loc_val[loc_idx[i]]);
          }
        }

        /**
         * \brief Gathers the values of a local matrix into an owner-neighbor buffer-vector
         *
         * \param[out] buffer
         * A buffer vector that receives the gathered values.
         *
         * \param[in] local_matrix
         * The local matrix whose values are to be gathered.
         *
         * \param[in] data_mirror
         * The owner-data mirror to be used.
         */
        void gather_owner_data(LAFEM::DenseVector<DT_, IT_>& buffer, const LAFEM::SparseMatrixCSR<DT_, IT_>& local_matrix) const
        {
          DT_* buf_val = buffer.elements();
          const DT_* loc_val = local_matrix.val();
          const IT_* buf_idx = this->_buf_mir.indices();
          const IT_* loc_idx = this->_loc_mir.indices();
          XASSERT(this->_buf_mir.num_indices() == this->_loc_mir.num_indices());
          Index n = this->_buf_mir.num_indices();
          for(Index i = 0; i < n; ++i)
          {
            ASSERT(buf_idx[i] < buffer.size());
            ASSERT(loc_idx[i] < local_matrix.used_elements());
            buf_val[buf_idx[i]] = loc_val[loc_idx[i]];
          }
        }

      }; // class ADPDOMM<LAFEM::SparseMatrixCSR<DT_, IT_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_, int bh_, int bw_>
      class ADPDOMM<LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>>
      {
      public:
        LAFEM::VectorMirror<DT_, IT_> _buf_mir, _loc_mir;

        Index num_indices() const
        {
          return _buf_mir.num_indices();
        }

        Index buf_size() const
        {
          return _buf_mir.size();
        }

        Index loc_size() const
        {
          return _loc_mir.size();
        }

        String dump() const
        {
          String s;
          XASSERT(_buf_mir.num_indices() == _loc_mir.num_indices());
          const IT_* buf_idx = _buf_mir.indices();
          const IT_* loc_idx = _loc_mir.indices();
          for(Index i = 0; i < _buf_mir.num_indices(); ++i)
            s += "  " + stringify(buf_idx[i]) + ":" + stringify(loc_idx[i]);
          return "[" + s + "]*" + stringify(_buf_mir.num_indices());
        }

        Index asm_neighbor_owner_data_mir(
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buffer,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bw_>& global_dof_idx,
          Index row_offset)
        {
          // get our mirror and matrix arrays
          const IT_* row_idx = row_mirror.indices();
          const IT_* row_ptr = local_matrix.row_ptr();
          const IT_* col_idx = local_matrix.col_ind();
          const IT_* buf_ptr = buffer.row_ptr();
          const IT_* buf_idx = buffer.col_ind();
          const auto* glob_dof_idx = global_dof_idx.elements();
          const Index row_mir_num_idx = row_mirror.num_indices();

          // count the number of non-zeros that this matrix contributes to the buffer
          Index data_mir_num_idx = 0u;
          for(Index i = 0; i < row_mir_num_idx; ++i)
            data_mir_num_idx += row_ptr[row_idx[i]+1u] - row_ptr[row_idx[i]];
          data_mir_num_idx *= Index(bh_*bw_);

          XASSERT(data_mir_num_idx <= buffer.used_elements());

          // allocate data mirror
          std::vector<std::pair<IT_, IT_>> aux_data_mir;
          aux_data_mir.reserve(data_mir_num_idx);

          // vector for row pointers sorted by global indices
          std::vector<IT_> loc_it;
          loc_it.reserve(local_matrix.columns() * std::size_t(bh_*bw_));

          // declare a lambda that sorts row pointers by global column indices
          auto sort_rel_loc = [&glob_dof_idx, &col_idx] (IT_ a, IT_ b)
          {
            return glob_dof_idx[col_idx[a]][0] < glob_dof_idx[col_idx[b]][0];
          };

          // loop over all buffer matrix rows
          for(Index i(0); i < row_mir_num_idx; ++i)
          {
            const IT_ k_beg = row_ptr[row_idx[i]];
            const IT_ k_end = row_ptr[row_idx[i] + 1];
            loc_it.clear();
            for(IT_ k = k_beg; k < k_end; ++k)
              loc_it.push_back(k);
            std::sort(loc_it.begin(), loc_it.end(), sort_rel_loc);

            // loop over the block height
            for(int bi = 0; bi < bh_; ++bi)
            {
              // get local matrix row index
              IT_ j = buf_ptr[row_offset + i*IT_(bh_) + IT_(bi)];
              const IT_ j_end = buf_ptr[row_offset + i*IT_(bh_) + IT_(bi) + 1];

              auto kit = loc_it.begin();

              // loop over local matrix and buffer matrix row in a merge-style fashion
              while((j < j_end) && (kit != loc_it.end()))
              {
                if(buf_idx[j] < glob_dof_idx[col_idx[*kit]][0])
                  ++j;
                else if(buf_idx[j] > glob_dof_idx[col_idx[*kit]][0])
                  ++kit;
                else
                {
                  // global column indices match
                  for(int bj = 0; bj < bw_; ++bj)
                  {
                    aux_data_mir.push_back(std::make_pair(j, (*kit)*IT_(bh_*bw_) + IT_(bi*bw_) + IT_(bj)));
                    ++j;
                  }
                  ++kit;
                }
              }
            }
          }

          // create data mirrors
          Index num_idx = Index(aux_data_mir.size());
          this->_buf_mir = LAFEM::VectorMirror<DT_, IT_>(buffer.used_elements(), num_idx);
          this->_loc_mir = LAFEM::VectorMirror<DT_, IT_>(
            local_matrix.template used_elements<LAFEM::Perspective::pod>(), num_idx);

          IT_* bm_idx = this->_buf_mir.indices();
          IT_* lm_idx = this->_loc_mir.indices();
          for(Index i = 0; i < num_idx; ++i)
          {
            bm_idx[i] = aux_data_mir[i].first;
            lm_idx[i] = aux_data_mir[i].second;
          }

          return row_mir_num_idx * Index(bh_);
        }

        Index asm_neighbor_owner_data_mir(
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, 1>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buffer,
          const LAFEM::DenseVector<IT_, IT_>& global_dof_idx,
          Index row_offset)
        {
          // get our mirror and matrix arrays
          const IT_* row_idx = row_mirror.indices();
          const IT_* row_ptr = local_matrix.row_ptr();
          const IT_* col_idx = local_matrix.col_ind();
          const IT_* buf_ptr = buffer.row_ptr();
          const IT_* buf_idx = buffer.col_ind();
          const IT_* glob_dof_idx = global_dof_idx.elements();
          const Index row_mir_num_idx = row_mirror.num_indices();

          // count the number of non-zeros that this matrix contributes to the buffer
          Index data_mir_num_idx = 0u;
          for(Index i = 0; i < row_mir_num_idx; ++i)
            data_mir_num_idx += row_ptr[row_idx[i]+1u] - row_ptr[row_idx[i]];
          data_mir_num_idx *= Index(bh_);

          XASSERT(data_mir_num_idx <= buffer.used_elements());

          // allocate data mirror
          std::vector<std::pair<IT_, IT_>> aux_data_mir;
          aux_data_mir.reserve(data_mir_num_idx);

          // vector for row pointers sorted by global indices
          std::vector<IT_> loc_it;
          loc_it.reserve(local_matrix.columns() * std::size_t(bh_));

          // declare a lambda that sorts row pointers by global column indices
          auto sort_rel_loc = [&glob_dof_idx, &col_idx] (IT_ a, IT_ b)
          {
            return glob_dof_idx[col_idx[a]] < glob_dof_idx[col_idx[b]];
          };

          // loop over all buffer matrix rows
          for(Index i(0); i < row_mir_num_idx; ++i)
          {
            const IT_ k_beg = row_ptr[row_idx[i]];
            const IT_ k_end = row_ptr[row_idx[i] + 1];
            loc_it.clear();
            for(IT_ k = k_beg; k < k_end; ++k)
              loc_it.push_back(k);
            std::sort(loc_it.begin(), loc_it.end(), sort_rel_loc);

            // loop over the block height
            for(int bi = 0; bi < bh_; ++bi)
            {
              // get local matrix row index
              IT_ j = buf_ptr[row_offset + i*IT_(bh_) + IT_(bi)];
              const IT_ j_end = buf_ptr[row_offset + i*IT_(bh_) + IT_(bi) + 1];
              auto kit = loc_it.begin();

              // loop over local matrix and buffer matrix row in a merge-style fashion
              while((j < j_end) && (kit != loc_it.end()))
              {
                if(buf_idx[j] < glob_dof_idx[col_idx[*kit]])
                  ++j;
                else if(buf_idx[j] > glob_dof_idx[col_idx[*kit]])
                  ++kit;
                else
                {
                  // global column indices match
                  aux_data_mir.push_back(std::make_pair(j, (*kit)*IT_(bh_) + IT_(bi)));
                  ++j;
                  ++kit;
                }
              }
            }
          }

          // create data mirrors
          Index num_idx = Index(aux_data_mir.size());
          this->_buf_mir = LAFEM::VectorMirror<DT_, IT_>(buffer.used_elements(), num_idx);
          this->_loc_mir = LAFEM::VectorMirror<DT_, IT_>(
            local_matrix.template used_elements<LAFEM::Perspective::pod>(), num_idx);

          IT_* bm_idx = this->_buf_mir.indices();
          IT_* lm_idx = this->_loc_mir.indices();
          for(Index i = 0; i < num_idx; ++i)
          {
            bm_idx[i] = aux_data_mir[i].first;
            lm_idx[i] = aux_data_mir[i].second;
          }

          return row_mir_num_idx * Index(bh_);
        }

        Index asm_owned_data_mir(
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& local_matrix,
          const Adjacency::Graph& owned_graph,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bw_>& global_dof_idx,
          Index row_offset)
        {
          // get matrix arrays
          const Index* own_row_ptr = owned_graph.get_domain_ptr();
          const Index* own_col_idx = owned_graph.get_image_idx();
          const IT_* loc_row_ptr = local_matrix.row_ptr();
          const IT_* loc_col_idx = local_matrix.col_ind();

          // allocate data mirror
          std::vector<std::pair<IT_, IT_>> aux_data_mir;
          aux_data_mir.reserve(owned_graph.get_num_indices() * std::size_t(bh_*bw_));

          const Index num_owned_dofs = row_mirror.num_indices();
          const IT_* loc_dof_idx = row_mirror.indices();
          const auto* glob_dof_idx = global_dof_idx.elements();

          //XASSERT(num_owned_dofs == adp_matrix.rows());

          // row iterator deque
          std::deque<IT_> row_it;

          // declare a lambda that sorts row pointers by global column indices
          auto sort_rel = [&glob_dof_idx, &loc_col_idx] (IT_ a, IT_ b)
          {
            return glob_dof_idx[loc_col_idx[a]][0] < glob_dof_idx[loc_col_idx[b]][0];
          };

          // loop over all our owned DOFS
          for(Index own_dof(0); own_dof < num_owned_dofs; ++own_dof)
          {
            const IT_ k_beg = loc_row_ptr[loc_dof_idx[own_dof]];
            const IT_ k_end = loc_row_ptr[loc_dof_idx[own_dof] + 1];
            row_it.clear();
            for(IT_ k = k_beg; k < k_end; ++k)
              row_it.push_back(k);
            std::sort(row_it.begin(), row_it.end(), sort_rel);

            // loop over the block height
            for(int bi = 0; bi < bh_; ++bi)
            {
              IT_       j     = IT_(own_row_ptr[row_offset + own_dof*IT_(bh_) + IT_(bi)]);
              const IT_ j_end = IT_(own_row_ptr[row_offset + own_dof*IT_(bh_) + IT_(bi) + 1]);
              auto kit = row_it.begin();

              // loop over local matrix and owned matrix row in a merge-style fashion
              while((j < j_end) && (kit != row_it.end()))
              {
                if(own_col_idx[j] < glob_dof_idx[loc_col_idx[*kit]][0])
                  ++j;
                else if(own_col_idx[j] > glob_dof_idx[loc_col_idx[*kit]][0])
                  ++kit;
                else
                {
                  // global column indices match
                  for(int bj = 0; bj < bw_; ++bj)
                  {
                    aux_data_mir.push_back(std::make_pair(j, (*kit)*IT_(bh_*bw_) + IT_(bi*bw_) + IT_(bj)));
                    ++j;
                  }
                  ++kit;
                }
              }
            }
          }

          // create data mirrors
          Index num_idx = Index(aux_data_mir.size());
          this->_buf_mir = LAFEM::VectorMirror<DT_, IT_>(owned_graph.get_num_indices(), num_idx);
          this->_loc_mir = LAFEM::VectorMirror<DT_, IT_>(
            local_matrix.template used_elements<LAFEM::Perspective::pod>(), num_idx);

          IT_* bm_idx = this->_buf_mir.indices();
          IT_* lm_idx = this->_loc_mir.indices();
          for(Index i = 0; i < num_idx; ++i)
          {
            bm_idx[i] = aux_data_mir[i].first;
            lm_idx[i] = aux_data_mir[i].second;
          }

          return num_owned_dofs * Index(bh_);
        }

        Index asm_owned_data_mir(
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, 1>& local_matrix,
          const Adjacency::Graph& owned_graph,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVector<IT_, IT_>& global_dof_idx,
          Index row_offset)
        {
          // get matrix arrays
          const Index* own_row_ptr = owned_graph.get_domain_ptr();
          const Index* own_col_idx = owned_graph.get_image_idx();
          const IT_* loc_row_ptr = local_matrix.row_ptr();
          const IT_* loc_col_idx = local_matrix.col_ind();

          // allocate data mirror
          std::vector<std::pair<IT_, IT_>> aux_data_mir;
          aux_data_mir.reserve(owned_graph.get_num_indices() * std::size_t(bh_));

          const Index num_owned_dofs = row_mirror.num_indices();
          const IT_* loc_dof_idx = row_mirror.indices();
          const IT_* glob_dof_idx = global_dof_idx.elements();

          //XASSERT(num_owned_dofs == adp_matrix.rows());

          // row iterator deque
          std::deque<IT_> row_it;

          // declare a lambda that sorts row pointers by global column indices
          auto sort_rel = [&glob_dof_idx, &loc_col_idx] (IT_ a, IT_ b)
          {
            return glob_dof_idx[loc_col_idx[a]] < glob_dof_idx[loc_col_idx[b]];
          };

          // loop over all our owned DOFS
          for(Index own_dof(0); own_dof < num_owned_dofs; ++own_dof)
          {
            const IT_ k_beg = loc_row_ptr[loc_dof_idx[own_dof]];
            const IT_ k_end = loc_row_ptr[loc_dof_idx[own_dof] + 1];
            row_it.clear();
            for(IT_ k = k_beg; k < k_end; ++k)
              row_it.push_back(k);
            std::sort(row_it.begin(), row_it.end(), sort_rel);

            // loop over the block height
            for(int bi = 0; bi < bh_; ++bi)
            {
              IT_       j     = IT_(own_row_ptr[row_offset + own_dof*IT_(bh_) + IT_(bi)]);
              const IT_ j_end = IT_(own_row_ptr[row_offset + own_dof*IT_(bh_) + IT_(bi) + 1]);
              auto kit = row_it.begin();

              // loop over local matrix and owned matrix row in a merge-style fashion
              while((j < j_end) && (kit != row_it.end()))
              {
                if(own_col_idx[j] < glob_dof_idx[loc_col_idx[*kit]])
                  ++j;
                else if(own_col_idx[j] > glob_dof_idx[loc_col_idx[*kit]])
                  ++kit;
                else
                {
                  // global column indices match
                  aux_data_mir.push_back(std::make_pair(j, (*kit)*IT_(bh_) + IT_(bi)));
                  ++j;
                  ++kit;
                }
              }
            }
          }

          // create data mirrors
          Index num_idx = Index(aux_data_mir.size());
          this->_buf_mir = LAFEM::VectorMirror<DT_, IT_>(owned_graph.get_num_indices(), num_idx);
          this->_loc_mir = LAFEM::VectorMirror<DT_, IT_>(
            local_matrix.template used_elements<LAFEM::Perspective::pod>(), num_idx);

          IT_* bm_idx = this->_buf_mir.indices();
          IT_* lm_idx = this->_loc_mir.indices();
          for(Index i = 0; i < num_idx; ++i)
          {
            bm_idx[i] = aux_data_mir[i].first;
            lm_idx[i] = aux_data_mir[i].second;
          }

          return num_owned_dofs * Index(bh_);
        }

        template<typename DTA_>
        void upload_owned_data(DTA_* adp_val, const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& local_matrix) const
        {
          const DT_* loc_val = local_matrix.template val<LAFEM::Perspective::pod>();
          const IT_* adp_idx = _buf_mir.indices();
          const IT_* loc_idx = _loc_mir.indices();
          XASSERT(_buf_mir.num_indices() == _loc_mir.num_indices());
          Index n = _buf_mir.num_indices();
          for(Index i = 0; i < n; ++i)
          {
            ASSERT(loc_idx[i] < local_matrix.template used_elements<LAFEM::Perspective::pod>());
            adp_val[adp_idx[i]] = DTA_(loc_val[loc_idx[i]]);
          }
        }

        void gather_owner_data(
          LAFEM::DenseVector<DT_, IT_>& buffer,
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& local_matrix) const
        {
          DT_* buf_val = buffer.elements();
          const DT_* loc_val = local_matrix.template val<LAFEM::Perspective::pod>();
          const IT_* buf_idx = this->_buf_mir.indices();
          const IT_* loc_idx = this->_loc_mir.indices();
          XASSERT(this->_buf_mir.num_indices() == this->_loc_mir.num_indices());
          Index n = this->_buf_mir.num_indices();
          for(Index i = 0; i < n; ++i)
          {
            ASSERT(buf_idx[i] < buffer.size());
            ASSERT(loc_idx[i] < local_matrix.template used_elements<LAFEM::Perspective::pod>());
            buf_val[buf_idx[i]] = loc_val[loc_idx[i]];
          }
        }

      }; // class ADPDOMM<LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstMatrix_, typename... RestMatrix_>
      class ADPDOMM<LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>>
      {
      public:
        typedef typename FirstMatrix_::DataType DT_;
        typedef typename FirstMatrix_::IndexType IT_;

        ADPDOMM<FirstMatrix_> first;
        ADPDOMM<LAFEM::TupleMatrixRow<RestMatrix_...>> rest;

        Index num_indices() const
        {
          return first.num_indices() + rest.num_indices();
        }

        Index buf_size() const
        {
          return first.buf_size() + rest.buf_size();
        }

        Index loc_size() const
        {
          return first.loc_size() + rest.loc_size();
        }

        String dump_raw() const
        {
          return first.dump() + "," + rest.dump_raw();
        }

        String dump() const
        {
          return "[" + dump_raw() + "]";
        }

        template<
          typename RowMirror_,
          typename FirstIdxVec_, typename... RestIdxVec_>
        Index asm_neighbor_owner_data_mir(
          const LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>& local_matrix,
          const RowMirror_& row_mirror,
          const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buffer,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          Index row_offset)
        {
          Index n1 = first.asm_neighbor_owner_data_mir(local_matrix.first(), row_mirror, buffer, global_dof_idx.first(), row_offset);
          Index n2 = rest.asm_neighbor_owner_data_mir(local_matrix.rest(), row_mirror, buffer, global_dof_idx.rest(), row_offset);
          XASSERT(n1 == n2);
          return n1;
        }

        template<
          typename RowMirror_,
          typename FirstIdxVec_, typename... RestIdxVec_>
        Index asm_owned_data_mir(
          const LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>& local_matrix,
          const Adjacency::Graph& owned_graph,
          const RowMirror_& row_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          Index row_offset)
        {
          Index n1 = first.asm_owned_data_mir(local_matrix.first(), owned_graph, row_mirror, global_dof_idx.first(), row_offset);
          Index n2 = rest.asm_owned_data_mir(local_matrix.rest(), owned_graph, row_mirror, global_dof_idx.rest(), row_offset);
          XASSERT(n1 == n2);
          return n1;
        }

        template<typename DTA_>
        void upload_owned_data(DTA_* adp_val, const LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>& local_matrix) const
        {
          first.upload_owned_data(adp_val, local_matrix.first());
          rest.upload_owned_data(adp_val, local_matrix.rest());
        }

        void gather_owner_data(
          LAFEM::DenseVector<DT_, IT_>& buffer,
          const LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>& local_matrix) const
        {
          first.gather_owner_data(buffer, local_matrix.first());
          rest.gather_owner_data(buffer, local_matrix.rest());
        }
      }; // class ADPDOMM<LAFEM::TupleMatrixRow<First_, Rest_...>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstMatrix_>
      class ADPDOMM<LAFEM::TupleMatrixRow<FirstMatrix_>>
      {
      public:
        typedef typename FirstMatrix_::DataType DT_;
        typedef typename FirstMatrix_::IndexType IT_;

        ADPDOMM<FirstMatrix_> first;

        Index num_indices() const
        {
          return first.num_indices();
        }

        Index buf_size() const
        {
          return first.buf_size();
        }

        Index loc_size() const
        {
          return first.loc_size();
        }

        String dump_raw() const
        {
          return first.dump();
        }

        String dump() const
        {
          return "[" + dump_raw() + "]";
        }

        template<
          typename RowMirror_,
          typename FirstIdxVec_>
        Index asm_neighbor_owner_data_mir(
          const LAFEM::TupleMatrixRow<FirstMatrix_>& local_matrix,
          const RowMirror_& row_mirror,
          const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buffer,
          const LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx,
          Index row_offset)
        {
          return first.asm_neighbor_owner_data_mir(local_matrix.first(), row_mirror, buffer, global_dof_idx.first(), row_offset);
        }

        template<
          typename RowMirror_,
          typename FirstIdxVec_>
        Index asm_owned_data_mir(
          const LAFEM::TupleMatrixRow<FirstMatrix_>& local_matrix,
          const Adjacency::Graph& owned_graph,
          const RowMirror_& row_mirror,
          const LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx,
          Index row_offset)
        {
          return first.asm_owned_data_mir(local_matrix.first(), owned_graph, row_mirror, global_dof_idx.first(), row_offset);
        }

        template<typename DTA_>
        void upload_owned_data(DTA_* adp_val, const LAFEM::TupleMatrixRow<FirstMatrix_>& local_matrix) const
        {
          first.upload_owned_data(adp_val, local_matrix.first());
        }

        void gather_owner_data(
          LAFEM::DenseVector<DT_, IT_>& buffer,
          const LAFEM::TupleMatrixRow<FirstMatrix_>& local_matrix) const
        {
          first.gather_owner_data(buffer, local_matrix.first());
        }
      }; // class ADPDOMM<LAFEM::TupleMatrixRow<First_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstMatRow_, typename... RestMatRow_>
      class ADPDOMM<LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>>
      {
      public:
        typedef typename FirstMatRow_::DataType DT_;
        typedef typename FirstMatRow_::IndexType IT_;

        ADPDOMM<FirstMatRow_> first;
        ADPDOMM<LAFEM::TupleMatrix<RestMatRow_...>> rest;

        Index num_indices() const
        {
          return first.num_indices() + rest.num_indices();
        }

        Index buf_size() const
        {
          return first.buf_size() + rest.buf_size();
        }

        Index loc_size() const
        {
          return first.loc_size() + rest.loc_size();
        }

        String dump_raw() const
        {
          return first.dump() + "," + rest.dump_raw();
        }

        String dump() const
        {
          return "[" + dump_raw() + "]";
        }

        template<
          typename FirstRowMir_, typename... RestRowMir_,
          typename FirstIdxVec_, typename... RestIdxVec_>
        Index asm_neighbor_owner_data_mir(
          const LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>& local_matrix,
          const LAFEM::TupleMirror<FirstRowMir_, RestRowMir_...>& row_mirror,
          const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buffer,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          Index row_offset)
        {
          Index n1 = first.asm_neighbor_owner_data_mir(local_matrix.first(), row_mirror.first(), buffer, global_dof_idx, row_offset);
          Index n2 = rest.asm_neighbor_owner_data_mir(local_matrix.rest(), row_mirror.rest(), buffer, global_dof_idx, row_offset + n1);
          return n1 + n2;
        }

        template<
          typename FirstRowMir_, typename... RestRowMir_,
          typename FirstIdxVec_, typename... RestIdxVec_>
        Index asm_owned_data_mir(
          const LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>& local_matrix,
          const Adjacency::Graph& owned_graph,
          const LAFEM::TupleMirror<FirstRowMir_, RestRowMir_...>& row_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          Index row_offset)
        {
          Index n1 = first.asm_owned_data_mir(local_matrix.first(), owned_graph, row_mirror.first(), global_dof_idx, row_offset);
          Index n2 = rest.asm_owned_data_mir(local_matrix.rest(), owned_graph, row_mirror.rest(), global_dof_idx, row_offset + n1);
          return n1 + n2;
        }

        template<typename DTA_>
        void upload_owned_data(DTA_* adp_val, const LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>& local_matrix) const
        {
          first.upload_owned_data(adp_val, local_matrix.first());
          rest.upload_owned_data(adp_val, local_matrix.rest());
        }

        void gather_owner_data(
          LAFEM::DenseVector<DT_, IT_>& buffer,
          const LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>& local_matrix) const
        {
          first.gather_owner_data(buffer, local_matrix.first());
          rest.gather_owner_data(buffer, local_matrix.rest());
        }
      }; // class ADPDOMM<LAFEM::TupleMatrix<First_, Rest_...>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstMatRow_>
      class ADPDOMM<LAFEM::TupleMatrix<FirstMatRow_>>
      {
      public:
        typedef typename FirstMatRow_::DataType DT_;
        typedef typename FirstMatRow_::IndexType IT_;

        ADPDOMM<FirstMatRow_> first;

        Index num_indices() const
        {
          return first.num_indices();
        }

        Index buf_size() const
        {
          return first.buf_size();
        }

        Index loc_size() const
        {
          return first.loc_size();
        }

        String dump_raw() const
        {
          return first.dump();
        }

        String dump() const
        {
          return "[" + dump_raw() + "]";
        }

        template<
          typename FirstRowMir_,
          typename FirstIdxVec_, typename... RestIdxVec_>
        Index asm_neighbor_owner_data_mir(
          const LAFEM::TupleMatrix<FirstMatRow_>& local_matrix,
          const LAFEM::TupleMirror<FirstRowMir_>& row_mirror,
          const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buffer,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          Index row_offset)
        {
          return first.asm_neighbor_owner_data_mir(local_matrix.first(), row_mirror.first(), buffer, global_dof_idx, row_offset);
        }

        template<
          typename FirstRowMir_,
          typename FirstIdxVec_, typename... RestIdxVec_>
        Index asm_owned_data_mir(
          const LAFEM::TupleMatrix<FirstMatRow_>& local_matrix,
          const Adjacency::Graph& owned_graph,
          const LAFEM::TupleMirror<FirstRowMir_>& row_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          Index row_offset)
        {
          return first.asm_owned_data_mir(local_matrix.first(), owned_graph, row_mirror.first(), global_dof_idx, row_offset);
        }

        template<typename DTA_>
        void upload_owned_data(DTA_* adp_val, const LAFEM::TupleMatrix<FirstMatRow_>& local_matrix) const
        {
          first.upload_owned_data(adp_val, local_matrix.first());
        }

        void gather_owner_data(
          LAFEM::DenseVector<DT_, IT_>& buffer,
          const LAFEM::TupleMatrix<FirstMatRow_>& local_matrix) const
        {
          first.gather_owner_data(buffer, local_matrix.first());
        }
      }; // class ADPDOMM<LAFEM::TupleMatrix<First_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
      class ADPDOMM<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>>
      {
      public:
        typedef typename MatrixA_::DataType DT_;
        typedef typename MatrixA_::IndexType IT_;

        ADPDOMM<MatrixA_> block_a;
        ADPDOMM<MatrixB_> block_b;
        ADPDOMM<MatrixD_> block_d;

        Index num_indices() const
        {
          return block_a.num_indices() + block_b.num_indices() + block_d.num_indices();
        }

        Index buf_size() const
        {
          return block_a.buf_size() + block_b.buf_size() + block_d.buf_size();
        }

        Index loc_size() const
        {
          return block_a.loc_size() + block_b.loc_size() + block_d.loc_size();
        }

        String dump() const
        {
          return "[[" + block_a.dump() + "," + block_b.dump() + "],[" + block_d.dump() + ",0]";
        }

        template<
          typename MirrorV_, typename MirrorP_,
          typename IdxVecV_, typename IdxVecP_>
        Index asm_neighbor_owner_data_mir(
          const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& local_matrix,
          const LAFEM::TupleMirror<MirrorV_, MirrorP_>& row_mirror,
          const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buffer,
          const LAFEM::TupleVector<IdxVecV_, IdxVecP_>& global_dof_idx,
          Index row_offset)
        {
          Index na = block_a.asm_neighbor_owner_data_mir(local_matrix.block_a(), row_mirror.template at<0>(), buffer, global_dof_idx.template at<0>(), row_offset);
          Index nb = block_b.asm_neighbor_owner_data_mir(local_matrix.block_b(), row_mirror.template at<0>(), buffer, global_dof_idx.template at<1>(), row_offset);
          XASSERT(na == nb);
          Index nd = block_d.asm_neighbor_owner_data_mir(local_matrix.block_d(), row_mirror.template at<1>(), buffer, global_dof_idx.template at<0>(), row_offset + na);
          return na + nd;
        }

        template<
          typename MirrorV_, typename MirrorP_,
          typename IdxVecV_, typename IdxVecP_>
        Index asm_owned_data_mir(
          const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& local_matrix,
          const Adjacency::Graph& owned_graph,
          const LAFEM::TupleMirror<MirrorV_, MirrorP_>& row_mirror,
          const LAFEM::TupleVector<IdxVecV_, IdxVecP_>& global_dof_idx,
          Index row_offset)
        {
          Index na = block_a.asm_owned_data_mir(local_matrix.block_a(), owned_graph, row_mirror.template at<0>(), global_dof_idx.template at<0>(), row_offset);
          Index nb = block_b.asm_owned_data_mir(local_matrix.block_b(), owned_graph, row_mirror.template at<0>(), global_dof_idx.template at<1>(), row_offset);
          XASSERT(na == nb);
          Index nd = block_d.asm_owned_data_mir(local_matrix.block_d(), owned_graph, row_mirror.template at<1>(), global_dof_idx.template at<0>(), row_offset + na);
          return na + nd;
        }

        template<typename DTA_>
        void upload_owned_data(DTA_* adp_val, const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& local_matrix) const
        {
          block_a.upload_owned_data(adp_val, local_matrix.block_a());
          block_b.upload_owned_data(adp_val, local_matrix.block_b());
          block_d.upload_owned_data(adp_val, local_matrix.block_d());
        }

        void gather_owner_data(
          LAFEM::DenseVector<DT_, IT_>& buffer,
          const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& local_matrix) const
        {
          block_a.gather_owner_data(buffer, local_matrix.block_a());
          block_b.gather_owner_data(buffer, local_matrix.block_b());
          block_d.gather_owner_data(buffer, local_matrix.block_d());
        }
      }; // class ADPDOMM<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_>
      class ADPMatAux<LAFEM::SparseMatrixCSR<DT_, IT_>>
      {
      public:
        static Index calc_mat_buf_num_rows(
          const LAFEM::SparseMatrixCSR<DT_, IT_>& DOXY(local_matrix),
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror)
        {
          return row_mirror.num_indices();
        }

        static Index calc_mat_buf_row_nze(
          IT_* buf_row_nze,
          const LAFEM::SparseMatrixCSR<DT_, IT_>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror)
        {
          const Index num_idx = row_mirror.num_indices();
          const IT_* mir_idx = row_mirror.indices();
          const IT_* row_ptr = local_matrix.row_ptr();

          for(Index i(0); i < num_idx; ++i)
          {
            buf_row_nze[i] += (row_ptr[mir_idx[i]+1] - row_ptr[mir_idx[i]]);
          }

          return num_idx;
        }

        static Index gather_mat_buf_col_idx(
          IT_* aux_row_ptr,
          IT_* buf_col_idx,
          const LAFEM::SparseMatrixCSR<DT_, IT_>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVector<IT_, IT_>& global_dof_idx)
        {
          const Index num_idx = row_mirror.num_indices();
          const IT_* mir_idx = row_mirror.indices();
          const IT_* row_ptr = local_matrix.row_ptr();
          const IT_* col_idx = local_matrix.col_ind();
          const IT_* glob_dof_idx = global_dof_idx.elements();

          for(Index i(0); i < num_idx; ++i)
          {
            // k must be a reference as we need to update the auxiliary row pointer
            IT_& k = aux_row_ptr[i];
            const IT_ row = mir_idx[i];

            // map local DOFs to global DOFs
            for(IT_ j(row_ptr[row]); j < row_ptr[row+1]; ++j, ++k)
              buf_col_idx[k] = glob_dof_idx[col_idx[j]];
          }

          return num_idx;
        }

        static Index gather_owned_struct(
          Adjacency::DynamicGraph& graph,
          const LAFEM::SparseMatrixCSR<DT_, IT_>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVector<IT_, IT_>& global_dof_idx,
          const Index row_offset)
        {
          const Index num_rows = row_mirror.num_indices();
          const IT_* loc_row_ptr = local_matrix.row_ptr();
          const IT_* loc_col_idx = local_matrix.col_ind();
          const IT_* row_mir_idx = row_mirror.indices();
          const IT_* glob_dof_idx = global_dof_idx.elements();

          for(Index i(0); i < num_rows; ++i)
          {
            const IT_ lrow = row_mir_idx[i];
            for(IT_ j(loc_row_ptr[lrow]); j < loc_row_ptr[lrow + 1]; ++j)
            {
              graph.insert(row_offset + i, Index(glob_dof_idx[loc_col_idx[j]]));
            }
          }
          return num_rows;
        }
      }; // class ADPMatAux<LAFEM::SparseMatrixCSR<DT_, IT_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_, int bh_, int bw_>
      class ADPMatAux<LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>>
      {
      public:
        static Index calc_mat_buf_num_rows(
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& DOXY(local_matrix),
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror)
        {
          return row_mirror.num_indices() * Index(bh_);
        }

        static Index calc_mat_buf_row_nze(
          IT_* buf_row_nze,
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror)
        {
          const Index num_idx = row_mirror.num_indices();
          const IT_* mir_idx = row_mirror.indices();
          const IT_* row_ptr = local_matrix.row_ptr();

          for(Index i(0), j(0); i < num_idx; ++i)
          {
            IT_ nn = IT_(bw_) * (row_ptr[mir_idx[i]+1] - row_ptr[mir_idx[i]]);
            for(int k(0); k < bh_; ++k, ++j)
              buf_row_nze[j] += nn;
          }

          return num_idx * Index(bh_);
        }

        static Index gather_mat_buf_col_idx(
          IT_* aux_row_ptr,
          IT_* buf_col_idx,
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bw_>& global_dof_idx)
        {
          const Index num_idx = row_mirror.num_indices();
          const IT_* mir_idx = row_mirror.indices();
          const IT_* row_ptr = local_matrix.row_ptr();
          const IT_* col_idx = local_matrix.col_ind();
          const auto* glob_dof_idx = global_dof_idx.elements();

          for(Index i(0); i < num_idx; ++i)
          {
            const IT_ row = mir_idx[i];

            // loop over the rows of the block
            for(int bi(0); bi < bh_; ++bi)
            {
              // k must be a reference as we need to update the auxiliary row pointer
              IT_& k = aux_row_ptr[i*IT_(bh_) + IT_(bi)];

              // map local DOFs to global DOFs
              for(IT_ j(row_ptr[row]); j < row_ptr[row+1]; ++j)
              {
                // loop over the columns of this block
                for(int bj(0); bj < bw_; ++bj, ++k)
                {
                  buf_col_idx[k] = glob_dof_idx[col_idx[j]][bj];
                }
              }
            }
          }

          return num_idx * Index(bh_);
        }

        static Index gather_mat_buf_col_idx(
          IT_* aux_row_ptr,
          IT_* buf_col_idx,
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, 1>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVector<IT_, IT_>& global_dof_idx)
        {
          const Index num_idx = row_mirror.num_indices();
          const IT_* mir_idx = row_mirror.indices();
          const IT_* row_ptr = local_matrix.row_ptr();
          const IT_* col_idx = local_matrix.col_ind();
          const auto* glob_dof_idx = global_dof_idx.elements();

          for(Index i(0); i < num_idx; ++i)
          {
            const IT_ row = mir_idx[i];

            // loop over the rows of the block
            for(int bi(0); bi < bh_; ++bi)
            {
              // k must be a reference as we need to update the auxiliary row pointer
              IT_& k = aux_row_ptr[i*IT_(bh_) + IT_(bi)];

              // map local DOFs to global DOFs
              for(IT_ j(row_ptr[row]); j < row_ptr[row+1]; ++j, ++k)
              {
                buf_col_idx[k] = glob_dof_idx[col_idx[j]];
              }
            }
          }

          return num_idx * Index(bh_);
        }

        static Index gather_owned_struct(
          Adjacency::DynamicGraph& graph,
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bw_>& global_dof_idx,
          const Index row_offset)
        {
          const Index num_rows = row_mirror.num_indices();
          const IT_* loc_row_ptr = local_matrix.row_ptr();
          const IT_* loc_col_idx = local_matrix.col_ind();
          const IT_* row_mir_idx = row_mirror.indices();
          const auto* glob_dof_idx = global_dof_idx.elements();

          for(Index i(0); i < num_rows; ++i)
          {
            const Index lrow = row_mir_idx[i];
            for(IT_ j(loc_row_ptr[lrow]); j < loc_row_ptr[lrow + 1]; ++j)
            {
              const auto& gci = glob_dof_idx[loc_col_idx[j]];
              for(int bi = 0; bi < bh_; ++bi)
              {
                for(int bj = 0; bj < bw_; ++bj)
                {
                  graph.insert(row_offset + i*Index(bh_) + Index(bi), Index(gci[bj]));
                }
              }
            }
          }
          return num_rows * Index(bh_);
        }

        static Index gather_owned_struct(
          Adjacency::DynamicGraph& graph,
          const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, 1>& local_matrix,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          const LAFEM::DenseVector<IT_, IT_>& global_dof_idx,
          const Index row_offset)
        {
          const Index num_rows = row_mirror.num_indices();
          const IT_* loc_row_ptr = local_matrix.row_ptr();
          const IT_* loc_col_idx = local_matrix.col_ind();
          const IT_* row_mir_idx = row_mirror.indices();
          const IT_* glob_dof_idx = global_dof_idx.elements();

          for(Index i(0); i < num_rows; ++i)
          {
            const Index lrow = row_mir_idx[i];
            for(IT_ j(loc_row_ptr[lrow]); j < loc_row_ptr[lrow + 1]; ++j)
            {
              const IT_ gci = glob_dof_idx[loc_col_idx[j]];
              for(int bi = 0; bi < bh_; ++bi)
              {
                graph.insert(row_offset + i*Index(bh_) + Index(bi), Index(gci));
              }
            }
          }
          return num_rows * Index(bh_);
        }
      }; // class ADPMatAux<LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstMatrix_, typename... RestMatrix_>
      class ADPMatAux<LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>>
      {
      public:
        typedef typename FirstMatrix_::IndexType IT_;
        typedef ADPMatAux<FirstMatrix_> ADPMatAuxFirst;
        typedef ADPMatAux<LAFEM::TupleMatrixRow<RestMatrix_...>> ADPMatAuxRest;

        template<typename RowMirror_>
        static Index calc_mat_buf_num_rows(
          const LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>& local_matrix,
          const RowMirror_& row_mirror)
        {
          Index n1 = ADPMatAuxFirst::calc_mat_buf_num_rows(local_matrix.first(), row_mirror);
          Index n2 = ADPMatAuxRest::calc_mat_buf_num_rows(local_matrix.rest(), row_mirror);
          XASSERT(n1 == n2);
          return n1;
        }

        template<typename RowMirror_>
        static Index calc_mat_buf_row_nze(
          IT_* buf_row_nze,
          const LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>& local_matrix,
          const RowMirror_& row_mirror)
        {
          Index n1 = ADPMatAuxFirst::calc_mat_buf_row_nze(buf_row_nze, local_matrix.first(), row_mirror);
          Index n2 = ADPMatAuxRest::calc_mat_buf_row_nze(buf_row_nze, local_matrix.rest(), row_mirror);
          XASSERT(n1 == n2);
          return n1;
        }

        template<
          typename RowMirror_,
          typename FirstIdxVec_, typename... RestIdxVec_>
        static Index gather_mat_buf_col_idx(
          IT_* aux_row_ptr,
          IT_* buf_col_idx,
          const LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>& local_matrix,
          const RowMirror_& row_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx)
        {
          Index n1 = ADPMatAuxFirst::gather_mat_buf_col_idx(aux_row_ptr, buf_col_idx, local_matrix.first(), row_mirror, global_dof_idx.first());
          Index n2 = ADPMatAuxRest::gather_mat_buf_col_idx(aux_row_ptr, buf_col_idx, local_matrix.rest(), row_mirror, global_dof_idx.rest());
          XASSERT(n1 == n2);
          return n1;
        }

        template<
          typename RowMirror_,
          typename FirstIdxVec_, typename... RestIdxVec_>
        static Index gather_owned_struct(
          Adjacency::DynamicGraph& graph,
          const LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>& local_matrix,
          const RowMirror_& row_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          const Index row_offset)
        {
          Index n1 = ADPMatAuxFirst::gather_owned_struct(graph, local_matrix.first(), row_mirror, global_dof_idx.first(), row_offset);
          Index n2 = ADPMatAuxRest::gather_owned_struct(graph, local_matrix.rest(), row_mirror, global_dof_idx.rest(), row_offset);
          XASSERT(n1 == n2);
          return n1;
        }
      }; // class ADPMatAux<LAFEM::TupleMatrixRow<FirstMatrix_, RestMatrix_...>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstMatrix_>
      class ADPMatAux<LAFEM::TupleMatrixRow<FirstMatrix_>>
      {
      public:
        typedef typename FirstMatrix_::IndexType IT_;
        typedef ADPMatAux<FirstMatrix_> ADPMatAuxFirst;

        template<typename RowMirror_>
        static Index calc_mat_buf_num_rows(
          const LAFEM::TupleMatrixRow<FirstMatrix_>& local_matrix,
          const RowMirror_& row_mirror)
        {
          return ADPMatAuxFirst::calc_mat_buf_num_rows(local_matrix.first(), row_mirror);
        }

        template<typename RowMirror_>
        static Index calc_mat_buf_row_nze(
          IT_* buf_row_nze,
          const LAFEM::TupleMatrixRow<FirstMatrix_>& local_matrix,
          const RowMirror_& row_mirror)
        {
          return ADPMatAuxFirst::calc_mat_buf_row_nze(buf_row_nze, local_matrix.first(), row_mirror);
        }

        template<
          typename RowMirror_,
          typename FirstIdxVec_>
        static Index gather_mat_buf_col_idx(
          IT_* aux_row_ptr,
          IT_* buf_col_idx,
          const LAFEM::TupleMatrixRow<FirstMatrix_>& local_matrix,
          const RowMirror_& row_mirror,
          const LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx)
        {
          return ADPMatAuxFirst::gather_mat_buf_col_idx(aux_row_ptr, buf_col_idx, local_matrix.first(), row_mirror, global_dof_idx.first());
        }

        template<
          typename RowMirror_,
          typename FirstIdxVec_>
        static Index gather_owned_struct(
          Adjacency::DynamicGraph& graph,
          const LAFEM::TupleMatrixRow<FirstMatrix_>& local_matrix,
          const RowMirror_& row_mirror,
          const LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx,
          const Index row_offset)
        {
          return ADPMatAuxFirst::gather_owned_struct(graph, local_matrix.first(), row_mirror, global_dof_idx.first(), row_offset);
        }
      }; // class ADPMatAux<LAFEM::TupleMatrixRow<FirstMatrix_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstMatRow_, typename... RestMatRow_>
      class ADPMatAux<LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>>
      {
      public:
        typedef typename FirstMatRow_::IndexType IT_;
        typedef ADPMatAux<FirstMatRow_> ADPMatAuxFirst;
        typedef ADPMatAux<LAFEM::TupleMatrix<RestMatRow_...>> ADPMatAuxRest;

        template<typename FirstMirror_, typename... RestMirror_>
        static Index calc_mat_buf_num_rows(
          const LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>& local_matrix,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& row_mirror)
        {
          Index n1 = ADPMatAuxFirst::calc_mat_buf_num_rows(local_matrix.first(), row_mirror.first());
          Index n2 = ADPMatAuxRest::calc_mat_buf_num_rows(local_matrix.rest(), row_mirror.rest());
          return n1 + n2;
        }

        template<typename FirstMirror_, typename... RestMirror_>
        static Index calc_mat_buf_row_nze(
          IT_* buf_row_nze,
          const LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>& local_matrix,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& row_mirror)
        {
          Index n1 = ADPMatAuxFirst::calc_mat_buf_row_nze(buf_row_nze, local_matrix.first(), row_mirror.first());
          Index n2 = ADPMatAuxRest::calc_mat_buf_row_nze(&buf_row_nze[n1], local_matrix.rest(), row_mirror.rest());
          return n1 + n2;
        }

        template<
          typename FirstMirror_, typename... RestMirror_,
          typename IndexVector_>
        static Index gather_mat_buf_col_idx(
          IT_* aux_row_ptr,
          IT_* buf_col_idx,
          const LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>& local_matrix,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& row_mirror,
          const IndexVector_& global_dof_idx)
        {
          Index n1 = ADPMatAuxFirst::gather_mat_buf_col_idx(aux_row_ptr, buf_col_idx, local_matrix.first(), row_mirror.first(), global_dof_idx);
          Index n2 = ADPMatAuxRest::gather_mat_buf_col_idx(&aux_row_ptr[n1], buf_col_idx, local_matrix.rest(), row_mirror.rest(), global_dof_idx);
          return n1 + n2;
        }

        template<
          typename FirstMirror_, typename... RestMirror_,
          typename IndexVector_>
        static Index gather_owned_struct(
          Adjacency::DynamicGraph& graph,
          const LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>& local_matrix,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& row_mirror,
          const IndexVector_& global_dof_idx,
          const Index row_offset)
        {
          Index n1 = ADPMatAuxFirst::gather_owned_struct(graph, local_matrix.first(), row_mirror.first(), global_dof_idx, row_offset);
          Index n2 = ADPMatAuxRest::gather_owned_struct(graph, local_matrix.rest(), row_mirror.rest(), global_dof_idx, row_offset + n1);
          return n1 + n2;
        }
      }; // class ADPMatAux<LAFEM::TupleMatrix<FirstMatRow_, RestMatRow_...>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstMatRow_>
      class ADPMatAux<LAFEM::TupleMatrix<FirstMatRow_>>
      {
      public:
        typedef typename FirstMatRow_::IndexType IT_;
        typedef ADPMatAux<FirstMatRow_> ADPMatAuxFirst;

        template<typename FirstMirror_>
        static Index calc_mat_buf_num_rows(
          const LAFEM::TupleMatrix<FirstMatRow_>& local_matrix,
          const LAFEM::TupleMirror<FirstMirror_>& row_mirror)
        {
          return ADPMatAuxFirst::calc_mat_buf_num_rows(local_matrix.first(), row_mirror.first());
        }

        template<typename FirstMirror_>
        static Index calc_mat_buf_row_nze(
          IT_* buf_row_nze,
          const LAFEM::TupleMatrix<FirstMatRow_>& local_matrix,
          const LAFEM::TupleMirror<FirstMirror_>& row_mirror)
        {
          return ADPMatAuxFirst::calc_mat_buf_row_nze(buf_row_nze, local_matrix.first(), row_mirror.first());
        }

        template<
          typename FirstMirror_,
          typename IndexVector_>
        static Index gather_mat_buf_col_idx(
          IT_* aux_row_ptr,
          IT_* buf_col_idx,
          const LAFEM::TupleMatrix<FirstMatRow_>& local_matrix,
          const LAFEM::TupleMirror<FirstMirror_>& row_mirror,
          const IndexVector_& global_dof_idx)
        {
          return ADPMatAuxFirst::gather_mat_buf_col_idx(aux_row_ptr, buf_col_idx, local_matrix.first(), row_mirror.first(), global_dof_idx);
        }

        template<
          typename FirstMirror_,
          typename IndexVector_>
        static Index gather_owned_struct(
          Adjacency::DynamicGraph& graph,
          const LAFEM::TupleMatrix<FirstMatRow_>& local_matrix,
          const LAFEM::TupleMirror<FirstMirror_>& row_mirror,
          const IndexVector_& global_dof_idx,
          const Index row_offset)
        {
          return ADPMatAuxFirst::gather_owned_struct(graph, local_matrix.first(), row_mirror.first(), global_dof_idx, row_offset);
        }
      }; // class ADPMatAux<LAFEM::TupleMatrix<FirstMatRow_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename MatrixA_, typename MatrixB_, typename MatrixD_>
      class ADPMatAux<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>>
      {
      public:
        typedef typename MatrixA_::IndexType IT_;
        typedef ADPMatAux<MatrixA_> ADPMatAuxA;
        typedef ADPMatAux<MatrixB_> ADPMatAuxB;
        typedef ADPMatAux<MatrixD_> ADPMatAuxD;

        template<typename MirrorV_, typename MirrorP_>
        static Index calc_mat_buf_num_rows(
          const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& local_matrix,
          const LAFEM::TupleMirror<MirrorV_, MirrorP_>& row_mirror)
        {
          Index n1 = ADPMatAuxB::calc_mat_buf_num_rows(local_matrix.block_b(), row_mirror.template at<0>());
          Index n2 = ADPMatAuxD::calc_mat_buf_num_rows(local_matrix.block_d(), row_mirror.template at<1>());
          return n1 + n2;
        }

        template<typename MirrorV_, typename MirrorP_>
        static Index calc_mat_buf_row_nze(
          IT_* buf_row_nze,
          const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& local_matrix,
          const LAFEM::TupleMirror<MirrorV_, MirrorP_>& row_mirror)
        {
          Index na = ADPMatAuxA::calc_mat_buf_row_nze(buf_row_nze, local_matrix.block_a(), row_mirror.template at<0>());
          Index nb = ADPMatAuxB::calc_mat_buf_row_nze(buf_row_nze, local_matrix.block_b(), row_mirror.template at<0>());
          XASSERT(na == nb);
          Index nd = ADPMatAuxD::calc_mat_buf_row_nze(&buf_row_nze[na], local_matrix.block_d(), row_mirror.template at<1>());
          return na + nd;
        }

        template<
          typename MirrorV_, typename MirrorP_,
          typename IdxVecV_, typename IdxVecP_>
        static Index gather_mat_buf_col_idx(
          IT_* aux_row_ptr,
          IT_* buf_col_idx,
          const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& local_matrix,
          const LAFEM::TupleMirror<MirrorV_, MirrorP_>& row_mirror,
          const LAFEM::TupleVector<IdxVecV_, IdxVecP_>& global_dof_idx)
        {
          Index na = ADPMatAuxA::gather_mat_buf_col_idx(aux_row_ptr, buf_col_idx, local_matrix.block_a(), row_mirror.template at<0>(), global_dof_idx.template at<0>());
          Index nb = ADPMatAuxB::gather_mat_buf_col_idx(aux_row_ptr, buf_col_idx, local_matrix.block_b(), row_mirror.template at<0>(), global_dof_idx.template at<1>());
          XASSERT(na == nb);
          Index nd = ADPMatAuxD::gather_mat_buf_col_idx(&aux_row_ptr[na], buf_col_idx, local_matrix.block_d(), row_mirror.template at<1>(), global_dof_idx.template at<0>());
          return na + nd;
        }

        template<
          typename MirrorV_, typename MirrorP_,
          typename IdxVecV_, typename IdxVecP_>
        static Index gather_owned_struct(
          Adjacency::DynamicGraph& graph,
          const LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>& local_matrix,
          const LAFEM::TupleMirror<MirrorV_, MirrorP_>& row_mirror,
          const LAFEM::TupleVector<IdxVecV_, IdxVecP_>& global_dof_idx,
          const Index row_offset)
        {
          Index na = ADPMatAuxA::gather_owned_struct(graph, local_matrix.block_a(), row_mirror.template at<0>(), global_dof_idx.template at<0>(), row_offset);
          Index nb = ADPMatAuxB::gather_owned_struct(graph, local_matrix.block_b(), row_mirror.template at<0>(), global_dof_idx.template at<1>(), row_offset);
          XASSERT(na == nb);
          Index nd = ADPMatAuxD::gather_owned_struct(graph, local_matrix.block_d(), row_mirror.template at<1>(), global_dof_idx.template at<0>(), row_offset + na);
          return na + nd;
        }
      }; // class ADPMatAux<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_>
      class ADPFilAux<LAFEM::NoneFilter<DT_, IT_>>
      {
      public:
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& DOXY(unit_rows),
          const LAFEM::NoneFilter<DT_, IT_>& DOXY(filter),
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          Index DOXY(row_offset))
        {
          // nothing to do here
          return row_mirror.num_indices();
        }
      }; // class ADPFilAux<LAFEM::NoneFilter<DT_, IT_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_, int bs_>
      class ADPFilAux<LAFEM::NoneFilterBlocked<DT_, IT_, bs_>>
      {
      public:
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& DOXY(unit_rows),
          const LAFEM::NoneFilterBlocked<DT_, IT_, bs_>& DOXY(filter),
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          Index DOXY(row_offset))
        {
          // nothing to do here
          return row_mirror.num_indices() * Index(bs_);
        }
      }; // class ADPFilAux<LAFEM::NoneFilterBlocked<DT_, IT_, bs_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_>
      class ADPFilAux<LAFEM::UnitFilter<DT_, IT_>>
      {
      public:
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& unit_rows,
          const LAFEM::UnitFilter<DT_, IT_>& filter,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          Index row_offset)
        {
          if(filter.used_elements() <= Index(0))
            return row_mirror.num_indices();

          const IT_* fil_idx = filter.get_indices();
          const IT_* mir_idx = row_mirror.indices();
          const Index num_fil = filter.used_elements();
          const Index num_mir = row_mirror.num_indices();

          Index i = 0u, j = 0u;
          while((i < num_fil) && (j < num_mir))
          {
            if(fil_idx[i] < Index(mir_idx[j]))
              ++i;
            else if(fil_idx[i] > Index(mir_idx[j]))
              ++j;
            else
            {
              unit_rows.add(IT_(row_offset + j), DT_(0));
              ++i;
              ++j;
            }
          }

          return row_mirror.num_indices();
        }
      }; // class ADPFilAux<LAFEM::UnitFilter<DT_, IT_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_, int bs_>
      class ADPFilAux<LAFEM::UnitFilterBlocked<DT_, IT_, bs_>>
      {
      public:
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& unit_rows,
          const LAFEM::UnitFilterBlocked<DT_, IT_, bs_>& filter,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          Index row_offset)
        {
          if(filter.used_elements() <= Index(0))
            return row_mirror.num_indices() * Index(bs_);

          const bool ignore_nans = filter.get_ignore_nans();
          const auto* fil_val = filter.get_values();
          const IT_* fil_idx = filter.get_indices();
          const IT_* mir_idx = row_mirror.indices();
          const Index num_fil = filter.used_elements();
          const Index num_mir = row_mirror.num_indices();

          Index i = 0u, j = 0u;
          while((i < num_fil) && (j < num_mir))
          {
            if(fil_idx[i] < Index(mir_idx[j]))
              ++i;
            else if(fil_idx[i] > Index(mir_idx[j]))
              ++j;
            else
            {
              // process the entire block
              for(int k = 0; k < bs_; ++k)
              {
                if(!ignore_nans || !Math::isnan(fil_val[i][k]))
                  unit_rows.add(IT_(row_offset + j*Index(bs_) + Index(k)), DT_(0));
              }
              ++i;
            }
          }

          return row_mirror.num_indices() * Index(bs_);
        }
      }; // class ADPFilAux<LAFEM::UnitFilterBlocked<DT_, IT_, bs_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_>
      class ADPFilAux<LAFEM::MeanFilter<DT_, IT_>>
      {
      public:
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& DOXY(unit_rows),
          const LAFEM::MeanFilter<DT_, IT_>& filter,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          Index DOXY(row_offset))
        {
          XASSERTM(filter.empty(), "AlgDofParti does not support LAFEM::MeanFilter yet!");
          return row_mirror.num_indices();
        }
      }; // class ADPFilAux<LAFEM::MeanFilter<DT_, IT_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_, int bs_>
      class ADPFilAux<LAFEM::MeanFilterBlocked<DT_, IT_, bs_>>
      {
      public:
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& DOXY(unit_rows),
          const LAFEM::MeanFilterBlocked<DT_, IT_, bs_>& filter,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          Index DOXY(row_offset))
        {
          XASSERTM(filter.empty(), "AlgDofParti does not support LAFEM::MeanFilterBlocked yet!");
          return row_mirror.num_indices() * Index(bs_);
        }
      }; // class ADPFilAux<LAFEM::MeanFilterBlocked<DT_, IT_, bs_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_>
      class ADPFilAux<Global::MeanFilter<DT_, IT_>>
      {
      public:
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& DOXY(unit_rows),
          const Global::MeanFilter<DT_, IT_>& filter,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          Index DOXY(row_offset))
        {
          XASSERTM(filter.empty(), "AlgDofParti does not support Global::MeanFilter yet!");
          return row_mirror.num_indices();
        }
      }; // class ADPFilAux<Global::MeanFilter<DT_, IT_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename DT_, typename IT_, int bs_>
      class ADPFilAux<LAFEM::SlipFilter<DT_, IT_, bs_>>
      {
      public:
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& DOXY(unit_rows),
          const LAFEM::SlipFilter<DT_, IT_, bs_>& filter,
          const LAFEM::VectorMirror<DT_, IT_>& row_mirror,
          Index DOXY(row_offset))
        {
          XASSERTM(filter.empty(), "AlgDofParti does not support LAFEM::SlipFilter yet!");
          return row_mirror.num_indices() * Index(bs_);
        }
      }; // class ADPFilAux<LAFEM::SlipFilter<DT_, IT_, bs_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstFilter_, typename... RestFilter_>
      class ADPFilAux<LAFEM::FilterChain<FirstFilter_, RestFilter_...>>
      {
      public:
        typedef typename FirstFilter_::DataType DT_;
        typedef typename FirstFilter_::IndexType IT_;
        typedef ADPFilAux<FirstFilter_> ADPFilAuxFirst;
        typedef ADPFilAux<LAFEM::FilterChain<RestFilter_...>> ADPFilAuxRest;

        template<typename RowMirror_>
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& unit_rows,
          const LAFEM::FilterChain<FirstFilter_, RestFilter_...>& filter,
          const RowMirror_& row_mirror,
          Index row_offset)
        {
          Index n1 = ADPFilAuxFirst::upload_filter(unit_rows, filter.first(), row_mirror, row_offset);
          Index n2 = ADPFilAuxRest::upload_filter(unit_rows, filter.rest(), row_mirror, row_offset);
          XASSERT(n1 == n2);
          return n1;
        }
      }; // class ADPFilAux<LAFEM::FilterChain<FirstFilter_, RestFilter_...>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstFilter_>
      class ADPFilAux<LAFEM::FilterChain<FirstFilter_>>
      {
      public:
        typedef typename FirstFilter_::DataType DT_;
        typedef typename FirstFilter_::IndexType IT_;
        typedef ADPFilAux<FirstFilter_> ADPFilAuxFirst;

        template<typename RowMirror_>
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& unit_rows,
          const LAFEM::FilterChain<FirstFilter_>& filter,
          const RowMirror_& row_mirror,
          Index row_offset)
        {
          return ADPFilAuxFirst::upload_filter(unit_rows, filter.first(), row_mirror, row_offset);
        }
      }; // class ADPFilAux<LAFEM::FilterChain<FirstFilter_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename SubFilter_>
      class ADPFilAux<LAFEM::FilterSequence<SubFilter_>>
      {
      public:
        typedef typename SubFilter_::DataType DT_;
        typedef typename SubFilter_::IndexType IT_;
        typedef ADPFilAux<SubFilter_> ADPFilAuxSub;

        template<typename RowMirror_>
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& unit_rows,
          const LAFEM::FilterSequence<SubFilter_>& filter,
          const RowMirror_& row_mirror,
          Index row_offset)
        {
          for(const auto& x : filter)
          {
            ADPFilAuxSub::upload_filter(unit_rows, x.second, row_mirror, row_offset);
          }
          return row_mirror.num_indices();
        }
      }; // class ADPFilAux<LAFEM::FilterSequence<SubFilter_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstFilter_, typename... RestFilter_>
      class ADPFilAux<LAFEM::TupleFilter<FirstFilter_, RestFilter_...>>
      {
      public:
        typedef typename FirstFilter_::DataType DT_;
        typedef typename FirstFilter_::IndexType IT_;
        typedef ADPFilAux<FirstFilter_> ADPFilAuxFirst;
        typedef ADPFilAux<LAFEM::TupleFilter<RestFilter_...>> ADPFilAuxRest;

        template<typename FirstMirror_, typename... RestMirror_>
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& unit_rows,
          const LAFEM::TupleFilter<FirstFilter_, RestFilter_...>& filter,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& row_mirror,
          Index row_offset)
        {
          Index n1 = ADPFilAuxFirst::upload_filter(unit_rows, filter.first(), row_mirror.first(), row_offset);
          Index n2 = ADPFilAuxRest::upload_filter(unit_rows, filter.rest(), row_mirror.rest(), row_offset + n1);
          return n1 + n2;
        }
      }; // class ADPFilAux<LAFEM::TupleFilter<FirstFilter_, RestFilter_...>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstFilter_>
      class ADPFilAux<LAFEM::TupleFilter<FirstFilter_>>
      {
      public:
        typedef typename FirstFilter_::DataType DT_;
        typedef typename FirstFilter_::IndexType IT_;
        typedef ADPFilAux<FirstFilter_> ADPFilAuxFirst;

        template<typename FirstMirror_>
        static Index upload_filter(
          LAFEM::UnitFilter<DT_, IT_>& unit_rows,
          const LAFEM::TupleFilter<FirstFilter_>& filter,
          const LAFEM::TupleMirror<FirstMirror_>& row_mirror,
          Index row_offset)
        {
          return ADPFilAuxFirst::upload_filter(unit_rows, filter.first(), row_mirror.first(), row_offset);
        }
      }; // class ADPFilAux<LAFEM::TupleFilter<FirstFilter_, RestFilter_...>>
    } // namespace Intern
    /// \endcond
  } // namespace Global
} // namespace FEAT
