// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/arch/axpy_dense.hpp>
#include <kernel/lafem/arch/max_rel_diff_dense.hpp>
#include <kernel/lafem/arch/norm2_dense.hpp>
#include <kernel/lafem/arch/matvecmult_cscr_dense.hpp>
#include <kernel/lafem/arch/scale_dense.hpp>

// includes, system
#include <fstream>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief CSCR based sparse matrix.
     *
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse compressed row format.\n
     * This format is very similar to the standard CSR format but stores explicit row numbers for each non empty row.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: row start index (including matrix end index)\n
     * _indices[1]: row number of each non empty row\n
     * _indices[2]: column index per non zero element \n
     *
     * _scalar_index[0]: row count \n
     * _scalar_index[1]: column count \n
     * _scalar_index[2]: non zero row count
     * _scalar_index[3]: non zero element count (used elements) \n
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename DT_, typename IT_ = Index>
    class SparseMatrixCSCR : public Container<DT_, IT_>
    {
    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Value type, meaning the type of each 'block'
      typedef DT_ ValueType;
      /// Compatible L-vector type
      typedef DenseVector<DT_, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<DT_, IT_> VectorTypeR;
      /// Our used layout type
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_cscr;
      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = SparseMatrixCSCR<DT2_, IT2_>;

      /// this typedef lets you create a matrix container with different Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

    protected:
      Index & _rows()
      {
        return this->_scalar_index.at(0);
      }

      Index & _cols()
      {
        return this->_scalar_index.at(1);
      }

      Index & _nonzero_rows()
      {
        return this->_scalar_index.at(2);
      }

      Index & _nzes()
      {
        return this->_scalar_index.at(3);
      }

      void _alloc(Index rows, Index cols, Index nz_rows, Index nzes)
      {
        XASSERT(this->_scalar_index.empty());
        XASSERT(this->_indices.empty());
        XASSERT(this->_elements.empty());

        this->_scalar_index.push_back(rows);
        this->_scalar_index.push_back(cols);
        this->_scalar_index.push_back(nz_rows);
        this->_scalar_index.push_back(nzes);

        if(rows <= Index(0))
          return;

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * (nz_rows + 1)));
        this->_indices_size.push_back(nz_rows + 1);

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * nz_rows));
        this->_indices_size.push_back(nz_rows);

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * nzes));
        this->_indices_size.push_back(nzes);

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * nzes));
        this->_elements_size.push_back(nzes);
      }

    public:
      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      SparseMatrixCSCR()
      {
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       *
       * Creates an empty matrix.
       * Because SparseMatrixCSCR is a read-only container, it stays empty.
       *
       * \note This matrix does not allocate any memory
       */
      explicit SparseMatrixCSCR(Index rows_in, Index columns_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows The row count of the created matrix.
       * \param[in] cols The column count of the created matrix.
       * \param[in] nz_rows The amount of non zero rows of the created matrix.
       * \param[in] nzes The amount of non zero elements of the created matrix.
       *
       * Creates an empty (but allocated) matrix.
       *
       * \note The allocated memory will not be initialized.
       */
      explicit SparseMatrixCSCR(Index rows, Index cols, Index nz_rows, Index nzes)
      {
        XASSERT(rows != Index(0) && cols != Index(0));
        this->_alloc(rows, cols, nz_rows, nzes);
      }

      /**
       * \brief Constructor
       *
       * \param[in] layout_in The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixCSCR(const SparseLayout<IT_, layout_id> & layout_in)
      {
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());

        for(const auto& idx : layout_in.get_indices())
          this->_indices.push_back(idx.attach());

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * _nzes()));
        this->_elements_size.push_back(_nzes());
      }

      /**
       * \brief Constructor
       *
       * \param[in] other The source matrix.
       *
       * Creates a CSCR matrix based on the source matrix.
       */
      template <typename MT_>
      explicit SparseMatrixCSCR(const MT_ & other)
      {
        this->convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a CSCR matrix based on the source file.
       */
      explicit SparseMatrixCSCR(FileMode mode, const String& filename)
      {
        this->read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a CSCR matrix based on the source filestream.
       */
      explicit SparseMatrixCSCR(FileMode mode, std::istream& file)
      {
        this->read_from(mode, file);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] col_idx_in Vector with column indices.
       * \param[in] val_in Vector with non zero elements.
       * \param[in] row_ptr_in Vector with start indices of all rows into the val/col_idx arrays.
       * Note that this vector must also contain the end index of the last row and thus has a size of row_count + 1.
       * \param[in] row_idx_in Vector with the row number of each used row.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixCSCR(const Index rows_in, const Index columns_in,
        DenseVector<IT_, IT_> & row_ptr_in, DenseVector<IT_, IT_> & row_idx_in,
        DenseVector<IT_, IT_> & col_idx_in, DenseVector<DT_, IT_> & val_in)
      {
        /// \todo maybe create empty matrix if col_idx and val and row_ptr inputs are all three empty
        XASSERT(col_idx_in.size() > 0);
        XASSERT(val_in.size() == col_idx_in.size());
        XASSERT(row_ptr_in.size() > 0);
        XASSERT(row_idx_in.size() > 0);
        XASSERT(row_idx_in.size() + 1 == row_ptr_in.size());

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(row_idx_in.size());
        this->_scalar_index.push_back(val_in.size());

        this->_indices.push_back(row_ptr_in.elements_arbiter().attach());
        this->_indices_size.push_back(row_ptr_in.size());
        this->_indices.push_back(row_idx_in.elements_arbiter().attach());
        this->_indices_size.push_back(row_idx_in.size());
        this->_indices.push_back(col_idx_in.elements_arbiter().attach());
        this->_indices_size.push_back(col_idx_in.size());
        this->_elements.push_back(val_in.elements_arbiter().attach());
        this->_elements_size.push_back(val_in.size());
      }

      /**
       * \brief Constructor
       *
       * \param[in] csr The source csr matrix.
       * \param[in] non_zero_rows The mirror describing the rows to be selected
       *
       * Creates a matrix with selected rows from a given csr matrix.
       */
      explicit SparseMatrixCSCR(const SparseMatrixCSR<DT_, IT_> & csr, const VectorMirror<DT_, IT_> & non_zero_rows)
      {
        XASSERT(non_zero_rows.num_indices() > 0);

        const Memory::TypedView<IT_> row_ptr_in = csr.row_ptr_view_r();
        const Memory::TypedView<DT_> val_in = csr.val_view_r();
        const Memory::TypedView<IT_> col_idx_in = csr.col_idx_view_r();
        const Memory::TypedView<IT_> mir_idx = non_zero_rows.indices_view_r();

        // count number of non-zeros
        Index nz_rows = non_zero_rows.num_indices();
        Index nnzes(0);
        for (Index i(0) ; i < nz_rows ; ++i)
          nnzes += row_ptr_in[mir_idx[i]+1] - row_ptr_in[mir_idx[i]];

        this->_alloc(csr.num_rows(), csr.num_cols(), nz_rows, nnzes);

        Memory::TypedView<DT_> val = this->val_view_w();
        Memory::TypedView<IT_> row_ptr = this->row_ptr_view_w();
        Memory::TypedView<IT_> col_idx = this->col_idx_view_w();
        Memory::TypedView<IT_> row_idx = this->row_idx_view_w();

        row_ptr[0] = IT_(0);
        for(Index i(0); i < nz_rows; ++i)
        {
          IT_ row = mir_idx[i];
          IT_ j = row_ptr[i];
          row_idx[i] = row;
          for(IT_ k(row_ptr_in[row]); k < row_ptr_in[row+1]; ++k, ++j)
          {
            val[j] = val_in[k];
            col_idx[j] = col_idx_in[k];
          }
          row_ptr[i+1] = j;
        }

        XASSERT(row_ptr[nz_rows] == nnzes);
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix from
       *
       * Creates a CSCR matrix based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit SparseMatrixCSCR(const Adjacency::Graph & graph)
      {
        // get number of rows, columns and indices
        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nzes = graph.get_num_indices();

        if (num_nzes <= Index(0))
        {
          this->_alloc(num_rows, num_cols, Index(0), Index(0));
          return;
        }

        // get graph arrays
        const Index* dom_ptr = graph.get_domain_ptr();
        const Index* img_idx = graph.get_image_idx();

        // count number of non-empty rows
        Index num_nzrs = Index(0);
        for(Index i(0); i < num_rows; ++i)
        {
          num_nzrs += Index(dom_ptr[i] < dom_ptr[i+1] ? 1 : 0);
        }

        // allocate output matrix
        this->_alloc(num_rows, num_cols, num_nzrs, num_nzes);

        // get matrix arrays
        Memory::TypedView<IT_> trow_ptr = this->row_ptr_view_w();
        Memory::TypedView<IT_> trow_idx = this->row_idx_view_w();
        Memory::TypedView<IT_> tcol_idx = this->col_idx_view_w();

        // fill arrays
        trow_ptr[0] = IndexType(dom_ptr[0]);
        for(Index i(0), j(0); i < num_rows; ++i)
        {
          if(dom_ptr[i] < dom_ptr[i + 1])
          {
            ASSERT(j < num_nzrs);
            trow_idx[  j] = IndexType(i);
            trow_ptr[++j] = IndexType(dom_ptr[i+1]);
          }
        }

        for(Index k(0); k < num_nzes; ++k)
          tcol_idx[k] = IndexType(img_idx[k]);
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixCSCR(std::vector<char> input)
      {
        this->deserialize<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixCSCR(SparseMatrixCSCR && other) :
        Container<DT_, IT_>(std::forward<SparseMatrixCSCR>(other))
      {
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixCSCR & operator= (SparseMatrixCSCR && other)
      {
        this->move(std::forward<SparseMatrixCSCR>(other));

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a clone of this container.
       *
       * \param[in] clone_mode The actual cloning procedure.
       * \returns The created clone.
       *
       */
      SparseMatrixCSCR clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        SparseMatrixCSCR t;
        t.clone(*this, clone_mode);
        return t;
      }

      /** \brief Clone operation
       *
       * Create a clone of another container.
       *
       * \param[in] other The source container to create the clone from.
       * \param[in] clone_mode The actual cloning procedure.
       *
       */
      template<typename DT2_, typename IT2_>
      void clone(const SparseMatrixCSCR<DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Weak)
      {
        Container<DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename DT2_, typename IT2_>
      void convert(const SparseMatrixCSCR<DT2_, IT2_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] a The input matrix.
       *
       * Converts any matrix to SparseMatrixCSCR-format
       */
      template <typename MT_>
      void convert(const MT_ & a)
      {
        XASSERT((void*)this != (void*)&a);

        this->clear();

        const Index num_rows_a = a.num_rows_raw();
        const Index num_cols_a = a.num_cols_raw();
        const Index num_nzes_a = a.num_nzes_raw();
        if(num_nzes_a <= Index(0))
        {
          this->_alloc(num_rows_a, num_cols_a, Index(0), Index(0));
          return;
        }

        Index num_nz_rows_a(0);
        for (Index i(0); i < num_rows_a; ++i)
        {
          num_nz_rows_a += Index(a.row_degree(i) > Index(0) ? 1 : 0);
        }

        this->_alloc(num_rows_a, num_cols_a, num_nz_rows_a, num_nzes_a);

        Memory::TypedView<IT_> row_ptr = this->row_ptr_view_w();
        Memory::TypedView<IT_> row_idx = this->row_idx_view_w();
        Memory::TypedView<IT_> col_idx = this->col_idx_view_w();
        Memory::TypedView<DT_> val = this->val_view_w();

        row_ptr[0] = Index(0);
        for(Index row(0), k(0); row < num_rows_a; ++row)
        {
          Index nz = a.row_degree(row);
          if(nz > Index(0))
          {
            a.get_row_col_indices(row, &col_idx[row_ptr[k]], IT_(0));
            a.get_row_values(row, &val[row_ptr[k]]);
            row_idx[k] = IT_(row);
            row_ptr[k+1] = row_ptr[k] + IT_(nz);
            ++k;
          }
        }
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixCSCR & operator= (const SparseLayout<IT_, layout_id> & layout_in)
      {
        this->clear();

        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());

        for(const auto& idx : layout_in.get_indices())
          this->_indices.push_back(idx.attach());

        this->_elements.push_back(Memory::Arbiter(_nzes() * sizeof(DT_)));
        this->_elements_size.push_back(_nzes());

        return *this;
      }

      /**
       * \brief Deserialization of complete container entity.
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialize(std::vector<char> input)
      {
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_cscr, input);
      }

      /**
       * \brief Serialization of complete container entity.
       *
       * \param[in] config LAFEM::SerialConfig, a struct describing the serialize configuration.
       * \note the corresponding configure flags 'zlib' and/or 'zfp' need to be added in the build-id at the configure call.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialize for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialize(const LAFEM::SerialConfig& config = SerialConfig()) const
      {
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_cscr, config);
      }

      /**
       * \brief Read in matrix from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, const String& filename)
      {
        std::ios_base::openmode bin = std::ifstream::in | std::ifstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ifstream::in;
        std::ifstream file(filename.c_str(), bin);
        if (! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);
        read_from(mode, file);
        file.close();
      }

      /**
       * \brief Read in matrix from stream.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be read in.
       */
      void read_from(FileMode mode, std::istream& file)
      {
        switch(mode)
        {
        case FileMode::fm_cscr:
          [[fallthrough]];

        case FileMode::fm_binary:
          this->template _deserialize<double, std::uint64_t>(FileMode::fm_cscr, file);
          break;

        default:
          XABORTM("Filemode not supported!");
        }
      }

      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, const String& filename) const
      {
        std::ios_base::openmode bin = std::ofstream::out | std::ofstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ofstream::out;
        std::ofstream file;
        char* buff = nullptr;
        if(mode == FileMode::fm_mtx)
        {
          buff = new char[LAFEM::file_out_stream_buffer_size];
          file.rdbuf()->pubsetbuf(buff, LAFEM::file_out_stream_buffer_size);
        }
        file.open(filename.c_str(), bin);
        if(! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);
        write_out(mode, file);
        file.close();
        delete[] buff;
      }

      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be written to.
       */
      void write_out(FileMode mode, std::ostream& file) const
      {
        switch(mode)
        {
        case FileMode::fm_cscr:
          [[fallthrough]];

        case FileMode::fm_binary:
          this->template _serialize<double, std::uint64_t>(FileMode::fm_cscr, file);
          break;

        default:
          XABORTM("Filemode not supported!");
        }
      }

      /**
       * \brief Retrieve convenient sparse matrix layout object.
       *
       * \return An object containing the sparse matrix layout.
       */
      SparseLayout<IT_, layout_id> layout() const
      {
        return SparseLayout<IT_, layout_id>(this->_indices, this->_indices_size, this->_scalar_index);
      }

      /// Checks whether the matrix is empty, i.e. whether it is a 0x0 matrix
      bool empty() const
      {
        return this->_scalar_index.empty() || (this->_scalar_index.at(0) == Index(0));
      }

      /// Checks whether the matrix has no non-zero entries, i.e. whether it is a null matrix
      bool hollow() const
      {
        return this->_scalar_index.empty() || (this->_scalar_index.at(3) == Index(0));
      }

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count.
       */
      Index num_rows() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(0);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      Index num_cols() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve used matrix non zero row count.
       *
       * \returns Used matrix non zero row count.
       */
      Index num_nonzero_rows() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      Index num_nzes() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(3);
      }

      /// Returns the scalar matrix row count
      Index num_rows_raw() const
      {
        return this->num_rows();
      }

      /// Returns the scalar matrix column count
      Index num_cols_raw() const
      {
        return this->num_cols();
      }

      /// Returns the scalar matrix nonzero element count
      Index num_nzes_raw() const
      {
        return this->num_nzes();
      }

      /// Returns a reference to the row-pointer array arbiter
      Memory::Arbiter& row_ptr_arbiter()
      {
        return this->_indices.at(0);
      }

      /// Returns a reference to the row-pointer array arbiter
      const Memory::Arbiter& row_ptr_arbiter() const
      {
        return this->_indices.at(0);
      }

      /// Returns a reference to the row-index array arbiter
      Memory::Arbiter& row_idx_arbiter()
      {
        return this->_indices.at(1);
      }

      /// Returns a reference to the row-index array arbiter
      const Memory::Arbiter& row_idx_arbiter() const
      {
        return this->_indices.at(1);
      }

      /// Returns a reference to the column-index array arbiter
      Memory::Arbiter& col_idx_arbiter()
      {
        return this->_indices.at(2);
      }

      /// Returns a reference to the column-index array arbiter
      const Memory::Arbiter& col_idx_arbiter() const
      {
        return this->_indices.at(2);
      }

      /// Returns a reference to the values array arbiter
      Memory::Arbiter& val_arbiter()
      {
        return this->_elements.front();
      }

      /// Returns a reference to the values array arbiter
      const Memory::Arbiter& val_arbiter() const
      {
        return this->_elements.front();
      }

      Memory::TypedView<DT_> val_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<DT_> val_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<DT_> val_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<DT_> val_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, acc));
      }

      Memory::TypedView<IT_> row_ptr_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<IT_> row_ptr_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<IT_> row_ptr_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<IT_> row_ptr_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, acc));
      }

      Memory::TypedView<IT_> row_idx_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::read));
      }

      Memory::TypedView<IT_> row_idx_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::write));
      }

      Memory::TypedView<IT_> row_idx_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<IT_> row_idx_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, acc));
      }

      Memory::TypedView<IT_> col_idx_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(2).view(loc, Memory::Access::read));
      }

      Memory::TypedView<IT_> col_idx_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(2).view(loc, Memory::Access::write));
      }

      Memory::TypedView<IT_> col_idx_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(2).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<IT_> col_idx_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(2).view(loc, acc));
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SparseMatrixCSCR";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const SparseMatrixCSCR & x, bool full = false)
      {
        this->_copy_content(x, full);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow y + \alpha~ x\f$
       *
       * \param[in] x The first summand matrix to be scaled.
       * \param[in] y The second summand matrix
       * \param[in] alpha A scalar to multiply x with.
       *
       * \warning All three matrices must have the same non zero layout. This operation assumes this silently and does not check this on its own!
       */
      void axpy(const SparseMatrixCSCR & x, const DT_ alpha = DT_(1))
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Matrix rows do not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Matrix columns do not match!");
        XASSERTM(x.num_nzes() == this->num_nzes(), "Matrix nonzeros do not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::AxpyDense::template exec<DT_>(this->val_arbiter(), alpha, x.val_arbiter(), this->num_nzes());

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const SparseMatrixCSCR & x, const DT_ alpha)
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Row count does not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Column count does not match!");
        XASSERTM(x.num_nzes() == this->num_nzes(), "Nonzero count does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::ScaleDense::template exec<DT_>(this->val_arbiter(), x.val_arbiter(), alpha, this->num_nzes());

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes());
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::Norm2Dense::template exec<DT_>(this->val_arbiter(), this->num_nzes());

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the maximum relative difference of this matrix and another one
       * y.max_rel_diff(x) returns  \f$ \max_{0\leq i < n}\frac{|x_i-y_i|}{\max{|x_i|+|y_i|, eps}} \f$
       *
       * \return The largest relative difference.
       */
      DT_ max_rel_diff(const SparseMatrixCSCR& x) const
      {
        XASSERTM(x.num_nzes() == this->num_nzes(), "Nonzero count does not match!");
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MaxRelDiffDense::template exec<DT_>(this->val_arbiter(), x.val_arbiter(), this->num_nzes());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Checks if the structural layout of this matrix matches that of another matrix.
       * This excludes comparison of the actual data values.
       *
       * \param[in] x The matrix to compare this matrix to
       *
       * \returns true if the layouts match, false otherwise.
       */
      bool same_layout(const SparseMatrixCSCR& x) const
      {
        if (this->num_rows() != x.num_rows())
          return false;
        if (this->num_cols() != x.num_cols())
          return false;
        if (this->num_nzes() != x.num_nzes())
          return false;
        if (this->num_nonzero_rows() != x.num_nonzero_rows())
          return false;

        if(this->num_nzes() == Index(0))
          return true;

        // check if the arbiters for row_ptr, row_idx and col_idx are the same
        if((this->row_ptr_arbiter() == x.row_ptr_arbiter()) &&
          (this->row_idx_arbiter() == x.row_idx_arbiter()) &&
          (this->col_idx_arbiter() == x.col_idx_arbiter()))
          return true;

        const Memory::TypedView<IT_> row_ptr_a = this->row_ptr_view_r();
        const Memory::TypedView<IT_> row_idx_a = this->row_idx_view_r();
        const Memory::TypedView<IT_> col_idx_a = this->col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_b = x.row_ptr_view_r();
        const Memory::TypedView<IT_> row_idx_b = x.row_idx_view_r();
        const Memory::TypedView<IT_> col_idx_b = x.col_idx_view_r();

        const Index nzr = this->num_nonzero_rows();
        for (Index i(0) ; i < nzr; ++i)
        {
          if(row_ptr_a[i] != row_ptr_b[i])
            return false;
          if(row_idx_a[i] != row_idx_b[i])
            return false;
        }
        if(row_ptr_a[nzr] != row_ptr_b[nzr])
          return false;

        const Index nze = this->num_nzes();
        for (Index i(0) ; i < nze ; ++i)
        {
          if (col_idx_a[i] != col_idx_b[i])
            return false;
        }

        return true;
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<DT_, IT_> & r, const DenseVector<DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");

        if(this->empty())
          return;
        else if(this->hollow())
        {
          r.format();
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultCSCRDense::template exec<DT_, IT_>(r.elements_arbiter(), DT_(1), x.elements_arbiter(), Memory::Arbiter(),
          this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(), this->row_idx_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nonzero_rows(), this->num_nzes(), false);

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
      * \brief Calculate \f$ r \leftarrow this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      */
      void apply_transposed(DenseVector<DT_, IT_> & r, const DenseVector<DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");

        if(this->empty())
          return;
        else if(this->hollow())
        {
          r.format();
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultCSCRDense::template exec<DT_, IT_>(r.elements_arbiter(), DT_(1), x.elements_arbiter(), Memory::Arbiter(),
          this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(), this->row_idx_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nonzero_rows(), this->num_nzes(), true);

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
        DenseVector<DT_, IT_> & r,
        const DenseVector<DT_, IT_> & x,
        const DenseVector<DT_, IT_> & y,
        const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows(), "Vector size of y does not match!");

        if(this->empty())
          return;
        // we check alpha strictly for equality to 0, because testing if |alpha| < eps may
        // lead to false triggers if the matrix/vector contents are also < eps
        else if(this->hollow() || (alpha == DT_(0)))
        {
          r.copy(y);
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultCSCRDense::template exec<DT_, IT_>(r.elements_arbiter(), alpha, x.elements_arbiter(), y.elements_arbiter(),
          this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(), this->row_idx_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nonzero_rows(), this->num_nzes(), false);

        TimeStamp ts_stop;

        Statistics::add_flops( (this->num_nzes() + this->num_rows()) * 2 );
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha~ this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      * \param[in] y The summand vector.
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(
        DenseVector<DT_, IT_> & r,
        const DenseVector<DT_, IT_> & x,
        const DenseVector<DT_, IT_> & y,
        const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_cols(), "Vector size of y does not match!");

        if(this->empty())
          return;
        // we check alpha strictly for equality to 0, because testing if |alpha| < eps may
        // lead to false triggers if the matrix/vector contents are also < eps
        else if(this->hollow() || (alpha == DT_(0)))
        {
          r.copy(y);
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultCSCRDense::template exec<DT_, IT_>(r.elements_arbiter(), alpha, x.elements_arbiter(), y.elements_arbiter(),
          this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(), this->row_idx_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nonzero_rows(), this->num_nzes(), true);

        TimeStamp ts_stop;

        Statistics::add_flops( (this->num_nzes() + this->num_rows()) * 2 );
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }
      ///@}

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(this->num_rows());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->num_cols());
      }

      /**
       * \brief Prints this sparse matrix in human-readable format into an output stream
       *
       * \param[inout] os
       * The output stream to print to
       *
       * \param[in]
       * If \c true, then the matrix will be printed as a dense matrix including all implicit zero entries,
       * otherwise only the non-zero entries of each non-zero row prefixed by the corresponding column index are printed.
       */
      void print(std::ostream & os, bool print_dense) const
      {
        if(this->empty())
        {
          os << "[]";
          return;
        }

        const Memory::TypedView<IT_> row_ptr = this->row_ptr_view_r();
        const Memory::TypedView<IT_> row_idx = this->row_idx_view_r();
        const Memory::TypedView<IT_> col_idx = this->col_idx_view_r();
        const Memory::TypedView<DT_> val = this->val_view_r();
        const Index nrows = this->num_rows();
        const Index ncols = this->num_cols();
        const Index nzrows = this->num_nonzero_rows();

        IT_ l = 0;
        os << "[\n";
        for (Index i(0) ; i < nzrows ; ++i)
        {
          if(print_dense)
          {
            // write leading zero rows
            for(; l < row_idx[i]; ++l)
            {
              os << "[";
              for(Index j = 0; j < ncols; ++j)
                os << "  " << DT_(0);
              os << "]\n";
            }
          }
          else
          {
            // write row index
            os << row_idx[i] << ':';
          }

          os << "[";
          IT_ k = 0;
          for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
          {
            // write leading zeros
            if(print_dense)
            {
              for(; k < col_idx[j]; ++k)
                os << "  " << DT_(0);
            }
            else
            {
              os << "  " << col_idx[j] << ':';
            }

            // write value
            os << "  " << val[j];
            ++k;
          }

          // write trailing zeros
          if(print_dense)
          {
            for(; k < ncols; ++k)
              os << "  " << DT_(0);
            ++l;
          }
          os << "]\n";
        }

        // write trailing zero rows
        if(print_dense)
        {
          for(; l < nrows; ++l)
          {
            os << "[";
            for(Index j = 0; j < ncols; ++j)
              os << "  " << DT_(0);
            os << "]\n";
          }
        }

        os << "]";
      }

      /**
      * \brief SparseMatrixCSR streaming operator
      *
      * \param[in] os The target stream.
      * \param[in] b The matrix to be streamed.
      */
      friend std::ostream & operator<< (std::ostream & os, const SparseMatrixCSCR & b)
      {
        b.print(os, true);
        return os;
      }
    }; //SparseMatrixCSCR

#ifdef FEAT_EICKT
    extern template class SparseMatrixCSCR<float, std::uint32_t>;
    extern template class SparseMatrixCSCR<double, std::uint32_t>;
    extern template class SparseMatrixCSCR<float, std::uint64_t>;
    extern template class SparseMatrixCSCR<double, std::uint64_t>;
#endif

  } // namespace LAFEM
} // namespace FEAT
