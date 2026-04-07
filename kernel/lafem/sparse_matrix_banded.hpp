// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/axpy_dense.hpp>
#include <kernel/lafem/arch/matvecmult_banded_dense.hpp>
#include <kernel/lafem/arch/norm2_dense.hpp>
#include <kernel/lafem/arch/scale_dense.hpp>

// includes, system
#include <fstream>
#include <set>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Sparse banded matrix
     *
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its diagonal entries
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: vector of offsets (bottom-left diagonal has offset 0,
     *                                 main diagonal has offset rows - 1 and
     *                                 top-right diagonal has offset row + columns - 2)\n
     *
     * _scalar_index[0]: row count \n
     * _scalar_index[1]: column count \n
     * _scalar_index[2]: non zero element count (used elements) \n
     * _scalar_index[3]: number of offsets \n
     *
     * This class saves a sparse-matrix with a banded structure. For each diagonal of
     * the matrix with non-zero elements there must be reserved memory for the whole
     * diagonal. For faster access on the matrix-elements each diagonal get the virtual
     * length of the row-count of the matrix. They are enlarged to the left and right
     * side of the matrix as shown in the following layout.
     \verbatim
           +--                  --+
      \    | \           \      \ |
       \   |\ \           \      \|
        \  | \ \           \      |
         \ |  \ \           \     |\
          \|   \ \           \    | \
           |    \ \           \   |  \
           |\    \ \           \  |   \
           +--                  --+
     \endverbatim
     * To get the position of the diagonals in the matrix, the matching offsets are
     * saved from left to right in the offsets-array.
     * - The first diagonal is the one at the bottom-left and gets the offset = 0,
     * - the main diagonal has the offset = rows - 1
     * - and the last offset at the top-right has the offset = rows + columns - 2.
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Christoph Lohmann
     */
    template <typename DT_, typename IT_ = Index>
    class SparseMatrixBanded : public Container<DT_, IT_>
    {
    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Compatible L-vector type
      typedef DenseVector<DataType, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<DataType, IT_> VectorTypeR;
      /// Our used layout type
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_banded;
      /// our value type
      typedef DT_ ValueType;
      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = SparseMatrixBanded<DT2_, IT2_>;

      /// this typedef lets you create a matrix container with different Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

    protected:
      Index & _rows()
      {
        return this->_scalar_index.at(0);
      }

      Index & _cols()
      {
        return this->_scalar_index.at(1);
      }

      Index & _bands()
      {
        return this->_scalar_index.at(2);
      }

      Index & _nzes()
      {
        return this->_scalar_index.at(3);
      }


      /// Returns first row-index of the diagonal matching to the offset i
      template <typename ITX_>
      inline Index _start_offset(const Index i, const ITX_ * const offsets,
        const Index rows_in, const Index columns_in, const Index noo) const
      {
        if (i == Index(-1))
        {
          return rows_in;
        }
        else if (i == noo)
        {
          return Index(0);
        }
        else
        {
          return Math::max(columns_in + Index(1), rows_in + columns_in - Index(offsets[i])) - columns_in - Index(1);
        }
      }

      /// Returns last row-index of the diagonal matching to the offset i
      template <typename ITX_>
      inline Index _end_offset(const Index i, const ITX_ * const offsets,
        const Index rows_in, const Index columns_in, const Index noo) const
      {
        if (i == Index (-1))
        {
          return rows_in - 1;
        }
        else if (i == noo)
        {
          return Index(-1);
        }
        else
        {
          return Math::min(rows_in, columns_in + rows_in - Index(offsets[i]) - Index(1)) - Index(1);
        }
      }

    public:
      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      SparseMatrixBanded()
      {
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] bands_in The band count of the created matrix
       * \param[in] nzes_in The non-zero entry count of the created matrix
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixBanded(const Index rows_in, const Index columns_in, const Index bands_in, const Index nzes_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(bands_in);
        this->_scalar_index.push_back(nzes_in);

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * rows_in * bands_in,
          Memory::Location::none, Memory::Init::format_to_zero));
        this->_elements_size.push_back(rows_in * bands_in);
        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * rows_in * bands_in));
        this->_indices_size.push_back(rows_in * bands_in);
      }

      /**
       * \brief Constructor
       *
       * \param[in] layout_in The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixBanded(const SparseLayout<IT_, layout_id> & layout_in)
      {
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());

        for(const auto& idx : layout_in.get_indices())
          this->_indices.push_back(idx.attach());

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * _rows() * _bands(),
          Memory::Location::none, Memory::Init::format_to_zero));
        this->_elements_size.push_back(_rows() * _bands());
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] val_in The vector with non zero elements.
       * \param[in] bands_in The vector of band offsets.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixBanded(const Index rows_in, const Index columns_in,
                                  DenseVector<DT_, IT_> & val_in,
                                  DenseVector<IT_, IT_> & bands_in)
      {
        if (val_in.size() != rows_in * bands_in.size())
        {
          XABORTM("Size of values does not match to number of offsets and row count!");
        }

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);

        Index tnonzeros(0);

        Memory::TypedView<IT_> vbands(bands_in.elements_view_r());
        for (Index i(0); i < bands_in.size(); ++i)
        {
          const Index toffset(Index(vbands(i)));

          if (toffset + Index(2) > rows_in + columns_in)
          {
            XABORTM("Offset out of matrix!");
          }

          tnonzeros += columns_in + Math::min(rows_in, columns_in + rows_in - toffset - 1) - Math::max(columns_in + rows_in - toffset - 1, columns_in);
        }
        vbands.release();

        this->_scalar_index.push_back(bands_in.size());
        this->_scalar_index.push_back(tnonzeros);

        this->_elements.push_back(val_in.elements_arbiter().attach());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(bands_in.elements_arbiter().attach());
        this->_indices_size.push_back(bands_in.size());
      }

      /**
       * \brief Constructor
       *
       * \param[in] other The source matrix.
       *
       * Creates a Banded matrix based on the source matrix.
       */
      template <typename MT_>
      explicit SparseMatrixBanded(const MT_ & other)
      {
        this->convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix from
       *
       * Creates a matrix based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit SparseMatrixBanded(const Adjacency::Graph & graph)
      {
        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();

        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());

        std::set<IT_> moffsets;

        for(Index row(0); row < num_rows; ++row)
        {
          for(Index j(dom_ptr[row]); j < dom_ptr[row+1]; ++j)
          {
            moffsets.insert(IT_(num_rows - 1 + img_idx[j] - row));
          }
        }

        DenseVector<IT_, IT_> toffsets(Index(moffsets.size()));
        Memory::TypedView<IT_> ptoffsets = toffsets.elements_view_w();

        IT_ idx(0);
        for (auto off : moffsets)
        {
          ptoffsets[idx] = off;
          ++idx;
        }

        DenseVector<DT_, IT_> tval(Index(moffsets.size()) * num_rows, DT_(0));

        this->move(SparseMatrixBanded(num_rows, num_cols, tval, toffsets));
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a banded matrix based on the source filestream.
       */
      explicit SparseMatrixBanded(FileMode mode, String filename)
      {
        this->read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a banded matrix based on the source filestream.
       */
      explicit SparseMatrixBanded(FileMode mode, std::istream& file)
      {
        this->read_from(mode, file);
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixBanded(std::vector<char> input)
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
      SparseMatrixBanded(SparseMatrixBanded && other) :
        Container<DT_, IT_>(std::forward<SparseMatrixBanded>(other))
      {
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixBanded & operator= (SparseMatrixBanded && other)
      {
        this->move(std::forward<SparseMatrixBanded>(other));
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
      SparseMatrixBanded clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        SparseMatrixBanded t;
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
      void clone(const SparseMatrixBanded<DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Weak)
      {
        Container<DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixBanded & operator= (const SparseLayout<IT_, layout_id> & layout_in)
      {
        this->_elements.clear();
        this->_indices.clear();
        this->_elements_size.clear();
        this->_indices_size.clear();
        this->_scalar_index.clear();

        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());

        for(const auto& idx : layout_in.get_indices())
          this->_indices.push_back(idx.attach());

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * _rows() * _bands(),
          Memory::Location::none, Memory::Init::format_to_zero));
        this->_elements_size.push_back(_rows() * _bands());

        return *this;
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename DT2_, typename IT2_>
      void convert(const SparseMatrixBanded<DT2_, IT2_> & other)
      {
        this->assign(other);
      }

      template <typename DT2_, typename IT2_>
      void convert(const SparseMatrixCSR<DT2_, IT2_> & csr)
      {
        XASSERT(csr.num_nzes() > 0);

        const Memory::TypedView<IT_> row_ptr = csr.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx = csr.col_idx_view_r();
        const Memory::TypedView<DT_> csr_val = csr.val_view_r();

        std::set<IT_> offset_set;
        Index nrows(csr.num_rows());
        Index ncolumns(csr.num_cols());
        for (Index row(0) ; row < nrows ; ++row)
        {
          for (Index i(row_ptr[row]) ; i < row_ptr[row+1] ; ++i)
          {
            //compute band offset for each non zero entry
            IT_ off = IT_(col_idx[i] - row + nrows - 1);
            offset_set.insert(off);
          }
        }

        DenseVector<DT_, IT_> val_new(Index(offset_set.size()) * nrows, DT_(0));
        Memory::TypedView<DT_ > pval = val_new.elements_view_w();
        for (Index row(0) ; row < nrows ; ++row)
        {
          for (Index i(row_ptr[row]) ; i < row_ptr[row+1] ; ++i)
          {
            Index col = col_idx[i];
            Index offset = (col - row + nrows - 1);
            Index band = Index(std::distance(offset_set.begin(), offset_set.find(IT_(offset))));
            pval[band * nrows + row] = csr_val[i];
          }
        }

        DenseVector<IT_, IT_> offsets_new(Index(offset_set.size()));
        Memory::TypedView<IT_> poffsets = offsets_new.elements_view_w();
        auto it_offset = offset_set.begin();
        for (Index i(0) ; i < offset_set.size() ; ++i, ++it_offset)
        {
          poffsets[i] = *it_offset;
        }
        offset_set.clear();

        pval.release();
        poffsets.release();

        SparseMatrixBanded<DT_, IT_> temp(nrows, ncolumns, val_new, offsets_new);
        this->move(std::move(temp));
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
        case FileMode::fm_bm:
          [[fallthrough]];

        case FileMode::fm_binary:
          if (! std::is_same<DT_, double>::value)
            std::cout<<"Warning: You are writing out a banded matrix that is not double precision!\n";
          this->template _serialize<double, std::uint64_t>(FileMode::fm_bm, file);
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
       * \brief Retrieve number of bands.
       *
       * \returns Number of bands.
       */
      Index num_bands() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * This only contains the number of non-zero elements within the matrix region.
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

      /**
       * \brief Retrieve total element count.
       *
       * This contains the number of non-zero elements within and outside of the matrix region.
       *
       * \returns Total element count.
       */
      Index val_size() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(0) * this->_scalar_index.at(2);
      }

      /// Checks whether the matrix is empty, i.e. whether it is a 0x0 matrix
      bool empty() const
      {
        return this->_scalar_index.empty() || (this->_scalar_index.at(0) == Index(0));
      }

      /// Checks whether the matrix has no non-zero entries, i.e. whether it is a null matrix
      bool hollow() const
      {
        return this->_scalar_index.empty() || (this->_scalar_index.at(2) == Index(0));
      }

      /// Returns a reference to the band offsets array arbiter
      Memory::Arbiter& offsets_arbiter()
      {
        return this->_indices.at(0);
      }

      /// Returns a reference to the band offsets array arbiter
      const Memory::Arbiter& offsets_arbiter() const
      {
        return this->_indices.at(0);
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

      Memory::TypedView<IT_> offsets_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<IT_> offsets_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<IT_> offsets_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<IT_> offsets_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, acc));
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SparseMatrixBanded";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const SparseMatrixBanded & x, bool full = false)
      {
        this->_copy_content(x, full);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow this + \alpha~ x\f$
       *
       * \param[in] x The first summand matrix to be scaled.
       * \param[in] y The second summand matrix
       * \param[in] alpha A scalar to multiply x with.
       */
      void axpy(const SparseMatrixBanded & x, const DT_ alpha = DT_(1))
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Matrix rows do not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Matrix columns do not match!");
        XASSERTM(x.num_bands() == this->num_bands(), "Matrix bands do not match!");
        XASSERTM(x.num_nzes() == this->num_nzes(), "Matrix nonzeros do not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::AxpyDense::template exec<DT_>(this->val_arbiter(), alpha, x.val_arbiter(), this->val_size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->val_size() * 2);
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const SparseMatrixBanded & x, const DT_ alpha)
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Matrix rows do not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Matrix columns do not match!");
        XASSERTM(x.num_bands() == this->num_bands(), "Matrix bands do not match!");
        XASSERTM(x.num_nzes() == this->num_nzes(), "Matrix nonzeros do not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::ScaleDense::template exec<DT_>(this->val_arbiter(), x.val_arbiter(), alpha, this->val_size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->val_size());
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

        DT_ result = Arch::Norm2Dense::template exec<DT_>(this->val_arbiter(), this->val_size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->val_size() * 2);
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<DT_, IT_>& r, const DenseVector<DT_, IT_>& x) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");

        if(this->empty())
          return;

        if (this->num_nzes() == Index(0))
        {
          r.format();
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultBandedDense::template exec<DT_, IT_>(
          r.elements_arbiter(), DT_(1), x.elements_arbiter(), Memory::Arbiter(),
          this->val_arbiter(), this->offsets_arbiter(), this->num_bands(), this->num_rows(), this->num_cols());

        TimeStamp ts_stop;

        Statistics::add_flops(2 * this->val_size());
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
      * \brief Calculate \f$ r \leftarrow this^\top \cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      */
      void apply_transposed(DenseVector<DT_, IT_>& r, const DenseVector<DT_, IT_>& x) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");

        if(this->empty())
          return;

        if (this->num_nzes() == Index(0))
        {
          r.format();
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultBandedDense::template exec<DT_, IT_>(
          r.elements_arbiter(), DT_(1), x.elements_arbiter(), Memory::Arbiter(),
          this->val_arbiter(), this->offsets_arbiter(), this->num_bands(), this->num_rows(), this->num_cols(), true);

        TimeStamp ts_stop;
        Statistics::add_flops( 2 * this->val_size() );
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
      void apply(DenseVector<DT_, IT_>& r,
                 const DenseVector<DT_, IT_>& x,
                 const DenseVector<DT_, IT_>& y,
                 const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows(), "Vector size of y does not match!");

        if(this->empty())
          return;

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        // we check alpha strictly for equality to 0, because testing if |alpha| < eps may
        // lead to false triggers if the matrix/vector contents are also < eps
        if((this->num_nzes() == Index(0)) || (alpha == DT_(0)))
        {
          r.copy(y);
          return;
        }

        TimeStamp ts_start;

        Arch::MatVecMultBandedDense::template exec<DT_, IT_>(
          r.elements_arbiter(), alpha, x.elements_arbiter(), y.elements_arbiter(),
          this->val_arbiter(), this->offsets_arbiter(), this->num_bands(), this->num_rows(), this->num_cols());

        TimeStamp ts_stop;

        Statistics::add_flops( 2 * (this->val_size() + this->num_rows()) );
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha~ this^\top\cdot x \f$
      *
      * \param[out] r The vector that receives the result.
      * \param[in] x The vector to be multiplied by this matrix.
      * \param[in] y The summand vector.
      * \param[in] alpha A scalar to scale the product with.
      */
      void apply_transposed(DenseVector<DT_, IT_>& r,
        const DenseVector<DT_, IT_>& x,
        const DenseVector<DT_, IT_>& y,
        const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_cols(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_rows(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_cols(), "Vector size of y does not match!");

        if(this->empty())
          return;

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        // we check alpha strictly for equality to 0, because testing if |alpha| < eps may
        // lead to false triggers if the matrix/vector contents are also < eps
        if((this->num_nzes() == Index(0)) || (alpha == DT_(0)))
        {
          r.copy(y);
          return;
        }

        TimeStamp ts_start;
        Statistics::add_flops( 2 * (this->val_size() + this->num_rows()) );

        Arch::MatVecMultBandedDense::template exec<DT_, IT_>(
          r.elements_arbiter(), alpha, x.elements_arbiter(), y.elements_arbiter(),
          this->val_arbiter(), this->offsets_arbiter(), this->num_bands(), this->num_rows(), this->num_cols(), true);

        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }
      ///@}

      /// \copydoc extract_diag()
      void extract_diag(VectorTypeL & diag) const
      {
        XASSERTM(diag.size() == num_rows(), "diag size does not match matrix row count!");
        XASSERTM(num_rows() == num_cols(), "matrix is not square!");

        if(this->empty())
          return;

        const Memory::TypedView<IT_> offs = this->offsets_view_r();
        const Memory::TypedView<DT_> vals = this->val_view_r();

        const Index nrows = this->num_rows();
        for (Index i(0) ; i < num_bands() ; ++i)
        {
          if (offs[i]+1 == nrows)
          {
            diag.elements_arbiter().copy(&vals[i*nrows], Memory::Location::main);
            break;
          }
        }
      }

      /// extract main diagonal vector from matrix
      VectorTypeL extract_diag() const
      {
        VectorTypeL diag = create_vector_l();
        extract_diag(diag);
        return diag;
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
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_bm, input);
      }

      /**
       * \brief Retrieve the maximum relative difference of this matrix and another one
       * y.max_rel_diff(x) returns  \f$ \max_{0\leq i < n}\frac{|x_i-y_i|}{\max{|x_i|+|y_i|, eps}} \f$
       *
       * \return The largest relative difference.
       */
      DT_ max_rel_diff(const SparseMatrixBanded& x) const
      {
        XASSERTM(x.num_nzes() == this->num_nzes(), "Nonzero count does not match!");
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MaxRelDiffDense::template exec<DT_>(this->val_arbiter(), x.val_arbiter(), this->val_size());

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
      bool same_layout(const SparseMatrixBanded& x) const
      {
        if (this->num_rows() != x.num_rows())
          return false;
        if (this->num_cols() != x.num_cols())
          return false;
        if (this->num_bands() != x.num_bands())
          return false;
        if (this->num_nzes() != x.num_nzes())
          return false;

        if(this->num_nzes() == Index(0))
          return true;

        // check if the arbiters for offsets are the same
        if(this->offsets_arbiter() == x.offsets_arbiter())
          return true;

        // compare the offsets
        const Memory::TypedView<IT_> offsets_a = this->offsets_view_r();
        const Memory::TypedView<IT_> offsets_b = x.offsets_view_r();
        for (Index i(0); i < this->num_bands(); ++i)
        {
          if (offsets_a[i] != offsets_b[i])
            return false;
        }

        return true;
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
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_bm, config);
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
        case FileMode::fm_bm:
          [[fallthrough]];

        case FileMode::fm_binary:
          this->template _deserialize<double, std::uint64_t>(FileMode::fm_bm, file);
          break;

        default:
          XABORTM("Filemode not supported!");
        }
      }

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
       * otherwise only the non-zero entries of each row prefixed by the corresponding column index are printed.
       */
      void print(std::ostream& os, bool print_dense) const
      {
        if(this->num_rows() <= Index(0))
        {
          os << "[]";
          return;
        }

        const Memory::TypedView<IT_> off = this->offsets_view_r();
        const Memory::TypedView<DT_> val = this->val_view_r();
        const Index nrows = this->num_rows();
        const Index ncols = this->num_cols();
        const Index nbands = this->num_bands();

        os << "[\n";
        for (Index i(0) ; i < nrows ; ++i)
        {
          os << "[";

          Index k = 0;
          for(IT_ j(0); j < nbands; ++j)
          {
            if((off[j] + i + 2 > nrows) && (off[j] + i + 1 < 2 * nrows))
            {
              const Index col_idx = Index(off[j] + i + 1) - nrows;
              if(print_dense)
              {
                for(; k < col_idx; ++k)
                  os << "  " << DT_(0);
              }
              else
              {
                os << "  " << col_idx << ":";
              }
              os << "  " << val(j * nrows + i);
              ++k;
            }
          }

          // write trailing zeros
          if(print_dense)
          {
            for(; k < ncols; ++k)
              os << "  " << DT_(0);
          }

          os << "]\n";
        }
        os << "]";
      }

      /**
       * \brief SparseMatrixBanded streaming operator
       *
       * \param[in] os The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & os, const SparseMatrixBanded & b)
      {
        b.print(os, true);
        return os;
      }

    public:
      class Adjactor
      {
      public:
        /// ImageIterator class for Adjactor interface implementation
        class ImageIterator
        {
        private:
          IT_ _k;
          const IT_ _s;
          const IT_ * const _offsets;

        public:
          ImageIterator() : _k(IT_(0)), _s(IT_(0)), _offsets(nullptr) {}

          ImageIterator(IT_ k, IT_ s, const IT_ * offsets) : _k(k), _s(s), _offsets(offsets) {}

          ImageIterator& operator=(const ImageIterator& other)
          {
            _k       = other._k;
            _s       = other._s;
            _offsets = other._offsets;
            return *this;
          }

          bool operator!=(const ImageIterator& other) const
          {
            return this->_k != other._k;
          }

          ImageIterator& operator++()
          {
            ++_k;
            return *this;
          }

          Index operator*() const
          {
            return Index(_s + _offsets[_k]);
          }
        };

      private:
        const SparseMatrixBanded& _matrix;
        const Memory::TypedView<IT_> _offsets;

      public:
        explicit Adjactor(const SparseMatrixBanded& matrix) :
          _matrix(matrix),
          _offsets(matrix.offsets_view_r())
        {
        }

        /** \copydoc Adjactor::get_num_nodes_domain() */
        inline Index get_num_nodes_domain() const
        {
          return _matrix.num_rows();
        }

        /** \copydoc Adjactor::get_num_nodes_image() */
        inline Index get_num_nodes_image() const
        {
          return _matrix.num_cols();
        }

        /** \copydoc Adjactor::image_begin() */
        inline ImageIterator image_begin(Index domain_node) const
        {
          const IT_ tbands(_matrix.num_bands());
          const Index trows(_matrix.num_rows());

          XASSERT(domain_node < trows);

          for(IT_ k(0); k < tbands; ++k)
          {
            if(_offsets[k] + domain_node + 1 >= trows)
              return ImageIterator(k, domain_node + 1 - trows, _offsets.get_r());
          }
        }

        /** \copydoc Adjactor::image_end() */
        inline ImageIterator image_end(Index domain_node) const
        {
          const IT_ tbands(_matrix.num_bands());
          const Index trows(_matrix.num_rows());

          XASSERT(domain_node < trows);

          for(IT_ k(0); k < tbands; ++k)
          {
            if(_offsets[k] + domain_node + 1 >= 2 * trows)
            {
              return ImageIterator(k, domain_node + 1 - trows, _offsets.get_r());
            }
          }
          return ImageIterator(tbands, domain_node + 1 - trows, _offsets.get_r());
        }
      }; // class Adjactor

      /**
       * \brief Scatter-Axpy operation for SparseMatrixBanded
       *
       * \author Christoph Lohmann
       */
      class ScatterAxpy
      {
      public:
        typedef LAFEM::SparseMatrixBanded<DT_, IT_> MatrixType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
#ifdef DEBUG
        static constexpr IT_ _deadcode = ~IT_(0);
#endif
        Index _num_rows;
        Index _num_cols;
        Index _num_bands;
        Memory::Arbiter _col_ptr_arbiter;
        const Memory::TypedView<IT_> _offsets;
        Memory::TypedView<IT_> _col_ptr;
        Memory::TypedView<DT_> _data;

      public:
        explicit ScatterAxpy(MatrixType& matrix) :
          _num_rows(matrix.num_rows()),
          _num_cols(matrix.num_cols()),
          _num_bands(matrix.num_bands()),
          _col_ptr_arbiter(_num_cols * sizeof(IT_), Memory::Location::main, Memory::Init::format_to_one),
          _offsets(matrix.offsets_view_r()),
          _col_ptr(_col_ptr_arbiter.view(Memory::Location::main, Memory::Access::read_write)),
          _data(matrix.val_view(Memory::Location::main, Memory::Access::read_write | Memory::Access::overlap))
        {
        }

        template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
        void operator()(const LocalMatrix_& loc_mat, const RowMapping_& row_map,
          const ColMapping_& col_map, DT_ alpha = DT_(1))
        {
          // loop over all local row entries
          for(int i(0); i < row_map.get_num_local_dofs(); ++i)
          {
            // fetch row index
            const Index ix = row_map.get_index(i);

            // build column pointer for this row entry contribution
            for(IT_ k(0); k < _num_bands; ++k)
            {
              if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
              {
                _col_ptr[_offsets[k] + ix + 1 - _num_rows] = IT_(k * _num_rows + ix);
              }
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // fetch column index
              const Index jx = col_map.get_index(j);

              // ensure that the column pointer is valid for this index
              ASSERTM(_col_ptr[jx] != _deadcode, "invalid column index");

              // incorporate data into global matrix
              _data[_col_ptr[jx]] += alpha * loc_mat[i][j];

              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(IT_ k(0); k < _num_bands; ++k)
            {
              if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
              {
                _col_ptr[_offsets[k] + ix + 1 - _num_rows] = _deadcode;
              }
            }
#endif
            // continue with next row entry
          }
        }
      }; // class ScatterAxpy

      /**
       * \brief Gather-Axpy operation for SparseMatrixBanded
       *
       * \author Christoph Lohmann
       */
      class GatherAxpy
      {
      public:
        typedef LAFEM::SparseMatrixBanded<DT_, IT_> MatrixType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
#ifdef DEBUG
        static constexpr IT_ _deadcode = ~IT_(0);
#endif
        Index _num_rows;
        Index _num_cols;
        Index _num_bands;
        Memory::Arbiter _col_ptr_arbiter;
        const Memory::TypedView<IT_> _offsets;
        Memory::TypedView<IT_> _col_ptr;
        const Memory::TypedView<DT_> _data;

      public:
        explicit GatherAxpy(const MatrixType& matrix) :
          _num_rows(matrix.num_rows()),
          _num_cols(matrix.num_cols()),
          _num_bands(matrix.num_bands()),
          _col_ptr_arbiter(_num_cols * sizeof(IT_), Memory::Location::main, Memory::Init::format_to_one),
          _offsets(matrix.offsets_view_r()),
          _col_ptr(_col_ptr_arbiter.view(Memory::Location::main, Memory::Access::read_write)),
          _data(matrix.val_view_r())
        {
        }

        template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
        void operator()(LocalMatrix_& loc_mat, const RowMapping_& row_map,
          const ColMapping_& col_map, DT_ alpha = DT_(1))
        {
          // loop over all local row entries
          for(int i(0); i < row_map.get_num_local_dofs(); ++i)
          {
            // fetch row index
            const Index ix = row_map.get_index(i);

            // build column pointer for this row entry contribution
            for(IT_ k(0); k < _num_bands; ++k)
            {
              if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
              {
                _col_ptr[_offsets[k] + ix + 1 - _num_rows] = IT_(k * _num_rows + ix);
              }
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // fetch column index
              const Index jx = col_map.get_index(j);

              // ensure that the column pointer is valid for this index
              ASSERTM(_col_ptr[jx] != _deadcode, "invalid column index");

              // update local matrix data
              loc_mat[i][j] += alpha * _data[_col_ptr[jx]];

              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(IT_ k(0); k < _num_bands; ++k)
            {
              if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
              {
                _col_ptr[_offsets[k] + ix + 1 - _num_rows] = _deadcode;
              }
            }
#endif
            // continue with next row entry
          }
        }
      }; // class GatherAxpy
    }; // class SparseMatrixBanded
  } // namespace LAFEM
} // namespace FEAT
