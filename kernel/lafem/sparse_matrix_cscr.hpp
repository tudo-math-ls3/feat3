// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_CSCR_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_CSCR_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/scale_row_col.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/apply.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/lafem/arch/diagonal.hpp>
#include <kernel/lafem/arch/lumping.hpp>
#include <kernel/lafem/arch/row_norm.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/lafem/vector_mirror.hpp>

#include <fstream>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief CSCR based sparse matrix.
     *
     * \tparam Mem_ The \ref FEAT::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse compressed row format.\n
     * This format is very similar to the standard CSR format but stores explicit row numbers for each non empty row.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: column index per non zero element \n
     * _indices[1]: row start index (including matrix end index)\n
     * _indices[2]: row number of each non empty row\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_index[4]: non zero row count
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseMatrixCSCR : public Container<Mem_, DT_, IT_>
    {
    private:
      Index & _size()
      {
        return this->_scalar_index.at(0);
      }

      Index & _rows()
      {
        return this->_scalar_index.at(1);
      }

      Index & _columns()
      {
        return this->_scalar_index.at(2);
      }

      Index & _used_elements()
      {
        return this->_scalar_index.at(3);
      }

      Index & _used_rows()
      {
        return this->_scalar_index.at(4);
      }

    public:
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Compatible L-vector type
      typedef DenseVector<Mem_, DT_, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<Mem_, DT_, IT_> VectorTypeR;
      /// Our used layout type
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_cscr;
      /// ImageIterator typedef for Adjactor interface implementation
      typedef const IT_* ImageIterator;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = SparseMatrixCSCR<Mem2_, DT2_, IT2_>;

      /// this typedef lets you create a matrix container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit SparseMatrixCSCR() :
        Container<Mem_, DT_, IT_> (0)
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
      explicit SparseMatrixCSCR(Index rows_in, Index columns_in) :
        Container<Mem_, DT_, IT_> (rows_in * columns_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] used_elements_in The amount of non zero elements of the created matrix.
       * \param[in] used_rows_in The amount of non zero rows of the created matrix.
       *
       * Creates an empty (but allocated) matrix.
       *
       * \note The allocated memory will not be initialised.
       */
      explicit SparseMatrixCSCR(Index rows_in, Index columns_in, Index used_elements_in, Index used_rows_in) :
        Container<Mem_, DT_, IT_> (rows_in * columns_in)
      {
        XASSERT(rows_in != Index(0) && columns_in != Index(0));

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(used_elements_in);
        this->_scalar_index.push_back(used_rows_in);

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_rows() + 1));
        this->_indices_size.push_back(_used_rows() + 1);

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_rows()));
        this->_indices_size.push_back(_used_rows());

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
      }

      /**
       * \brief Constructor
       *
       * \param[in] layout_in The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixCSCR(const SparseLayout<Mem_, IT_, layout_id> & layout_in) :
        Container<Mem_, DT_, IT_> (layout_in._scalar_index.at(0))
      {
        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());

        for (auto i : this->_indices)
          MemoryPool<Mem_>::increase_memory(i);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
      }

      /**
       * \brief Constructor
       *
       * \param[in] other The source matrix.
       *
       * Creates a CSCR matrix based on the source matrix.
       */
      template <typename MT_>
      explicit SparseMatrixCSCR(const MT_ & other) :
        Container<Mem_, DT_, IT_>(other.size())
      {
        convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a CSCR matrix based on the source file.
       */
      explicit SparseMatrixCSCR(FileMode mode, String filename) :
        Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a CSCR matrix based on the source filestream.
       */
      explicit SparseMatrixCSCR(FileMode mode, std::istream& file) :
        Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, file);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] col_ind_in Vector with column indices.
       * \param[in] val_in Vector with non zero elements.
       * \param[in] row_ptr_in Vector with start indices of all rows into the val/col_ind arrays.
       * Note that this vector must also contain the end index of the last row and thus has a size of row_count + 1.
       * \param[in] row_numbers_in Vector with the row number of each used row.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixCSCR(const Index rows_in, const Index columns_in,
                               DenseVector<Mem_, IT_, IT_> & col_ind_in, DenseVector<Mem_, DT_, IT_> & val_in, DenseVector<Mem_, IT_, IT_> & row_ptr_in,
                               DenseVector<Mem_, IT_, IT_> & row_numbers_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        /// \todo maybe create empty matrix if col_ind and val and row_ptr inputs are all three empty
        XASSERT(col_ind_in.size() > 0);
        XASSERT(val_in.size() > 0);
        XASSERT(row_ptr_in.size() > 0);
        XASSERT(row_numbers_in.size() > 0);
        XASSERT(row_numbers_in.size() + 1 == row_ptr_in.size());

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(val_in.size());
        this->_scalar_index.push_back(row_numbers_in.size());

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(col_ind_in.elements());
        this->_indices_size.push_back(col_ind_in.size());
        this->_indices.push_back(row_ptr_in.elements());
        this->_indices_size.push_back(row_ptr_in.size());
        this->_indices.push_back(row_numbers_in.elements());
        this->_indices_size.push_back(row_numbers_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] csr The source csr matrix.
       * \param[in] non_zero_rows The mirror describing the rows to be selected
       *
       * Creates a matrix with selected rows from a given csr matrix.
       */
      explicit SparseMatrixCSCR(const SparseMatrixCSR<Mem_, DT_, IT_> & csr, const VectorMirror<Mem_, DT_, IT_> & non_zero_rows) :
        Container<Mem_, DT_, IT_>(csr.size())
      {
        XASSERT(non_zero_rows.num_indices() > 0);

        const IT_ * row_main(nullptr);
        const IT_ * indices(nullptr);
        IT_ * trow_main(nullptr);
        IT_ * tindices(nullptr);
        if (typeid(Mem_) == typeid(Mem::Main))
        {
          row_main = csr.row_ptr();
          indices = non_zero_rows.indices();
        }
        else
        {
          trow_main = new IT_[csr.rows() + 1];
          MemoryPool<Mem_>::download(trow_main, csr.row_ptr(), csr.rows() + 1);
          tindices = new IT_[non_zero_rows.num_indices()];
          MemoryPool<Mem_>::download(tindices, non_zero_rows.indices(), non_zero_rows.num_indices());
          row_main = trow_main;
          indices = tindices;
        }

        Index tused_elements(0);
        for (Index i(0) ; i < non_zero_rows.num_indices() ; ++i)
        {
          const Index row(indices[i]);
          tused_elements += row_main[row+1] - row_main[row];
        }

        this->_scalar_index.push_back(csr.rows());
        this->_scalar_index.push_back(csr.columns());
        this->_scalar_index.push_back(tused_elements);
        this->_scalar_index.push_back(non_zero_rows.num_indices());

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_rows() + 1));
        this->_indices_size.push_back(_used_rows() + 1);

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_rows()));
        this->_indices_size.push_back(_used_rows());

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());

        Index row_count(0);
        IT_ offset(0);
        MemoryPool<Mem_>::set_memory(this->row_ptr(), offset);
        for (Index i(0) ; i < non_zero_rows.num_indices() ; ++i)
        {
          const IT_ row(indices[i]);
          const IT_ start(row_main[row]);
          const IT_ end(row_main[row+1]);
          const IT_ row_size(end - start);
          MemoryPool<Mem_>::set_memory(this->row_ptr() + row_count + 1, offset + row_size);
          MemoryPool<Mem_>::set_memory(this->row_numbers() + row_count, row);
          MemoryPool<Mem_>::copy(this->val() + offset, csr.val() + start, row_size);
          MemoryPool<Mem_>::copy(this->col_ind() + offset, csr.col_ind() + start, row_size);
          ++row_count;
          offset += row_size;
        }

        delete[] trow_main;
        delete[] tindices;
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix from
       *
       * Creates a CSCR matrix based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit SparseMatrixCSCR(const Adjacency::Graph & graph) :
        Container<Mem_, DT_, IT_>(0)
      {
        // get number of rows, columns and indices
        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nzes = graph.get_num_indices();

        if (num_nzes == 0)
        {
          this->assign(SparseMatrixCSCR(num_rows, num_cols));
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
        SparseMatrixCSCR<Mem::Main, DT_, IT_> matrix(num_rows, num_cols, num_nzes, num_nzrs);

        // get matrix arrays
        IndexType* trow_ptr = matrix.row_ptr();
        IndexType* trow_idx = matrix.row_numbers();
        IndexType* tcol_idx = matrix.col_ind();

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

        this->convert(matrix);
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixCSCR(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        deserialise<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixCSCR(SparseMatrixCSCR && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixCSCR>(other))
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
      template<typename Mem2_, typename DT2_, typename IT2_>
      void clone(const SparseMatrixCSCR<Mem2_, DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Weak)
      {
        Container<Mem_, DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixCSCR<Mem2_, DT2_, IT2_> & other)
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
        typename MT_::template ContainerType<Mem::Main, DT_, IT_> ta;
        ta.convert(a);

        const Index arows(ta.template rows<Perspective::pod>());
        const Index acolumns(ta.template columns<Perspective::pod>());
        const Index aused_elements(ta.template used_elements<Perspective::pod>());
        Index aused_rows(0);

        for (Index i(0); i < arows; ++i)
        {
          aused_rows += ta.get_length_of_line(i) > 0;
        }

        DenseVector<Mem::Main, DT_, IT_> tval(aused_elements);
        DenseVector<Mem::Main, IT_, IT_> tcol_ind(aused_elements);
        DenseVector<Mem::Main, IT_, IT_> trow_ptr(aused_rows + 1);
        DenseVector<Mem::Main, IT_, IT_> trow_numbers(aused_rows);

        DT_ * pval(tval.elements());
        IT_ * pcol_ind(tcol_ind.elements());
        IT_ * prow_ptr(trow_ptr.elements());
        IT_ * prow_numbers(trow_numbers.elements());

        Index row_index(0);
        Index offset(0);
        prow_ptr[0] = IT_(0);
        for (Index i(0); i < arows; ++i)
        {
          Index length = ta.get_length_of_line(i);
          if (length > 0)
          {
            prow_ptr[row_index + 1] = IT_(offset + length);
            prow_numbers[row_index] = IT_(i);
            ++row_index;
            offset += length;
          }
        }


        for (Index i(0); i < arows; ++i)
        {
          ta.set_line(i, pval + prow_ptr[i], pcol_ind + prow_ptr[i], 0);
        }

        SparseMatrixCSCR<Mem::Main, DT_, IT_> ta_cscr(arows, acolumns, tcol_ind, tval, trow_ptr, trow_numbers);
        SparseMatrixCSCR<Mem_, DT_, IT_> a_cscr;
        a_cscr.convert(ta_cscr);

        this->assign(a_cscr);
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixCSCR & operator= (const SparseLayout<Mem_, IT_, layout_id> & layout_in)
      {
        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_indices.at(i));

        this->_elements.clear();
        this->_indices.clear();
        this->_elements_size.clear();
        this->_indices_size.clear();
        this->_scalar_index.clear();

        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());

        for (auto i : this->_indices)
          MemoryPool<Mem_>::increase_memory(i);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());

        return *this;
      }

      /**
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialise(std::vector<char> input)
      {
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_cscr, input);
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * \param[in] config LAFEM::SerialConfig, a struct describing the serialise configuration.
       * \note the corresponding configure flags 'zlib' and/or 'zfp' need to be added in the build-id at the configure call.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialise for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialise(const LAFEM::SerialConfig& config = SerialConfig())
      {
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_cscr, config);
      }

      /**
       * \brief Read in matrix from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, String filename)
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
        case FileMode::fm_binary:
          this->template _deserialise<double, std::uint64_t>(FileMode::fm_cscr, file);
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
      void write_out(FileMode mode, String filename) const
      {
        std::ios_base::openmode bin = std::ofstream::out | std::ofstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ofstream::out;
        std::ofstream file(filename.c_str(), bin);
        if (! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);
        write_out(mode, file);
        file.close();
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
        case FileMode::fm_binary:
          this->template _serialise<double, std::uint64_t>(FileMode::fm_cscr, file);
          break;
        default:
          XABORTM("Filemode not supported!");
        }
      }

      /**
       * \brief Retrieve specific matrix element.
       *
       * \param[in] row The row of the matrix element.
       * \param[in] col The column of the matrix element.
       *
       * \returns Specific matrix element.
       */
      DT_ operator()(Index row, Index col) const
      {
        ASSERT(row < rows());
        ASSERT(col < columns());

        Index row_index(0);
        while (row_index < this->used_rows() && Index(MemoryPool<Mem_>::get_element(this->row_numbers(), row_index)) < row)
        {
          ++row_index;
        }

        if (row_index == this->used_rows() || Index(MemoryPool<Mem_>::get_element(this->row_numbers(), row_index)) > row)
          return DT_(0.);

        //row_numbers[row_index] == row
        for (Index i(Index(MemoryPool<Mem_>::get_element(this->row_ptr(), row_index))) ; i < Index(MemoryPool<Mem_>::get_element(this->row_ptr(), row_index + 1)) ; ++i)
        {
          if (Index(MemoryPool<Mem_>::get_element(this->_indices.at(0), i)) == col)
            return MemoryPool<Mem_>::get_element(this->_elements.at(0), i);
          if (Index(MemoryPool<Mem_>::get_element(this->_indices.at(0), i)) > col)
            return DT_(0.);
        }
        return DT_(0.);
      }

      /**
       * \brief Retrieve convenient sparse matrix layout object.
       *
       * \return An object containing the sparse matrix layout.
       */
      SparseLayout<Mem_, IT_, layout_id> layout() const
      {
        return SparseLayout<Mem_, IT_, layout_id>(this->_indices, this->_indices_size, this->_scalar_index);
      }

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count.
       */
      template <Perspective = Perspective::native>
      Index rows() const
      {
        return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      template <Perspective = Perspective::native>
      Index columns() const
      {
        return this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve used matrix non zero row count.
       *
       * \returns Used matrix non zero row count.
       */
      template <Perspective = Perspective::native>
      Index used_rows() const
      {
        return this->_scalar_index.at(4);
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      template <Perspective = Perspective::native>
      Index used_elements() const
      {
        return this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve column indices array.
       *
       * \returns Column indices array.
       */
      IT_ * col_ind()
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(0);
      }

      IT_ const * col_ind() const
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(0);
      }

      /**
       * \brief Retrieve non zero element array.
       *
       * \returns Non zero element array.
       */
      DT_ * val()
      {
        if (this->_elements.size() == 0)
          return nullptr;

        return this->_elements.at(0);
      }

      DT_ const * val() const
      {
        if (this->_elements.size() == 0)
          return nullptr;

        return this->_elements.at(0);
      }

      /**
       * \brief Retrieve row start index array.
       *
       * \returns Row start index array.
       */
      IT_ * row_ptr()
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(1);
      }

      IT_ const * row_ptr() const
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(1);
      }

      /**
       * \brief Retrieve row numbers of non zero rows.
       *
       * \returns Row number array.
       */
      IT_ * row_numbers()
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(2);
      }

      IT_ const * row_numbers() const
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(2);
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

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      template <typename Mem2_>
      void copy(const SparseMatrixCSCR<Mem2_, DT_, IT_> & x, bool full = false)
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
      void axpy(
                const SparseMatrixCSCR & x,
                const SparseMatrixCSCR & y,
                const DT_ alpha = DT_(1))
      {
        XASSERTM(x.rows() == y.rows(), "Matrix rows do not match!");
        XASSERTM(x.rows() == this->rows(), "Matrix rows do not match!");
        XASSERTM(x.columns() == y.columns(), "Matrix columns do not match!");
        XASSERTM(x.columns() == this->columns(), "Matrix columns do not match!");
        XASSERTM(x.used_elements() == y.used_elements(), "Matrix used_elements do not match!");
        XASSERTM(x.used_elements() == this->used_elements(), "Matrix used_elements do not match!");

        if (Math::abs(alpha) < Math::eps<DT_>())
        {
          this->copy(y);
          //y.scale(beta);
          return;
        }

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements() * 2);
        Arch::Axpy<Mem_>::dv(this->val(), alpha, x.val(), y.val(), this->used_elements());

        TimeStamp ts_stop;
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
        XASSERTM(x.rows() == this->rows(), "Row count does not match!");
        XASSERTM(x.columns() == this->columns(), "Column count does not match!");
        XASSERTM(x.used_elements() == this->used_elements(), "Nonzero count does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements());
        Arch::Scale<Mem_>::value(this->val(), x.val(), alpha, this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements() * 2);
        DT_ result = Arch::Norm2<Mem_>::value(this->val(), this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] transposed Should the product use the transposed matrix?
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x, bool transposed = false) const
      {
        if (transposed)
        {
          XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
          XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");
        }
        else
        {
          XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
          XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        }

        TimeStamp ts_start;

        if (this->used_elements() == 0)
        {
          r.format();
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops(this->used_elements() * 2);
        Arch::Apply<Mem_>::cscr(r.elements(), DT_(1), x.elements(), DT_(0), r.elements(),
            this->val(), this->col_ind(), this->row_ptr(), this->row_numbers(), this->used_rows(), this->rows(), this->columns(), this->used_elements(), transposed);

        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       * \param[in] transposed Should the product use the transposed matrix?
       */
      void apply(
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ alpha = DT_(1),
                 const bool transposed = false) const
      {
        if (transposed)
        {
          XASSERTM(r.size() == this->columns(), "Vector size of r does not match!");
          XASSERTM(x.size() == this->rows(), "Vector size of x does not match!");
          XASSERTM(y.size() == this->columns(), "Vector size of y does not match!");
        }
        else
        {
          XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
          XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
          XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");
        }

        TimeStamp ts_start;

        if (this->used_elements() == 0 || Math::abs(alpha) < Math::eps<DT_>())
        {
          r.copy(y);
          //r.scale(beta);
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops( (this->used_elements() + this->rows()) * 2 );
        Arch::Apply<Mem_>::cscr(r.elements(), alpha, x.elements(), DT_(1.), y.elements(),
            this->val(), this->col_ind(), this->row_ptr(), this->row_numbers(), this->used_rows(), this->rows(), this->columns(), this->used_elements(), transposed);

        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }
      ///@}

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(this->rows());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->columns());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        const IT_ * prow_ptr(this->row_ptr());
        return Index(prow_ptr[row + 1] - prow_ptr[row]);
      }

      /// \cond internal

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DT_ * const pval_set, IT_ * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        const IT_ * prow_ptr(this->row_ptr());
        const IT_ * pcol_ind(this->col_ind());
        const DT_ * pval(this->val());

        Index row_idx(0);
        while (row_idx < this->used_rows() && row_idx < row)
        {
          ++row_idx;
        }
        XASSERTM(row_idx < this->used_rows(), "invalid row number provided!");

        const Index start((Index(prow_ptr[row_idx])));
        const Index end((Index(prow_ptr[row_idx + 1] - prow_ptr[row_idx])));
        for (Index i(0); i < end; ++i)
        {
          pval_set[i * stride] = pval[start + i];
          pcol_set[i * stride] = pcol_ind[start + i] + IT_(col_start);
        }
      }

      /// \endcond

      /**
       * \brief SparseMatrixCSCR comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_> friend bool operator== (const SparseMatrixCSCR & a, const SparseMatrixCSCR<Mem2_, DT_, IT_> & b)
      {
        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;
        if (a.used_elements() != b.used_elements())
          return false;
        if (a.used_rows() != b.used_rows())
          return false;

        if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && a.get_indices().size() == 0 && b.get_elements().size() == 0 && b.get_indices().size() == 0)
          return true;

        IT_ * col_ind_a;
        IT_ * col_ind_b;
        DT_ * val_a;
        DT_ * val_b;
        IT_ * row_ptr_a;
        IT_ * row_ptr_b;
        IT_ * row_numbers_a;
        IT_ * row_numbers_b;

        bool ret(true);

        if(std::is_same<Mem::Main, Mem_>::value)
        {
          col_ind_a = const_cast<IT_*>(a.col_ind());
          val_a = const_cast<DT_*>(a.val());
          row_ptr_a = const_cast<IT_*>(a.row_ptr());
          row_numbers_a = const_cast<IT_*>(a.row_numbers());
        }
        else
        {
          col_ind_a = new IT_[a.used_elements()];
          MemoryPool<Mem_>::template download<IT_>(col_ind_a, a.col_ind(), a.used_elements());
          val_a = new DT_[a.used_elements()];
          MemoryPool<Mem_>::template download<DT_>(val_a, a.val(), a.used_elements());
          row_ptr_a = new IT_[a.rows() + 1];
          MemoryPool<Mem_>::template download<IT_>(row_ptr_a, a.row_ptr(), a.used_rows() + 1);
          row_numbers_a = new IT_[a.rows()];
          MemoryPool<Mem_>::template download<IT_>(row_numbers_a, a.row_numbers(), a.used_rows());
        }
        if(std::is_same<Mem::Main, Mem2_>::value)
        {
          col_ind_b = const_cast<IT_*>(b.col_ind());
          val_b = const_cast<DT_*>(b.val());
          row_ptr_b = const_cast<IT_*>(b.row_ptr());
          row_numbers_b = const_cast<IT_*>(b.row_numbers());
        }
        else
        {
          col_ind_b = new IT_[b.used_elements()];
          MemoryPool<Mem2_>::template download<IT_>(col_ind_b, b.col_ind(), b.used_elements());
          val_b = new DT_[b.used_elements()];
          MemoryPool<Mem2_>::template download<DT_>(val_b, b.val(), b.used_elements());
          row_ptr_b = new IT_[b.rows() + 1];
          MemoryPool<Mem2_>::template download<IT_>(row_ptr_b, b.row_ptr(), b.used_rows() + 1);
          row_numbers_b = new IT_[b.rows()];
          MemoryPool<Mem2_>::template download<IT_>(row_numbers_b, b.row_numbers(), b.used_rows());
        }

        for (Index i(0) ; i < a.used_elements() ; ++i)
        {
          if (col_ind_a[i] != col_ind_b[i])
          {
            ret = false;
            break;
          }
          if (val_a[i] != val_b[i])
          {
            ret = false;
            break;
          }
        }
        if (ret)
        {
          for (Index i(0) ; i < a.used_rows() + 1; ++i)
          {
            if (row_ptr_a[i] != row_ptr_b[i])
            {
              ret = false;
              break;
            }
          }
        }
        if (ret)
        {
          for (Index i(0) ; i < a.used_rows(); ++i)
          {
            if (row_numbers_a[i] != row_numbers_b[i])
            {
              ret = false;
              break;
            }
          }
        }

        if(! std::is_same<Mem::Main, Mem_>::value)
        {
          delete[] col_ind_a;
          delete[] val_a;
          delete[] row_ptr_a;
          delete[] row_numbers_a;
        }
        if(! std::is_same<Mem::Main, Mem2_>::value)
        {
          delete[] col_ind_b;
          delete[] val_b;
          delete[] row_ptr_b;
          delete[] row_numbers_b;
        }

        return ret;
      }

      /**
       * \brief SparseMatrixCSCR streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const SparseMatrixCSCR & b)
      {

        lhs << "[" << std::endl;
        for (Index i(0) ; i < b.rows() ; ++i)
        {
          lhs << "[";
          for (Index j(0) ; j < b.columns() ; ++j)
          {
            lhs << "  " << b(i, j);
          }
          lhs << "]" << std::endl;
        }
        lhs << "]" << std::endl;

        return lhs;
      }
    }; //SparseMatrixCSCR

#ifdef FEAT_EICKT
    extern template class SparseMatrixCSCR<Mem::Main, float, unsigned int>;
    extern template class SparseMatrixCSCR<Mem::Main, double, unsigned int>;
    extern template class SparseMatrixCSCR<Mem::Main, float, unsigned long>;
    extern template class SparseMatrixCSCR<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
    extern template class SparseMatrixCSCR<Mem::CUDA, float, unsigned int>;
    extern template class SparseMatrixCSCR<Mem::CUDA, double, unsigned int>;
    extern template class SparseMatrixCSCR<Mem::CUDA, float, unsigned long>;
    extern template class SparseMatrixCSCR<Mem::CUDA, double, unsigned long>;
#endif
#endif

  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_SPARSE_MATRIX_CSCR_HPP
