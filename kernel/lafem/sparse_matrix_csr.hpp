// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/likwid_marker.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/axpy_dense.hpp>
#include <kernel/lafem/arch/diagonal_csr.hpp>
#include <kernel/lafem/arch/lumping_csr.hpp>
#include <kernel/lafem/arch/matvecmult_csr_block.hpp>
#include <kernel/lafem/arch/matvecmult_csr_dense.hpp>
#include <kernel/lafem/arch/norm2_dense.hpp>
#include <kernel/lafem/arch/row_norm2_csr.hpp>
#include <kernel/lafem/arch/scale_col_csr.hpp>
#include <kernel/lafem/arch/scale_dense.hpp>
#include <kernel/lafem/arch/scale_row_csr.hpp>

// includes, system
#include <fstream>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief CSR based sparse matrix.
     *
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse row format.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: row start index (including matrix end index)\n
     * _indices[1]: column index per non zero element \n
     *
     * _scalar_index[0]: row count \n
     * _scalar_index[1]: column count \n
     * _scalar_index[2]: non zero element count (used elements) \n
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename DT_, typename IT_ = Index>
    class SparseMatrixCSR : public Container<DT_, IT_>
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
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_csr;
      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = SparseMatrixCSR<DT2_, IT2_>;

      /// this typedef lets you create a matrix container with different Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

      /// this is not a global container
      static constexpr bool is_global = false;

      /// this is a local container
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

      Index & _nzes()
      {
        return this->_scalar_index.at(2);
      }

    public:
      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      SparseMatrixCSR()
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
       *
       * Creates an empty matrix.
       * Because SparseMatrixCSR is a read-only container, it stays empty.
       *
       * \note This matrix does not allocate any memory
       */
      explicit SparseMatrixCSR(Index rows_in, Index columns_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(0);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] nonzeros_in The amount of non zero elements of the created matrix.
       *
       * Creates an empty (but allocated) matrix.
       *
       * \note The allocated memory will not be initialized.
       */
      explicit SparseMatrixCSR(Index rows_in, Index columns_in, Index nonzeros_in)
      {
        XASSERT(rows_in > Index(0));
        XASSERT(columns_in > Index(0));

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(nonzeros_in);

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * (rows_in + 1)));
        this->_indices_size.push_back(rows_in + 1);

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * nonzeros_in));
        this->_indices_size.push_back(nonzeros_in);

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * nonzeros_in));
        this->_elements_size.push_back(nonzeros_in);
      }

      /**
       * \brief Constructor
       *
       * \param[in] layout_in The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixCSR(const SparseLayout<IT_, layout_id> & layout_in)
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
       * Creates a CSR matrix based on the source matrix.
       */
      template <typename MT_>
      explicit SparseMatrixCSR(const MT_ & other)
      {
        this->convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix from
       *
       * Creates a CSR matrix based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit SparseMatrixCSR(const Adjacency::Graph & graph) :
        SparseMatrixCSR(graph.get_num_nodes_domain(), graph.get_num_nodes_image(), graph.get_num_indices())
      {
        const Index num_nnze = graph.get_num_indices();
        const Index num_rows = graph.get_num_nodes_domain();
        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());

        Memory::TypedView<IT_> prow_ptr = this->row_ptr_view_w();
        Memory::TypedView<IT_> pcol_idx = this->col_idx_view_w();

        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i <= num_rows; ++i)
          prow_ptr[i] = IT_(dom_ptr[i]);

        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < num_nnze; ++i)
          pcol_idx[i] = IT_(img_idx[i]);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a CSR matrix based on the source file.
       */
      explicit SparseMatrixCSR(FileMode mode, const String& filename) :
        Container<DT_, IT_>()
      {
        this->read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a CSR matrix based on the source filestream.
       */
      explicit SparseMatrixCSR(FileMode mode, std::istream& file) :
        Container<DT_, IT_>()
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
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixCSR(const Index rows_in, const Index columns_in,
        DenseVector<IT_, IT_> & row_ptr_in, DenseVector<IT_, IT_> & col_idx_in, DenseVector<DT_, IT_> & val_in)
      {
        /// \todo maybe create empty matrix if col_idx and val and row_ptr inputs are all three empty
        XASSERT(row_ptr_in.size() > 0);
        XASSERT(col_idx_in.size() > 0);
        XASSERT(col_idx_in.size() == val_in.size());

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(val_in.size());

        this->_indices.push_back(row_ptr_in.elements_arbiter().attach());
        this->_indices_size.push_back(row_ptr_in.size());
        this->_indices.push_back(col_idx_in.elements_arbiter().attach());
        this->_indices_size.push_back(col_idx_in.size());
        this->_elements.push_back(val_in.elements_arbiter().attach());
        this->_elements_size.push_back(val_in.size());
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixCSR(std::vector<char> input)
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
      SparseMatrixCSR(SparseMatrixCSR && other) :
        Container<DT_, IT_>(std::forward<SparseMatrixCSR>(other))
      {
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixCSR & operator= (SparseMatrixCSR && other)
      {
        this->move(std::forward<SparseMatrixCSR>(other));

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
      SparseMatrixCSR clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        SparseMatrixCSR t;
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
      void clone(const SparseMatrixCSR<DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Weak)
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
      void convert(const SparseMatrixCSR<DT2_, IT2_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] graph The graph to create the matrix from
       *
       * Creates a CSR matrix based on a given adjacency graph, representing the sparsity pattern.
       */
      void convert(const Adjacency::Graph & graph)
      {
        this->move(SparseMatrixCSR(graph));
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
        this->clear();

        this->_scalar_index.push_back(other.num_rows());
        this->_scalar_index.push_back(other.num_cols());
        this->_scalar_index.push_back(other.num_nzes());

        if (other.num_nzes() == 0)
          return;

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * (_rows() + 1)));
        this->_indices_size.push_back(_rows() + 1);
        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * _nzes()));
        this->_indices_size.push_back(_nzes());
        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * _nzes()));
        this->_elements_size.push_back(_nzes());

        Memory::TypedView<IT_> row_ptr = this->row_ptr_view_w();
        Memory::TypedView<IT_> col_idx = this->col_idx_view_w();
        Memory::TypedView<DT_> val = this->val_view_w();

        const Memory::TypedView<DT_> cval = other.val_view_r();
        const Memory::TypedView<IT_> coffsets = other.offsets_view_r();
        const Index nbands(other.num_bands());
        const Index nrows(other.num_rows());
        const Index ncols(other.num_cols());
        //const Index nnzes(other.num_nzes());

        row_ptr[0] = 0;

        // Search first offset of the upper triangular matrix
        Index k(0);
        while (k < nbands && coffsets[k] + 1 < nrows)
        {
          ++k;
        }

        IT_ ue(0);
        // iteration over all offsets of the lower triangular matrix
        for (Index i(k + 1); i > 0;)
        {
          --i;

          // iteration over all offsets of the upper triangular matrix
          for (Index j(nbands + 1); j > 0;)
          {
            --j;

            // iteration over all rows which contain the offsets between offset i and offset j
            Index so1 = Index(0);
            Index so2 = nrows;
            Index eo1 = Index(0);
            Index eo2 = nrows;
            if(i < nbands)
              so1 = Math::max(ncols + Index(1), nrows + ncols - Index(coffsets[i])) - ncols - Index(1);
            if(j < nbands)
              eo1 = Math::min(nrows, ncols + nrows - Index(coffsets[j]) - Index(1));
            if(i > 0)
              so2 = Math::max(ncols + Index(1), nrows + ncols - Index(coffsets[i-1])) - ncols - Index(1);
            if(j > 0)
              eo2 = Math::min(nrows, ncols + nrows - Index(coffsets[j-1]) - Index(1));

            const Index start = Math::max(so1, eo1);
            const Index end   = Math::min(so2, eo2);

            for (Index l(start); l < end; ++l)
            {
              for (Index a(i); a < j; ++a)
              {
                val[ue] = cval[a * nrows + l];
                col_idx[ue] = IT_(l + coffsets[a] + 1 - nrows);
                ++ue;
              }
              row_ptr[l + 1] = ue;
            }
          }
        }
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename DT2_, typename IT2_, int block_height_, int block_width_>
      void convert(const SparseMatrixBCSR<DT2_, IT2_, block_height_, block_width_> & other)
      {
        this->clear();

        this->_scalar_index.push_back(other.num_rows_raw());
        this->_scalar_index.push_back(other.num_cols_raw());
        this->_scalar_index.push_back(other.num_nzes_raw());

        if (other.num_nzes_raw() == 0)
          return;

        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * (_rows() + 1)));
        this->_indices_size.push_back(_rows() + 1);
        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * _nzes()));
        this->_indices_size.push_back(_nzes());
        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * _nzes()));
        this->_elements_size.push_back(_nzes());

        Memory::TypedView<IT_> row_ptr = this->row_ptr_view_w();
        Memory::TypedView<IT_> col_idx = this->col_idx_view_w();
        Memory::TypedView<DT_> val = this->val_view_w();

        const Memory::TypedView<IT_> brow_ptr = other.row_ptr_view_r();
        const Memory::TypedView<IT_> bcol_idx = other.col_idx_view_r();
        const Memory::TypedView<Tiny::Matrix<DT_, block_height_, block_width_>> bval = other.val_view_r();
        const Index brows =  other.num_rows();

        row_ptr[0] = IT_(0);
        for (Index orow(0), ait(0) ; orow < brows ; ++orow)
        {
          for (int row(0) ; row < block_height_ ; ++row)
          {
            for(Index ocol(brow_ptr[orow]) ; ocol < brow_ptr[orow + 1] ; ++ocol)
            {
              for (int col(0) ; col < block_width_ ; ++col, ++ait)
              {
                val[ait] = bval[ocol](row,col);
                col_idx[ait] = bcol_idx[ocol] * IT_(block_width_) + IT_(col);
              }
            }
            row_ptr[orow * Index(block_height_) + Index(row) + 1] = IT_(ait);
          }
        }
      }

      /**
       * \brief Conversion method
       *
       * \param[in] a The input matrix.
       *
       * Converts any matrix to SparseMatrixCSR-format
       */
      template <typename MT_>
      void convert(const MT_ & a)
      {
        XASSERT(a.num_nzes_raw() > 0);

        const Index arows(a.num_rows_raw());
        const Index acols(a.num_cols_raw());
        const Index anzes(a.num_nzes_raw());

        SparseMatrixCSR<DT_, IT_> csr(arows, acols, anzes);

        Memory::TypedView<IT_> row_ptr = csr.row_ptr_view_w();
        Memory::TypedView<IT_> col_idx = csr.col_idx_view_w();
        Memory::TypedView<DT_> val = csr.val_view_w();

        row_ptr[0] = IT_(0);
        for (Index i(0); i < arows; ++i)
        {
          row_ptr[i + 1] = row_ptr[i] + IT_(a.row_degree(i));
        }

        for (Index i(0); i < arows; ++i)
        {
          a.get_row_col_indices(i, &col_idx[row_ptr[i]], IT_(0));
          a.get_row_values(i, &val[row_ptr[i]]);
        }

        val.release();
        col_idx.release();
        row_ptr.release();

        this->move(std::move(csr));
      }

      /**
      * \brief Copies the values of the input matrix into this matrix
      *
      * \param[in] source
      * A \transient reference to the source matrix.
      *
      * \attention
      * This function silently assumes that this matrix was created by using the convert function
      * from the source matrix, therefore assuring that the sparsity patterns are identical.
      */
      template <typename MT_>
      void copy(const MT_ & source)
      {
        XASSERT(source.num_rows_raw() == this->num_rows());
        XASSERT(source.num_cols_raw() == this->num_cols());
        XASSERT(source.num_nzes_raw() == this->num_nzes());

        const Memory::TypedView<IT_> row_ptr = this->row_ptr_view_r();
        Memory::TypedView<DT_> val = this->val_view_w();

        const Index nrows = this->num_rows();
        for (Index i(0); i < nrows; ++i)
        {
          source.get_row_values(i, &val[row_ptr[i]]);
        }
      }

      /**
       * \brief Copies the values of this matrix into the target matrix
       *
       * \param[out] target
       * A \transient reference to the target matrix.
       *
       * \attention
       * This function silently assumes that this matrix was created by using the convert function
       * from the target matrix, therefore assuring that the sparsity patterns are identical.
       */
      template <typename MT_>
      void copy_to(MT_ & target) const
      {
        XASSERT(target.num_rows_raw() == this->num_rows());
        XASSERT(target.num_cols_raw() == this->num_cols());
        XASSERT(target.num_nzes_raw() == this->num_nzes());

        const Memory::TypedView<IT_> row_ptr = this->row_ptr_view_r();
        const Memory::TypedView<DT_> val = this->val_view_r();

        const Index nrows = this->num_rows();
        for (Index i(0); i < nrows; ++i)
        {
          target.set_row_values(i, &val[row_ptr[i]]);
        }
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixCSR & operator= (const SparseLayout<IT_, layout_id> & layout_in)
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
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_csr, input);
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
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_csr, config);
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
        case FileMode::fm_mtx:
          {
            this->clear();
            this->_scalar_index.push_back(0);
            this->_scalar_index.push_back(0);
            this->_scalar_index.push_back(0);

            std::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >

            Index ue(0);
            String line;
            std::getline(file, line);
            const bool general((line.find("%%MatrixMarket matrix coordinate real general") != String::npos) ? true : false);
            const bool symmetric((line.find("%%MatrixMarket matrix coordinate real symmetric") != String::npos) ? true : false);

            if (symmetric == false && general == false)
            {
              XABORTM("Input-file is not a compatible mtx-file");
            }

            while(!file.eof())
            {
              std::getline(file,line);
              if (file.eof())
                XABORTM("Input-file is empty");

              String::size_type begin(line.find_first_not_of(" "));
              if (line.at(begin) != '%')
                break;
            }
            {
              String::size_type begin(line.find_first_not_of(" "));
              line.erase(0, begin);
              String::size_type end(line.find_first_of(" "));
              String srow(line, 0, end);
              Index row((Index)atol(srow.c_str()));
              line.erase(0, end);

              begin = line.find_first_not_of(" ");
              line.erase(0, begin);
              end = line.find_first_of(" ");
              String scol(line, 0, end);
              Index col((Index)atol(scol.c_str()));
              line.erase(0, end);
              _rows() = row;
              _cols() = col;
            }

            while(!file.eof())
            {
              std::getline(file, line);
              if (file.eof())
                break;

              String::size_type begin(line.find_first_not_of(" "));
              line.erase(0, begin);
              String::size_type end(line.find_first_of(" "));
              String srow(line, 0, end);
              IT_ row((IT_)atol(srow.c_str()));
              --row;
              line.erase(0, end);

              begin = line.find_first_not_of(" ");
              line.erase(0, begin);
              end = line.find_first_of(" ");
              String scol(line, 0, end);
              IT_ col((IT_)atol(scol.c_str()));
              --col;
              line.erase(0, end);

              begin = line.find_first_not_of(" ");
              line.erase(0, begin);
              end = line.find_first_of(" ");
              String sval(line, 0, end);
              DT_ tval((DT_)atof(sval.c_str()));

              entries[IT_(row)].insert(std::pair<IT_, DT_>(col, tval));
              ++ue;
              if (symmetric == true && row != col)
              {
                entries[IT_(col)].insert(std::pair<IT_, DT_>(row, tval));
                ++ue;
              }
            }
            _nzes() = ue;

            this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * (_rows() + 1)));
            this->_indices_size.push_back(num_rows() + 1);
            this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * _nzes()));
            this->_indices_size.push_back(_nzes());
            this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * _nzes()));
            this->_elements_size.push_back(_nzes());

            Memory::TypedView<IT_> trow_ptr(this->row_ptr_view_w());
            Memory::TypedView<IT_> tcol_idx(this->col_idx_view_w());
            Memory::TypedView<DT_> tval(this->val_view_w());

            IT_ idx(0);
            Index row_idx(0);
            for (auto row : entries)
            {
              trow_ptr[row_idx] = idx;
              for (auto col : row.second)
              {
                tcol_idx[idx] = col.first;
                tval[idx] = col.second;
                ++idx;
              }
              row.second.clear();
              ++row_idx;
            }
            trow_ptr[num_rows()] = IT_(ue);
            entries.clear();

            break;
          }

        case FileMode::fm_csr:
          [[fallthrough]];

        case FileMode::fm_binary:
          this->template _deserialize<double, std::uint64_t>(FileMode::fm_csr, file);
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
       * \param[in] symmetric Should we store only the lower half of the matrix in symmetric format?
       *
       * \warning This routine does no check on symmetric properties of the source matrix!
       * \warning Only uses symmetric properties if mode is mtx!
       */
      void write_out(FileMode mode, const String& filename, bool symmetric = false) const
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
        write_out(mode, file, symmetric);
        file.close();
        delete[] buff;
      }

      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be written to.
       * \param[in] symmetric Should we store only the lower half of the matrix in symmetric format?
       *
       * \warning This routine does no check on symmetric properties of the source matrix!
       * \warning Only uses symmetric properties if mode is mtx!
       */
      void write_out(FileMode mode, std::ostream& file, bool symmetric = false) const
      {
        switch(mode)
        {
        case FileMode::fm_csr:
          [[fallthrough]];

        case FileMode::fm_binary:
          this->template _serialize<double, std::uint64_t>(FileMode::fm_csr, file);
          break;

        case FileMode::fm_mtx:
          {
            const Memory::TypedView<IT_> trow_ptr = this->row_ptr_view_r();
            const Memory::TypedView<IT_> tcol_idx = this->col_idx_view_r();
            const Memory::TypedView<DT_> tval = this->val_view_r();

            if (symmetric)
            {
              file << "%%MatrixMarket matrix coordinate real symmetric\n";
              std::vector<IT_> rowv;
              std::vector<IT_> colv;
              std::vector<DT_> valv;

              for (Index row(0) ; row < num_rows() ; ++row)
              {
                const IT_ end(trow_ptr[row + 1]);
                for (IT_ i(trow_ptr[row]) ; i < end ; ++i)
                {
                  const IT_ col(tcol_idx[i]);
                  if (row >= col)
                  {
                    rowv.push_back(IT_(row + 1));
                    colv.push_back(col + 1);
                    valv.push_back(tval[i]);
                  }
                }
              }
              file << this->num_rows() << " " << this->num_cols() << " " << valv.size() << "\n";
              for (Index i(0) ; i < valv.size() ; ++i)
              {
                file << rowv.at(i) << " " << colv.at(i) << " " << valv.at(i) << "\n";
              }
            }
            else
            {
              file << "%%MatrixMarket matrix coordinate real general\n";
              file << this->num_rows() << " " << this->num_cols() << " " << this->num_nzes() << "\n";

              for (Index row(0) ; row < num_rows() ; ++row)
              {
                const IT_ end(trow_ptr[row + 1]);
                for (IT_ i(trow_ptr[row]) ; i < end ; ++i)
                {
                  file << row + 1 << " " << tcol_idx[i] + 1 << " " << tval[i] << "\n";
                }
              }
            }
            break;
          }

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
        return this->_scalar_index.empty() || (this->_scalar_index.at(2) == Index(0));
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
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      Index num_nzes() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(2);
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

      /// Returns a reference to the column-index array arbiter
      Memory::Arbiter& col_idx_arbiter()
      {
        return this->_indices.at(1);
      }

      /// Returns a reference to the column-index array arbiter
      const Memory::Arbiter& col_idx_arbiter() const
      {
        return this->_indices.at(1);
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

      Memory::TypedView<IT_> col_idx_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::read));
      }

      Memory::TypedView<IT_> col_idx_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::write));
      }

      Memory::TypedView<IT_> col_idx_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<IT_> col_idx_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(1).view(loc, acc));
      }

      /**
       * \brief Retrieve maximum bandwidth among all rows.
       *
       * \param[out] bandw The maximum bandwidth.
       * \param[out] bandw_i The row, where the bandwidth is maximal.
       */
      void bandwidth_row(Index & bandw, Index & bandw_i) const
      {
        XASSERT(!this->empty());

        const Memory::TypedView<IT_> row_ptr = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx = this->col_idx_view_r();
        bandw = 0;
        bandw_i = 0;
        for (Index row(0) ; row < num_rows() ; ++row)
        {
          if (row_ptr[row+1] == row_ptr[row])
            continue;

          Index temp = col_idx[row_ptr[row+1]-1] - col_idx[row_ptr[row]] + 1;
          if(temp > bandw)
          {
            bandw = temp;
            bandw_i = row;
          }
        }
      }

      /**
       * \brief Retrieve maximum bandwidth among all columns.
       *
       * \param[out] bandw The maximum bandwidth.
       * \param[out] bandw_i The column, where the bandwidth is maximal.
       */
      //void bandwidth_column(Index & bandw, Index & bandw_i) const
      //{
      //  /// \todo implement this properly
      //  SparseMatrixCSR<DT_, IT_> temp;
      //  temp.transpose(*this);
      //  temp.bandwidth_row(bandw, bandw_i);
      //}

      /**
       * \brief Retrieve maximum radius among all rows.
       *
       * The radius is defined as the maximal distance to the main diagonal
       *
       * \param[out] radius The maximum radius.
       * \param[out] radius_i The row, where the radius is maximal.
       */
      void radius_row(Index & radius, Index & radius_i) const
      {
        XASSERT(!this->empty());

        const Memory::TypedView<IT_> row_ptr = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx = this->col_idx_view_r();
        radius = 0;
        radius_i = 0;

        for (Index row(0) ; row < num_rows() ; ++row)
        {
          if (row_ptr[row+1] == row_ptr[row])
            continue;

          if (row_ptr[row+1] > 0)
          {
            Index temp1 = col_idx[row_ptr[row+1]-1];
            if(Math::max(temp1,row) - Math::min(temp1, row) > radius)
            {
              radius = Math::max(temp1,row) - Math::min(temp1, row);
              radius_i = row;
            }
          }
          Index temp2 = col_idx[row_ptr[row]];
          if(Math::max(temp2,row) - Math::min(temp2, row) > radius)
          {
            radius = Math::max(temp2,row) - Math::min(temp2, row);
            radius_i = row;
          }
        }
      }

      /**
       * \brief Retrieve maximum radius among all column.
       *
       * The radius is defined as the maximal distance to the main diagonal
       *
       * \param[out] radius The maximum radius.
       * \param[out] radius_i The column, where the radius is maximal.
       */
      //void radius_column(Index & radius, Index & radius_i) const
      //{
      //  /// \todo implement this properly
      //  SparseMatrixCSR<DT_, IT_> temp(this->clone());
      //  temp.transpose(temp);
      //  temp.radius_row(radius, radius_i);
      //}

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SparseMatrixCSR";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const SparseMatrixCSR & x, bool full = false)
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
       *
       * \warning All three matrices must have the same non zero layout. This operation assumes this silently and does not check this on its own!
       */
      void axpy(const SparseMatrixCSR & x, const DT_ alpha = DT_(1))
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
      void scale(const SparseMatrixCSR & x, const DT_ alpha)
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
       * \brief Computes the 2-norm for every row
       *
       * \param[in] row_norms
       * For every row, this left-vector will contain its 2-norm
       */
      void row_norm2(VectorTypeL& row_norms) const
      {
        ASSERTM(row_norms.size() == this->num_rows(), "Matrix/Vector dimension mismatch");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::RowNorm2CSR::template exec<DT_, IT_>(row_norms.elements_arbiter(), this->val_arbiter(),
          this->col_idx_arbiter(), this->row_ptr_arbiter(), this->num_rows(), false);

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes()*2);
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Computes the square of the 2-norm for every row
       *
       * \param[out] row_norms
       * For every row, this left-vector will contain the square of its 2-norm
       */
      void row_norm2sqr(VectorTypeL& row_norms) const
      {
        ASSERTM(row_norms.size() == this->num_rows(), "Matrix/Vector dimension mismatch");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::RowNorm2CSR::template exec<DT_, IT_>(row_norms.elements_arbiter(), this->val_arbiter(),
          this->col_idx_arbiter(), this->row_ptr_arbiter(), this->num_rows(), true);

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Retrieve the absolute maximum value of this matrix.
       *
       * \return The largest absolute value.
       */
      DT_ max_abs_element() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->val_arbiter(), this->num_nzes(), false, true);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the absolute minimum value of this matrix.
       *
       * \return The smallest absolute value
       */
      DT_ min_abs_element() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->val_arbiter(), this->num_nzes(), true, true);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the maximum value of this matrix.
       *
       * \return The largest value.
       */
      DT_ max_element() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->val_arbiter(), this->num_nzes(), false, false);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the minimum value of this matrix.
       *
       * \return The smallest value.
       */
      DT_ min_element() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->val_arbiter(), this->num_nzes(), true, false);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the maximum relative difference of this matrix and another one
       * y.max_rel_diff(x) returns  \f$ \max_{0\leq i < n}\frac{|x_i-y_i|}{\max{|x_i|+|y_i|, eps}} \f$
       *
       * \attention
       * This function silently assumes that both matrices have the same sparsity pattern!
       *
       * \return The largest relative difference.
       */
      DT_ max_rel_diff(const SparseMatrixCSR & x) const
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
      bool same_layout(const SparseMatrixCSR& x) const
      {
        if (this->num_rows() != x.num_rows())
          return false;
        if (this->num_cols() != x.num_cols())
          return false;
        if (this->num_nzes() != x.num_nzes())
          return false;

        if(this->num_nzes() == Index(0))
          return true;

        // check if the arbiters for row_ptr and col_idx are the same
        if((this->row_ptr_arbiter() == x.row_ptr_arbiter()) && (this->col_idx_arbiter() == x.col_idx_arbiter()))
          return true;

        const Memory::TypedView<IT_> row_ptr_a = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_a = this->col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_b = x.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_b = x.col_idx_view_r();

        const Index n = this->num_rows();
        for (Index i(0) ; i <= n; ++i)
        {
          if (row_ptr_a[i] != row_ptr_b[i])
          {
            return false;
          }
        }

        const Index nze = this->num_nzes();
        for (Index i(0) ; i < nze ; ++i)
        {
          if (col_idx_a[i] != col_idx_b[i])
          {
            return false;
          }
        }

        return true;
      }

      /**
       * \brief Calculate \f$this^\top \f$
       *
       * \return The transposed matrix
       */
      SparseMatrixCSR transpose() const
      {
        const Index txrows(this->num_rows());
        const Index txcolumns(this->num_cols());
        const Index txnonzeros(this->num_nzes());

        if(txnonzeros <= Index(0))
          return SparseMatrixCSR(txcolumns, txrows);

        SparseMatrixCSR<DT_, IT_> trans_mat(txcolumns, txrows, txnonzeros);

        {
          const Memory::TypedView<IT_> ptxcol_idx(this->col_idx_view_r());
          const Memory::TypedView<IT_> ptxrow_ptr(this->row_ptr_view_r());
          const Memory::TypedView<DT_> ptxval(this->val_view_r());

          Memory::TypedView<IT_> ptrow_ptr(trans_mat.row_ptr_view_w());
          Memory::TypedView<IT_> ptcol_idx(trans_mat.col_idx_view_w());
          Memory::TypedView<DT_> ptval(trans_mat.val_view_w());

          for (Index i(0); i <= txcolumns; ++i)
          {
            ptrow_ptr[i] = IT_(0);
          }

          for (Index i(0); i < txnonzeros; ++i)
          {
            ++ptrow_ptr[ptxcol_idx[i] + 1];
          }

          for (Index i(1); i < txcolumns; ++i)
          {
            ptrow_ptr[i + 1] += ptrow_ptr[i];
          }

          for (Index i(0); i < txrows; ++i)
          {
            for (IT_ k(ptxrow_ptr[i]); k < ptxrow_ptr[i+1]; ++k)
            {
              const IT_ l(ptxcol_idx[k]);
              const IT_ j(ptrow_ptr[l]);
              ptval[j] = ptxval[k];
              ptcol_idx[j] = IT_(i);
              ++ptrow_ptr[l];
            }
          }

          for (Index i(txcolumns); i > 0; --i)
          {
            ptrow_ptr[i] = ptrow_ptr[i - 1];
          }
          ptrow_ptr[0] = 0;
        }

        return trans_mat;
      }

      /**
       * \brief Calculate \f$this \leftarrow x^\top \f$
       *
       * \param[in] x The matrix to be transposed.
       */
      void transpose(const SparseMatrixCSR & x)
      {
        this->move(x.transpose());
      }

      /**
       * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_i\f$
       *
       * \param[in] x The matrix whose rows are to be scaled.
       * \param[in] s The vector to the scale the rows by.
       */
      void scale_rows(const SparseMatrixCSR & x, const DenseVector<DT_,IT_> & s)
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Row count does not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Column count does not match!");
        XASSERTM(x.num_nzes() == this->num_nzes(), "Nonzero count does not match!");
        XASSERTM(s.size() == this->num_rows(), "Vector size does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::ScaleRowCSR::template exec<DT_, IT_>(this->val_arbiter(), x.val_arbiter(), this->col_idx_arbiter(),
          this->row_ptr_arbiter(), s.elements_arbiter(), this->num_rows(), this->num_cols(), this->num_nzes());

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes());
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_j\f$
       *
       * \param[in] x The matrix whose columns are to be scaled.
       * \param[in] s The vector to the scale the columns by.
       */
      void scale_cols(const SparseMatrixCSR & x, const DenseVector<DT_,IT_> & s)
      {
        XASSERTM(x.num_rows() == this->num_rows(), "Row count does not match!");
        XASSERTM(x.num_cols() == this->num_cols(), "Column count does not match!");
        XASSERTM(x.num_nzes() == this->num_nzes(), "Nonzero count does not match!");
        XASSERTM(s.size() == this->num_cols(), "Vector size does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::ScaleColCSR::template exec<DT_, IT_>(this->val_arbiter(), x.val_arbiter(), this->col_idx_arbiter(),
          this->row_ptr_arbiter(), s.elements_arbiter(), this->num_rows(), this->num_cols(), this->num_nzes());

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes());
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \attention r and x must \b not refer to the same vector object!
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

        Arch::MatVecMultCSRDense::template exec<DT_, IT_>(r.elements_arbiter(), DT_(1), x.elements_arbiter(), Memory::Arbiter(),
          this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nzes(), false);

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2);
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
      * \brief Calculate \f$ r \leftarrow this^\top \cdot x \f$
      *
      * \attention r and x must \b not refer to the same vector object!
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

        Arch::MatVecMultCSRDense::template exec<DT_, IT_>(r.elements_arbiter(), DT_(1), x.elements_arbiter(), Memory::Arbiter(),
          this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nzes(), true);

        TimeStamp ts_stop;
        Statistics::add_flops(this->num_nzes() * 2);
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \attention r and x must \b not refer to the same vector object!
       *
       * \param[out] r The block vector that receives the result.
       * \param[in] x The block vector to be multiplied by this matrix.
       *
       * \note Every element of each block in the vector x is multiplied with the corresponding single scalar entry in the matrix.
       */
      template<int block_size_>
      void apply(DenseVectorBlocked<DT_, IT_, block_size_> & r, const DenseVectorBlocked<DT_, IT_, block_size_> & x) const
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

        Arch::MatVecMultCSRBlock::template exec<DT_, IT_, block_size_>(r.elements_arbiter(), DT_(1), x.elements_arbiter(),
          Memory::Arbiter(), this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nzes(), false);

        TimeStamp ts_stop;

        Statistics::add_flops(this->num_nzes() * 2 * block_size_);
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \attention r and x must \b not refer to the same vector object!
       * \note r and y are allowed to refer to the same vector object.
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
        else if((this->hollow()) || (alpha == DT_(0)))
        {
          r.copy(y);
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultCSRDense::template exec<DT_, IT_>(r.elements_arbiter(), alpha, x.elements_arbiter(), y.elements_arbiter(),
          this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nzes(), false);

        TimeStamp ts_stop;

        Statistics::add_flops( 2 * (this->num_nzes() + this->num_rows()) );
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
      * \brief Calculate \f$ r \leftarrow y + \alpha~ this^\top \cdot x \f$
      *
      * \attention r and x must \b not refer to the same vector object!
      * \note r and y are allowed to refer to the same vector object.
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
        else if((this->hollow()) || (alpha == DT_(0)))
        {
          r.copy(y);
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultCSRDense::template exec<DT_, IT_>(r.elements_arbiter(), alpha, x.elements_arbiter(), y.elements_arbiter(),
          this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nzes(), true);

        TimeStamp ts_stop;
        Statistics::add_flops( 2 * (this->num_nzes() + this->num_rows()) );
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \attention r and x must \b not refer to the same vector object!
       * \note r and y are allowed to refer to the same vector object.
       *
       * \param[out] r The block vector that receives the result.
       * \param[in] x The block vector to be multiplied by this matrix.
       * \param[in] y The summand block vector.
       * \param[in] alpha A scalar to scale the product with.
       *
       * \note Every element of each block in the vector x is multiplied with the corresponding single scalar entry in the matrix.
       */
      template<int block_size_>
      void apply(
        DenseVectorBlocked<DT_, IT_, block_size_> & r,
        const DenseVectorBlocked<DT_, IT_, block_size_> & x,
        const DenseVectorBlocked<DT_, IT_, block_size_> & y,
        const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->num_rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->num_cols(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->num_rows(), "Vector size of y does not match!");

        if(this->empty())
          return;
        // we check alpha strictly for equality to 0, because testing if |alpha| < eps may
        // lead to false triggers if the matrix/vector contents are also < eps
        else if((this->hollow()) || (alpha == DT_(0)))
        {
          r.copy(y);
          return;
        }

        XASSERTM(r.elements_arbiter() != x.elements_arbiter(), "Vector x and r must not share the same memory!");

        TimeStamp ts_start;

        Arch::MatVecMultCSRBlock::template exec<DT_, IT_, block_size_>(r.elements_arbiter(),alpha, x.elements_arbiter(),
          y.elements_arbiter(), this->val_arbiter(), this->col_idx_arbiter(), this->row_ptr_arbiter(),
          this->num_rows(), this->num_cols(), this->num_nzes(), false);

        TimeStamp ts_stop;

        Statistics::add_flops( (this->num_nzes() + this->num_rows()) * 2 * block_size_);
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Adds a double-matrix product onto this matrix
       *
       * This function performs the following computation:
       * \f[ X \leftarrow X + \alpha D\cdot A\cdot B\f]
       *
       * where
       * - \e X denotes this m-by-n matrix
       * - \e D denotes a m-by-k matrix
       * - \e A denotes a k-by-l matrix
       * - \e B denotes a l-by-n matrix
       *
       * \attention
       * This function assumes that the output matrix already contains the
       * required sparsity pattern. This function will throw an exception
       * if the sparsity pattern of the output matrix is incomplete unless
       * \p allow_incomplete is set to \c true.
       *
       * \note
       * This function currently only supports data in main memory.
       *
       * \param[in] d, a, b
       * The three matrices to be multiplied
       *
       * \param[in] alpha
       * The scaling factor for the product
       *
       * \param[in] allow_incomplete
       * Specifies whether the output matrix structure is allowed to be incomplete.
       * If set to \c false, this function will throw an exception on incompleteness,
       * otherwise the missing entries are ignored (dropped).
       */
      void add_double_mat_product(
        const LAFEM::SparseMatrixCSR<DT_, IT_>& d,
        const LAFEM::SparseMatrixCSR<DT_, IT_>& a,
        const LAFEM::SparseMatrixCSR<DT_, IT_>& b,
        const DT_ alpha = DT_(1),
        const bool allow_incomplete = false)
      {
        // validate matrix dimensions
        XASSERT(this->num_rows() == d.num_rows());
        XASSERT(d.num_cols() == a.num_rows());
        XASSERT(a.num_cols() == b.num_rows());
        XASSERT(b.num_cols() == this->num_cols());

        // fetch matrix arrays:
        Memory::TypedView<DT_> data_x = this->val_view_rw();
        const Memory::TypedView<DT_> data_d = d.val_view_r();
        const Memory::TypedView<DT_> data_a = a.val_view_r();
        const Memory::TypedView<DT_> data_b = b.val_view_r();
        const Memory::TypedView<IT_> row_ptr_x = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_x = this->col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_d = d.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_d = d.col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_a = a.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_a = a.col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_b = b.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_b = b.col_idx_view_r();

        // loop over all rows of D and X, resp.
        for(IT_ i(0); i < IT_(this->num_rows()); ++i)
        {
          // loop over all non-zeros D_ik in row i of D
          for(IT_ ik(row_ptr_d[i]); ik  < row_ptr_d[i+1]; ++ik)
          {
            // get column index k
            const IT_ k = col_idx_d[ik];

            // loop over all non-zeros A_kl in row k of A
            for(IT_ kl(row_ptr_a[k]); kl < row_ptr_a[k+1]; ++kl)
            {
              // get column index l
              const IT_ l = col_idx_a[kl];

              // pre-compute factor (alpha * D_ik * A_kl)
              const DT_ omega = alpha *  data_d[ik] * data_a[kl];

              // loop over all non-zeros B_lj in row j of B and
              // loop over all non-zeros X_ij in row i of X and
              // perform a "sparse axpy" of B_l onto X_i, i.e.:
              //   X_i. += (alpha * D_ik * A_kl) * B_l.
              IT_ ij = row_ptr_x[i];
              IT_ lj = row_ptr_b[l];
              while(lj < row_ptr_b[l+1])
              {
                if(ij >= row_ptr_x[i+1])
                {
                  // we have reached the end of row X_i, but there is at least
                  // one entry in row B_l left, so the pattern of X is incomplete
                  // We let the caller decide whether this is a valid case or not:
                  if(allow_incomplete)
                    break;  // continue with next row
                  else
                    XABORTM("Incomplete output matrix structure");
                }
                else if(col_idx_x[ij] == col_idx_b[lj])
                {
                  // okay: B_lj contributes to X_ij
                  data_x[ij] += omega * data_b[lj];
                  ++ij;
                  ++lj;
                }
                else if(col_idx_x[ij] < col_idx_b[lj])
                {
                  // entry X_ij exists, but B_lj is missing:
                  // this is a perfectly valid case, so continue with the next non-zero of X_i
                  ++ij;
                }
                else //if(col_idx_x[ij] > col_idx_b[lj])
                {
                  // If we come out here, then the sparsity pattern of X is incomplete:
                  // B_lj is meant to be added onto X_ij, but the entry X_ij is missing
                  // We let the caller decide whether this is a valid case or not:
                  if(allow_incomplete)
                    ++lj;
                  else
                    XABORTM("Incomplete output matrix structure");
                }
              }
            }
          }
        }
      }

      /**
       * \brief Adds a double-matrix product onto this matrix
       *
       * This function performs the following computation:
       * \f[ X \leftarrow X + \alpha D\cdot \textnormal{diag}(A)\cdot B\f]
       *
       * where
       * - \e X denotes this m-by-n matrix
       * - \e D denotes a m-by-l matrix
       * - \e A denotes a vector representing a l-by-l diagonal matrix
       * - \e B denotes a l-by-n matrix
       *
       * \attention
       * This function assumes that the output matrix already contains the
       * required sparsity pattern. This function will throw an exception
       * if the sparsity pattern of the output matrix is incomplete unless
       * \p allow_incomplete is set to \c true.
       *
       * \note
       * This function currently only supports data in main memory.
       *
       * \param[in] a
       * The vector representing the diagonal matrix A.
       *
       * \param[in] d, b
       * The left and right multiplicand matrices
       *
       * \param[in] alpha
       * The scaling factor for the product
       *
       * \param[in] allow_incomplete
       * Specifies whether the output matrix structure is allowed to be incomplete.
       * If set to \c false, this function will throw an exception on incompleteness,
       * otherwise the missing entries are ignored (dropped).
       */
      void add_double_mat_product(
        const LAFEM::SparseMatrixCSR<DT_, IT_>& d,
        const LAFEM::DenseVector<DT_, IT_>& a,
        const LAFEM::SparseMatrixCSR<DT_, IT_>& b,
        const DT_ alpha = DT_(1),
        const bool allow_incomplete = false)
      {
        // validate matrix dimensions
        XASSERT(this->num_rows() == d.num_rows());
        XASSERT(d.num_cols() == a.size());
        XASSERT(a.size() == b.num_rows());
        XASSERT(b.num_cols() == this->num_cols());

        // fetch matrix arrays:
        Memory::TypedView<DT_> data_x = this->val_view_rw();
        const Memory::TypedView<DT_> data_d = d.val_view_r();
        const Memory::TypedView<DT_> data_a = a.elements_view_r();
        const Memory::TypedView<DT_> data_b = b.val_view_r();
        const Memory::TypedView<IT_> row_ptr_x = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_x = this->col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_d = d.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_d = d.col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_b = b.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_b = b.col_idx_view_r();

        // loop over all rows of D and X, resp.
        for(IT_ i(0); i < IT_(this->num_rows()); ++i)
        {
          // loop over all non-zeros D_ik in row i of D
          for(IT_ ik(row_ptr_d[i]); ik  < row_ptr_d[i+1]; ++ik)
          {
            // get column index k
            const IT_ k = col_idx_d[ik];

            // pre-compute factor (alpha * D_ik * A_kk)
            const DT_ omega = alpha *  data_d[ik] * data_a[k];

            // loop over all non-zeros B_kj in row j of B and
            // loop over all non-zeros X_ij in row i of X and
            // perform a "sparse axpy" of B_l onto X_i, i.e.:
            //   X_i. += (alpha * D_ik * A_kk) * B_k.
            IT_ ij = row_ptr_x[i];
            IT_ kj = row_ptr_b[k];
            while(kj < row_ptr_b[k+1])
            {
              if(ij >= row_ptr_x[i+1])
              {
                // we have reached the end of row X_i, but there is at least
                // one entry in row B_l left, so the pattern of X is incomplete
                // We let the caller decide whether this is a valid case or not:
                if(allow_incomplete)
                  break; // continue with next row
                else
                  XABORTM("Incomplete output matrix structure");
              }
              else if(col_idx_x[ij] == col_idx_b[kj])
              {
                // okay: B_kj contributes to X_ij
                data_x[ij] += omega * data_b[kj];
                ++ij;
                ++kj;
              }
              else if(col_idx_x[ij] < col_idx_b[kj])
              {
                // entry X_ij exists, but B_kj is missing:
                // this is a perfectly valid case, so continue with the next non-zero of X_i
                ++ij;
              }
              else //if(col_idx_x[ij] > col_idx_b[kj])
              {
                // If we come out here, then the sparsity pattern of X is incomplete:
                // B_kj is meant to be added onto X_ij, but the entry X_ij is missing
                // We let the caller decide whether this is a valid case or not:
                if(allow_incomplete)
                  ++kj;
                else
                  XABORTM("Incomplete output matrix structure");
              }
            }
          }
        }
      }

      /**
       * \brief Adds a double-matrix product onto this matrix
       *
       * This function performs the following computation:
       * \f[ X \leftarrow X + \alpha D\cdot \textnormal{diag}(A)\cdot B\f]
       *
       * where
       * - \e X denotes this m-by-n matrix
       * - \e D denotes a m-by-l matrix
       * - \e A denotes a vector representing a l-by-l diagonal matrix
       * - \e B denotes a l-by-n matrix
       *
       * \attention
       * This function assumes that the output matrix already contains the
       * required sparsity pattern. This function will throw an exception
       * if the sparsity pattern of the output matrix is incomplete unless
       * \p allow_incomplete is set to \c true.
       *
       * \note
       * This function currently only supports data in main memory.
       *
       * \param[in] a
       * The vector representing the diagonal matrix A.
       *
       * \param[in] d, b
       * The left and right multiplicand matrices
       *
       * \param[in] alpha
       * The scaling factor for the product
       *
       * \param[in] allow_incomplete
       * Specifies whether the output matrix structure is allowed to be incomplete.
       * If set to \c false, this function will throw an exception on incompleteness,
       * otherwise the missing entries are ignored (dropped).
       */
      template<int bs_>
      void add_double_mat_product(
        const LAFEM::SparseMatrixBCSR<DT_, IT_, 1, bs_>& d,
        const LAFEM::DenseVectorBlocked<DT_, IT_, bs_>& a,
        const LAFEM::SparseMatrixBCSR<DT_, IT_, bs_, 1>& b,
        const DT_ alpha = DT_(1),
        const bool allow_incomplete = false)
      {
        // validate matrix dimensions
        XASSERT(this->num_rows() == d.num_rows());
        XASSERT(d.num_cols() == a.size());
        XASSERT(a.size() == b.num_rows());
        XASSERT(b.num_cols() == this->num_cols());

        // fetch matrix arrays:
        Memory::TypedView<DT_> data_x = this->val_view_rw();
        const Memory::TypedView<Tiny::Matrix<DT_, 1, bs_>> data_d = d.val_view_r();
        const Memory::TypedView<Tiny::Vector<DT_, bs_>> data_a = a.elements_view_r();
        const Memory::TypedView<Tiny::Matrix<DT_, bs_, 1>> data_b = b.val_view_r();
        const Memory::TypedView<IT_> row_ptr_x = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_x = this->col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_d = d.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_d = d.col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_b = b.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_b = b.col_idx_view_r();

        // loop over all rows of D and X, resp.
        for(IT_ i(0); i < IT_(this->num_rows()); ++i)
        {
          // loop over all non-zeros D_ik in row i of D
          for(IT_ ik(row_ptr_d[i]); ik  < row_ptr_d[i+1]; ++ik)
          {
            // get column index k
            const IT_ k = col_idx_d[ik];

            // pre-compute factor (alpha * D_ik * A_kk)
            Tiny::Vector<DT_, bs_> omega;
            for(int l = 0; l < bs_; ++l)
              omega[l] = alpha * data_d[ik](0,l) * data_a[k][l];

            // loop over all non-zeros B_kj in row j of B and
            // loop over all non-zeros X_ij in row i of X and
            // perform a "sparse axpy" of B_l onto X_i, i.e.:
            //   X_i. += (alpha * D_ik * A_kk) * B_k.
            IT_ ij = row_ptr_x[i];
            IT_ kj = row_ptr_b[k];
            while(kj < row_ptr_b[k+1])
            {
              if(ij >= row_ptr_x[i+1])
              {
                // we have reached the end of row X_i, but there is at least
                // one entry in row B_l left, so the pattern of X is incomplete
                // We let the caller decide whether this is a valid case or not:
                if(allow_incomplete)
                  break; // continue with next row
                else
                  XABORTM("Incomplete output matrix structure");
              }
              else if(col_idx_x[ij] == col_idx_b[kj])
              {
                // okay: B_kj contributes to X_ij
                for(int l = 0; l < bs_; ++l)
                  data_x[ij] += omega[l] * data_b[kj](l,0);
                ++ij;
                ++kj;
              }
              else if(col_idx_x[ij] < col_idx_b[kj])
              {
                // entry X_ij exists, but B_kj is missing:
                // this is a perfectly valid case, so continue with the next non-zero of X_i
                ++ij;
              }
              else //if(col_idx_x[ij] > col_idx_b[kj])
              {
                // If we come out here, then the sparsity pattern of X is incomplete:
                // B_kj is meant to be added onto X_ij, but the entry X_ij is missing
                // We let the caller decide whether this is a valid case or not:
                if(allow_incomplete)
                  ++kj;
                else
                  XABORTM("Incomplete output matrix structure");
              }
            }
          }
        }
      }

      /**
       * \brief Adds a matrix-matrix product onto this matrix
       *
       * This function performs the following computation:
       * \f[ X \leftarrow X + \alpha D\cdot B\f]
       *
       * where
       * - \e X denotes this m-by-n matrix
       * - \e D denotes a m-by-l matrix
       * - \e B denotes a l-by-n matrix
       *
       * \attention
       * This function assumes that the output matrix already contains the
       * required sparsity pattern. This function will throw an exception
       * if the sparsity pattern of the output matrix is incomplete unless
       * \p allow_incomplete is set to \c true.
       *
       * \note
       * This function currently only supports data in main memory.
       *
       * \param[in] d, b
       * The left and right multiplicand matrices
       *
       * \param[in] alpha
       * The scaling factor for the product
       *
       * \param[in] allow_incomplete
       * Specifies whether the output matrix structure is allowed to be incomplete.
       * If set to \c false, this function will throw an exception on incompleteness,
       * otherwise the missing entries are ignored (dropped).
       */
      void add_mat_mat_product(
        const LAFEM::SparseMatrixCSR<DT_, IT_>& d,
        const LAFEM::SparseMatrixCSR<DT_, IT_>& b,
        const DT_ alpha = DT_(1),
        const bool allow_incomplete = false)
      {
        // validate matrix dimensions
        XASSERT(this->num_rows() == d.num_rows());
        XASSERT(d.num_cols() == b.num_rows());
        XASSERT(b.num_cols() == this->num_cols());

        // fetch matrix arrays:
        Memory::TypedView<DT_> data_x = this->val_view_rw();
        const Memory::TypedView<DT_> data_d = d.val_view_r();
        const Memory::TypedView<DT_> data_b = b.val_view_r();
        const Memory::TypedView<IT_> row_ptr_x = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_x = this->col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_d = d.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_d = d.col_idx_view_r();
        const Memory::TypedView<IT_> row_ptr_b = b.row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx_b = b.col_idx_view_r();

        // loop over all rows of D and X, resp.
        for(IT_ i(0); i < IT_(this->num_rows()); ++i)
        {
          // loop over all non-zeros D_ik in row i of D
          for(IT_ ik(row_ptr_d[i]); ik  < row_ptr_d[i+1]; ++ik)
          {
            // get column index k
            const IT_ k = col_idx_d[ik];

            // pre-compute factor (alpha * D_ik)
            const DT_ omega = alpha *  data_d[ik];

            // loop over all non-zeros B_kj in row j of B and
            // loop over all non-zeros X_ij in row i of X and
            // perform a "sparse axpy" of B_l onto X_i, i.e.:
            //   X_i. += (alpha * D_ik) * B_k.
            IT_ ij = row_ptr_x[i];
            IT_ kj = row_ptr_b[k];
            while(kj < row_ptr_b[k+1])
            {
              if(ij >= row_ptr_x[i+1])
              {
                // we have reached the end of row X_i, but there is at least
                // one entry in row B_l left, so the pattern of X is incomplete
                // We let the caller decide whether this is a valid case or not:
                if(allow_incomplete)
                  break; // continue with next row
                else
                  XABORTM("Incomplete output matrix structure");
              }
              else if(col_idx_x[ij] == col_idx_b[kj])
              {
                // okay: B_kj contributes to X_ij
                data_x[ij] += omega * data_b[kj];
                ++ij;
                ++kj;
              }
              else if(col_idx_x[ij] < col_idx_b[kj])
              {
                // entry X_ij exists, but B_kj is missing:
                // this is a perfectly valid case, so continue with the next non-zero of X_i
                ++ij;
              }
              else //if(col_idx_x[ij] > col_idx_b[kj])
              {
                // If we come out here, then the sparsity pattern of X is incomplete:
                // B_kj is meant to be added onto X_ij, but the entry X_ij is missing
                // We let the caller decide whether this is a valid case or not:
                if(allow_incomplete)
                  ++kj;
                else
                  XABORTM("Incomplete output matrix structure");
              }
            }
          }
        }
      }
      ///@}

      /// \copydoc lump_rows()
      void lump_rows(VectorTypeL& lump) const
      {
        XASSERTM(lump.size() == num_rows(), "lump vector size does not match matrix row count!");
        if(this->empty())
          return;

        Arch::LumpingCSR::template exec<DT_, IT_>(lump.elements_arbiter(), this->val_arbiter(),
          this->col_idx_arbiter(), this->row_ptr_arbiter(), this->num_rows());
      }

      /**
       * \brief Returns the lumped rows vector
       *
       * Each entry in the returned lumped rows vector contains the
       * the sum of all matrix elements in the corresponding row.
       *
       * \returns
       * The lumped vector.
       */
      VectorTypeL lump_rows() const
      {
        VectorTypeL lump = create_vector_l();
        lump_rows(lump);
        return lump;
      }

      /**
       * \brief extract main diagonal vector from matrix
       *
       * \param[out] diag
       * The vector containing the diagonal entry values
       */
      void extract_diag(VectorTypeL & diag) const
      {
        if(this->empty())
          return;

        Arch::DiagonalCSR::template exec<DT_, IT_>(diag.elements_arbiter(), this->val_arbiter(),
          this->col_idx_arbiter(), this->row_ptr_arbiter(), this->num_rows());
      }

      /**
       * \brief extract main diagonal vector from matrix
       *
       * \returns The vector containing the diagonal entry values
       */
      VectorTypeL extract_diag() const
      {
        VectorTypeL diag = create_vector_l();
        extract_diag(diag);
        return diag;
      }

      /// Permute matrix rows and columns according to the given Permutations
      void permute(const Adjacency::Permutation & perm_row, const Adjacency::Permutation & perm_col)
      {
        if (perm_row.size() == 0 && perm_col.size() == 0)
          return;

        XASSERTM(perm_row.size() == this->num_rows(), "Container rows does not match permutation size");
        XASSERTM(perm_col.size() == this->num_cols(), "Container columns does not match permutation size");

        // http://de.mathworks.com/help/matlab/math/sparse-matrix-operations.html#f6-13070
        IT_ * temp_row_ptr = new IT_[num_rows() + 1];
        IT_ * temp_col_idx = new IT_[num_nzes()];
        DT_ * temp_val = new DT_[num_nzes()];

        Memory::TypedView<IT_> row_ptr = this->row_ptr_view_rw();
        Memory::TypedView<IT_> col_idx = this->col_idx_view_rw();
        Memory::TypedView<DT_> val = this->val_view_rw();

        const Index * perm_pos = perm_row.get_perm_pos();

        //permute rows from source to temp_*
        Index new_start(0);
        temp_row_ptr[0] = 0;
        for (Index row(0) ; row < this->num_rows() ; ++row)
        {
          Index row_size(row_ptr[perm_pos[row] + 1] - row_ptr[perm_pos[row]]);

          //iterate over all elements in single one new and old row
          for (Index i(new_start), j(row_ptr[perm_pos[row]]) ; i < new_start + row_size ; ++i, ++j)
          {
            temp_col_idx[i] = col_idx[j];
            temp_val[i] = val[j];
          }

          new_start += row_size;
          temp_row_ptr[row+1] = (IT_)new_start;
        }

        //use inverse col permutation as lookup table: i -> new location of i
        Adjacency::Permutation perm_col_inv = perm_col.inverse();
        perm_pos = perm_col_inv.get_perm_pos();

        //permute columns from temp_* to result
        Memory::memcopy_main(row_ptr.get_w(), temp_row_ptr, (num_rows() + 1) * sizeof(IT_));
        Memory::memcopy_main(val.get_w(), temp_val, num_nzes() * sizeof(DT_));
        for (Index i(0) ; i < num_nzes() ; ++i)
        {
          col_idx[i] = (IT_)perm_pos[temp_col_idx[i]];
        }

        delete[] temp_row_ptr;
        delete[] temp_col_idx;
        delete[] temp_val;

        //sort columns in every row by column index
        IT_ swap_key;
        DT_ swap_val;
        for (Index row(0) ; row < num_rows() ; ++row)
        {
          Index offset(row_ptr[row]);
          Index row_size(row_ptr[row+1] - row_ptr[row]);
          for (Index i(1), j ; i < row_size ; ++i)
          {
            swap_key = col_idx[i + offset];
            swap_val = val[i + offset];
            j = i;
            while (j > 0 && col_idx[j - 1 + offset] > swap_key)
            {
              col_idx[j + offset] = col_idx[j - 1 + offset];
              val[j + offset] = val[j - 1 + offset];
              --j;
            }
            col_idx[j + offset] = swap_key;
            val[j + offset] = swap_val;
          }
        }
      }

      /// Shrink matrix and drop small values
      /**
       *
       * Shrinks the matrix by dropping all values, that have a smaller absolute value than eps.
       *
       * \param[in] eps The dropping criterion
       *
       **/
      void shrink(DT_ eps)
      {
        const Memory::TypedView<IT_> row_ptr = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx = this->col_idx_view_r();
        const Memory::TypedView<DT_> val = this->val_view_r();

        std::vector<IT_> row_entries(this->num_rows(), IT_(0));
        for (Index row(0) ; row < this->num_rows() ; ++row)
        {
          for (Index el(row_ptr[row]) ; el < row_ptr[row + 1] ; ++el)
          {
            if (Math::abs(val[el]) >= eps)
            {
              row_entries.at(row) += IT_(1);
            }
          }
        }

        Index ue(0);
        for (auto& n : row_entries)
        {
          ue += n;
        }
        if (ue == Index(0))
        {
          this->move(SparseMatrixCSR<DT_, IT_>(this->num_rows(), this->num_cols()));
          return;
        }

        DenseVector<IT_, IT_> new_row_ptr(this->num_rows() + 1);
        {
          Memory::TypedView<IT_> nrp = new_row_ptr.elements_view_w();
          nrp[0] = IT_(0);
          for (Index row(0) ; row < this->num_rows() ; ++row)
          {
            //new_row_ptr(row + 1, new_row_ptr(row) + row_entries.at(row));
            nrp[row + 1] = nrp[row] + row_entries.at(row);
          }
        }

        row_entries.clear();

        DenseVector<IT_, IT_> new_col_idx(ue);
        DenseVector<DT_, IT_> new_val(ue);
        Index counter(0);
        {
          Memory::TypedView<IT_> nci = new_col_idx.elements_view_w();
          Memory::TypedView<DT_> nval = new_val.elements_view_w();
          for (Index row(0) ; row < this->num_rows() ; ++row)
          {
            for (Index el(row_ptr[row]) ; el < row_ptr[row + 1] ; ++el)
            {
              if (Math::abs(val[el]) >= eps)
              {
                nci[counter] = col_idx[el];
                nval[counter] = val[el];
                ++counter;
              }
            }
          }
        }

        // release views
        row_ptr.release();
        col_idx.release();
        val.release();

        this->move(SparseMatrixCSR<DT_, IT_>(this->num_rows(), this->num_cols(), new_row_ptr, new_col_idx, new_val));
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

      /// Returns the number of NNZ-elements of the selected row
      Index row_degree(const Index row) const
      {
        const Memory::TypedView<IT_> prow_ptr(this->row_ptr_view_r());
        return Index(prow_ptr[row + 1] - prow_ptr[row]);
      }

      /**
       * \brief Extracts the column indices for a given row
       *
       * \param[in] row
       * The index of the row whose column indices are to be extracted
       *
       * \param[out] pcol_idx
       * A \transient array that receives the column indices
       *
       * \param[in] col_offset
       * An offset that is to be added onto each column index
       *
       * \returns The number of column indices extracted
       */
      template<typename IT2_>
      Index get_row_col_indices(const Index row, IT2_ * const pcol_idx, const IT2_ col_offset) const
      {
        const Memory::TypedView<IT_> trow_ptr(this->row_ptr_view_r());
        const Memory::TypedView<IT_> tcol_idx(this->col_idx_view_r());

        const Index start = Index(trow_ptr[row]);
        const Index end = Index(trow_ptr[row + 1] - trow_ptr[row]);
        for (Index i(0); i < end; ++i)
        {
          pcol_idx[i] = IT2_(tcol_idx[start + i]) + col_offset;
        }
        return end;
      }

      /**
       * \brief Extracts the values for a given row
       *
       * \param[in] row
       * The index of the row whose values are to be extracted
       *
       * \param[out] pvals
       * A \transient array that receives the values
       *
       * \returns The number of values extracted
       */
      template<typename DT2_>
      Index get_row_values(const Index row, DT2_ * const pvals) const
      {
        const Memory::TypedView<IT_> prow_ptr(this->row_ptr_view_r());
        const Memory::TypedView<DT_> pval(this->val_view_r());

        const Index start = Index(prow_ptr[row]);
        const Index end = Index(prow_ptr[row + 1] - prow_ptr[row]);
        for (Index i(0); i < end; ++i)
        {
          pvals[i] = DT2_(pval[start + i]);
        }
        return end;
      }

      /**
       * \brief Overwrites the values for a given row
       *
       * \param[in] row
       * The index of the row whose values are to be written
       *
       * \param[out] pvals
       * A \transient array containing the values to write to the row
       *
       * \returns The number of values written
       */
      template<typename DT2_>
      Index set_row_values(const Index row, const DT2_ * const pvals)
      {
        const Memory::TypedView<IT_> prow_ptr(this->row_ptr_view_r());
        Memory::TypedView<DT_> pval(this->val_view_rw());

        const Index start = Index(prow_ptr[row]);
        const Index end = Index(prow_ptr[row + 1] - prow_ptr[row]);
        for (Index i(0); i < end; ++i)
        {
          pval[start + i] = DT_(pvals[i]);
        }
        return end;
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
      void print(std::ostream & os, bool print_dense) const
      {
        if(this->empty())
        {
          os << "[]";
          return;
        }

        const Memory::TypedView<IT_> row_ptr = this->row_ptr_view_r();
        const Memory::TypedView<IT_> col_idx = this->col_idx_view_r();
        const Memory::TypedView<DT_> val = this->val_view_r();
        const Index nrows = this->num_rows();
        const Index ncols = this->num_cols();

        os << "[\n";
        for (Index i(0) ; i < nrows ; ++i)
        {
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
          }
          os << "]\n";
        }
        os << "]";
      }

      /**
       * \brief SparseMatrixCSR streaming operator
       *
       * \param[in] os The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & os, const SparseMatrixCSR & b)
      {
        b.print(os, true);
        return os;
      }

    public:
      /**
       * \brief Adjactor Wrapper for SparseMatrixCSR
       *
       * \author Peter Zajac
       */
      class Adjactor
      {
      public:
        /// ImageIterator typedef for Adjactor interface implementation
        typedef const IT_* ImageIterator;

      private:
        const SparseMatrixCSR& _matrix;
        const Memory::TypedView<IT_> _row_ptr;
        const Memory::TypedView<IT_> _col_idx;

      public:
        explicit Adjactor(const SparseMatrixCSR& matrix) :
          _matrix(matrix),
          _row_ptr(matrix.row_ptr_view_r()),
          _col_idx(matrix.col_idx_view_r())
        {
        }

        Adjactor(const Adjactor&) = delete;
        Adjactor& operator=(const Adjactor&) = delete;

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
          XASSERTM(domain_node < _matrix.num_rows(), "Domain node index out of range");
          return &_col_idx.get_r()[_row_ptr[domain_node]];
        }

        /** \copydoc Adjactor::image_end() */
        inline ImageIterator image_end(Index domain_node) const
        {
          XASSERTM(domain_node < _matrix.num_rows(), "Domain node index out of range");
          return &_col_idx.get_r()[_row_ptr[domain_node + 1u]];
        }
      }; // class Adjactor

      /**
       * \brief Scatter-Axpy operation for SparseMatrixCSR
       *
       * \author Peter Zajac
       */
      class ScatterAxpy
      {
      public:
        typedef LAFEM::SparseMatrixCSR<DT_, IT_> MatrixType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
#ifdef DEBUG
        static constexpr IT_ _deadcode = ~IT_(0);
#endif
        const Index _num_rows;
        const Index _num_cols;
        Memory::Arbiter _col_ptr_arbiter;
        const Memory::TypedView<IT_> _row_ptr_view;
        const Memory::TypedView<IT_> _col_idx_view;
        Memory::TypedView<IT_> _col_ptr_view;
        Memory::TypedView<DT_> _data_view;

      public:
        explicit ScatterAxpy(MatrixType& matrix) :
          _num_rows(matrix.num_rows()),
          _num_cols(matrix.num_cols()),
          _col_ptr_arbiter(_num_cols * sizeof(IT_), Memory::Location::main, Memory::Init::format_to_one),
          _row_ptr_view(matrix.row_ptr_view_r()),
          _col_idx_view(matrix.col_idx_view_r()),
          _col_ptr_view(_col_ptr_arbiter.view(Memory::Location::main, Memory::Access::read_write)),
          _data_view(matrix.val_view(Memory::Location::main, Memory::Access::read_write | Memory::Access::overlap))
        {
        }

        template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
        void operator()(const LocalMatrix_& loc_mat, const RowMapping_& row_map,
          const ColMapping_& col_map, DT_ alpha = DT_(1))
        {
          // note that initially all bits of _col_ptr_arbiter are formatted to 1,
          // which corresponds to the '_deadcode' value, so we don't need to reformat it here

          // loop over all local row entries
          for(int i(0); i < row_map.get_num_local_dofs(); ++i)
          {
            // fetch row index
            const Index ix = row_map.get_index(i);

            // build column pointer for this row entry contribution
            for(IT_ k(_row_ptr_view(ix)); k < _row_ptr_view(ix + 1); ++k)
            {
              _col_ptr_view[_col_idx_view(k)] = k;
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // fetch column index
              const Index jx = col_map.get_index(j);

              // ensure that the column pointer is valid for this index
              ASSERTM(_col_ptr_view[jx] != _deadcode, "invalid column index");

              // incorporate data into global matrix
              _data_view[_col_ptr_view[jx]] += alpha * loc_mat[i][j];

              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(IT_ k(_row_ptr_view(ix)); k < _row_ptr_view(ix + 1); ++k)
            {
              _col_ptr_view[_col_idx_view(k)] = _deadcode;
            }
#endif
            // continue with next row entry
          }
        }
      }; // class ScatterAxpy

      /**
       * \brief Gather-Axpy operation for SparseMatrixCSR
       *
       * \author Peter Zajac
       */
      class GatherAxpy
      {
      public:
        typedef LAFEM::SparseMatrixCSR<DT_, IT_> MatrixType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
#ifdef DEBUG
        static constexpr IT_ _deadcode = ~IT_(0);
#endif
        const Index _num_rows;
        const Index _num_cols;
        Memory::Arbiter _col_ptr_arbiter;
        const Memory::TypedView<IT_> _row_ptr_view;
        const Memory::TypedView<IT_> _col_idx_view;
        Memory::TypedView<IT_> _col_ptr_view;
        const Memory::TypedView<DT_> _data_view;

      public:
        explicit GatherAxpy(const MatrixType& matrix) :
          _num_rows(matrix.num_rows()),
          _num_cols(matrix.num_cols()),
          _col_ptr_arbiter(_num_cols * sizeof(IT_), Memory::Location::main, Memory::Init::format_to_one),
          _row_ptr_view(matrix.row_ptr_view_r()),
          _col_idx_view(matrix.col_idx_view_r()),
          _col_ptr_view(_col_ptr_arbiter.view(Memory::Location::main, Memory::Access::read_write)),
          _data_view(matrix.val_view_r())
        {
        }

        template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
        void operator()(LocalMatrix_& loc_mat, const RowMapping_& row_map, const ColMapping_& col_map, DT_ alpha = DT_(1))
        {
          // note that initially all bits of _col_ptr_arbiter are formatted to 1,
          // which corresponds to the '_deadcode' value, so we don't need to reformat it here

          // loop over all local row entries
          for(int i(0); i < row_map.get_num_local_dofs(); ++i)
          {
            // fetch row index
            const Index ix = row_map.get_index(i);

            // build column pointer for this row entry contribution
            for(IT_ k(_row_ptr_view(ix)); k < _row_ptr_view(ix + 1); ++k)
            {
              _col_ptr_view[_col_idx_view(k)] = k;
            }

            // loop over all local column entries
            for(int j(0); j < col_map.get_num_local_dofs(); ++j)
            {
              // fetch column index
              const Index jx = col_map.get_index(j);

              // ensure that the column pointer is valid for this index
              ASSERTM(_col_ptr_view[jx] != _deadcode, "invalid column index");

              // update local matrix data
              loc_mat[i][j] += alpha * _data_view[_col_ptr_view[jx]];

              // continue with next column entry
            }

#ifdef DEBUG
            // reformat column-pointer array
            for(IT_ k(_row_ptr_view(ix)); k < _row_ptr_view(ix + 1); ++k)
            {
              _col_ptr_view[_col_idx_view(k)] = _deadcode;
            }
#endif

            // continue with next row entry
          }
        }
      }; // class GatherAxpy

      Adjactor adjactor() const
      {
        return Adjactor(*this);
      }
    }; // class SparseMatrixCSR<...>

#ifdef FEAT_EICKT
    extern template class SparseMatrixCSR<float, std::uint32_t>;
    extern template class SparseMatrixCSR<double, std::uint32_t>;
    extern template class SparseMatrixCSR<float, std::uint64_t>;
    extern template class SparseMatrixCSR<double, std::uint64_t>;
#endif
  } // namespace LAFEM
} // namespace FEAT
