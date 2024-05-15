// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_BCSR_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_BCSR_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/apply.hpp>
#include <kernel/lafem/arch/lumping.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/lafem/arch/scale_row_col.hpp>
#include <kernel/lafem/arch/row_norm.hpp>
#include <kernel/lafem/arch/diagonal.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/likwid_marker.hpp>

#include <fstream>

namespace FEAT
{
  namespace LAFEM
  {
    /// \cond internal
    namespace Intern
    {
      template<typename DT_, typename IT_, int size_>
      struct BCSRVectorHelper
      {
        static_assert(size_ > 1, "invalid block size");
        typedef DenseVectorBlocked<DT_, IT_, size_> VectorType;
      };

      template<typename DT_, typename IT_>
      struct BCSRVectorHelper<DT_, IT_, 1>
      {
        typedef DenseVector< DT_, IT_> VectorType;
      };

      template<typename DT_, int BlockHeight_, int BlockWidth_, Perspective perspective_>
      struct BCSRPerspectiveHelper
      {
        typedef Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> Type;
      };

      template<typename DT_, int BlockHeight_, int BlockWidth_>
      struct BCSRPerspectiveHelper<DT_, BlockHeight_, BlockWidth_, Perspective::pod>
      {
        typedef DT_ Type;
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief CSR based blocked sparse matrix.
     *
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse row format.
     * Every non zero element is a small dense matrix block on itself, represented as Tiny::Matrix objects.
     * To be consistent with the external interface, the layout information are stored in block-scoped coordinates.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: column index per non zero element \n
     * _indices[1]: row start index (including matrix end index)\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used ValueType elements)\n
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
    class SparseMatrixBCSR : public Container<DT_, IT_>
    {
      static_assert(BlockHeight_ > 0, "invalid block size");
      static_assert(BlockWidth_ > 0, "invalid block size");

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our block height
      static constexpr int BlockHeight = BlockHeight_;
      /// Our block width
      static constexpr int BlockWidth = BlockWidth_;
      /// Value type, meaning the type of each block
      typedef Tiny::Matrix<DataType, BlockHeight, BlockWidth> ValueType;

      /// Our used layout type
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_csr;

      /// ImageIterator typedef for Adjactor interface implementation
      typedef const IT_* ImageIterator;

      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = SparseMatrixBCSR<DT2_, IT2_, BlockHeight_, BlockWidth_>;

      /// this typedef lets you create a matrix container with new Datatape and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

      /// Compatible L-vector type
      typedef typename Intern::BCSRVectorHelper<DT_, IT_, BlockHeight_>::VectorType VectorTypeL;
      /// Compatible R-vector type
      typedef typename Intern::BCSRVectorHelper<DT_, IT_, BlockWidth_>::VectorType VectorTypeR;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

      /**
       * \brief Scatter-Axpy operation for SparseMatrixBCSR
       *
       * Apart from the MatrixType and usage of DataTypeBlocked c&p from SparseMatrixCSR.
       *
       * \author Jordi Paul
       */
      class ScatterAxpy
      {
      public:
        typedef LAFEM::SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_> MatrixType;
        typedef DT_ DataType;
        typedef IT_ IndexType;
        typedef Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> ValueType;

      private:
#ifdef DEBUG
        const IT_ _deadcode;
#endif
        Index _num_rows;
        Index _num_cols;
        IT_* _row_ptr;
        IT_* _col_idx;
        IT_* _col_ptr;
        ValueType* _data;

      public:
        explicit ScatterAxpy(MatrixType& matrix) :
#ifdef DEBUG
          _deadcode(~IT_(0)),
#endif
          _num_rows(matrix.rows()),
          _num_cols(matrix.columns()),
          _row_ptr(matrix.row_ptr()),
          _col_idx(matrix.col_ind()),
          _col_ptr(nullptr),
          _data(matrix.val())
        {
          // allocate column-pointer array
          _col_ptr = new IT_[_num_cols];
#ifdef DEBUG
          for(Index i(0); i < _num_cols; ++i)
          {
            _col_ptr[i] = _deadcode;
          }
#endif
        }

        virtual ~ScatterAxpy()
        {
          if(_col_ptr != nullptr)
          {
            delete [] _col_ptr;
          }
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
            for(IT_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
            {
              _col_ptr[_col_idx[k]] = k;
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
            for(IT_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
            {
              _col_ptr[_col_idx[k]] = _deadcode;
            }
#endif
            // continue with next row entry
          }
        }
      }; // class ScatterAxpy

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

    public:

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit SparseMatrixBCSR() :
        Container<DT_, IT_> (0)
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
      explicit SparseMatrixBCSR(Index rows_in, Index columns_in) :
        Container<DT_, IT_> (rows_in * columns_in)
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
       * \param[in] used_elements_in The amount of non zero elements of the created matrix.
       *
       * Creates an empty (but allocated) matrix.
       *
       * \note The allocated memory will not be initialized.
       */
      explicit SparseMatrixBCSR(Index rows_in, Index columns_in, Index used_elements_in) :
        Container<DT_, IT_> (rows_in * columns_in)
      {
        XASSERT(rows_in != Index(0) && columns_in != Index(0));

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(used_elements_in);

        this->_indices.push_back(MemoryPool::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());

        this->_indices.push_back(MemoryPool::template allocate_memory<IT_>(_rows() + 1));
        this->_indices_size.push_back(_rows() + 1);

        this->_elements.push_back(MemoryPool::template allocate_memory<DT_>(used_elements<Perspective::pod>()));
        this->_elements_size.push_back(used_elements<Perspective::pod>());

      }

      /**
       * \brief Constructor
       *
       * \param[in] layout_in The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixBCSR(const SparseLayout<IT_, layout_id> & layout_in) :
        Container<DT_, IT_> (layout_in._scalar_index.at(0))
      {
        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());

        for (auto i : this->_indices)
          MemoryPool::increase_memory(i);

        this->_elements.push_back(MemoryPool::template allocate_memory<DT_>(used_elements<Perspective::pod>()));
        this->_elements_size.push_back(used_elements<Perspective::pod>());
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix from
       *
       * Creates a BCSR blocked matrix based on a given adjacency graph representing the sparsity pattern.
       */
      explicit SparseMatrixBCSR(const Adjacency::Graph & graph) :
        Container<DT_, IT_>(0)
      {
        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nnze = graph.get_num_indices();

        if (num_nnze == 0)
        {
          this->move(SparseMatrixBCSR(num_rows, num_cols));
          return;
        }

        // Create temporary vectors. Row and column pointer are block wise
        LAFEM::DenseVector<IT_, IT_> vrow_ptr(num_rows+1);
        LAFEM::DenseVector<IT_, IT_> vcol_idx(num_nnze);
        // The data array has to account for the block size
        LAFEM::DenseVector<DT_, IT_> vdata(num_nnze*Index(BlockHeight_*BlockWidth_), DT_(0));

        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());
        IT_ * prow_ptr(vrow_ptr.elements());
        IT_ * pcol_idx(vcol_idx.elements());

        // build row-end
        prow_ptr[0] = IT_(dom_ptr[0]);
        for(Index i(0); i < num_rows; ++i)
          prow_ptr[i+1] = IT_(dom_ptr[i+1]);

        // build col-idx
        for(Index i(0); i < num_nnze; ++i)
          pcol_idx[i] = IT_(img_idx[i]);

        // build the matrix
        this->move(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>(num_rows, num_cols, vcol_idx, vdata, vrow_ptr));
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a BCSR matrix based on the source file.
       */
      explicit SparseMatrixBCSR(FileMode mode, const String& filename) :
        Container< DT_, IT_>(0)
      {
        read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a BCSR matrix based on the source filestream.
       */
      explicit SparseMatrixBCSR(FileMode mode, std::istream& file) :
        Container<DT_, IT_>(0)
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
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixBCSR(const Index rows_in, const Index columns_in,
                                      DenseVector< IT_, IT_> & col_ind_in, DenseVector< DT_, IT_> & val_in, DenseVector< IT_, IT_> & row_ptr_in) :
        Container<DT_, IT_>(rows_in * columns_in)
      {
        /// \todo maybe create empty matrix if col_ind and val and row_ptr inputs are all three empty
        XASSERT(col_ind_in.size() > 0);
        XASSERT(val_in.size() > 0);
        XASSERT(row_ptr_in.size() > 0);

        XASSERT(rows_in != Index(0) && columns_in != Index(0));
        XASSERTM(val_in.size() % (BlockHeight_ * BlockWidth_) == 0, "input values size is not a multiple of container blocksize!");

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(val_in.size() / Index(BlockHeight_ * BlockWidth_));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(col_ind_in.elements());
        this->_indices_size.push_back(col_ind_in.size());
        this->_indices.push_back(row_ptr_in.elements());
        this->_indices_size.push_back(row_ptr_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool::increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool::increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixBCSR(std::vector<char> input) :
        Container<DT_, IT_>(0)
      {
        deserialize<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixBCSR(SparseMatrixBCSR && other) :
        Container<DT_, IT_>(std::forward<SparseMatrixBCSR>(other))
      {
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixBCSR & operator= (SparseMatrixBCSR && other)
      {
        this->move(std::forward<SparseMatrixBCSR>(other));

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
      SparseMatrixBCSR clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        SparseMatrixBCSR t;
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
      void clone(const SparseMatrixBCSR<DT2_, IT2_, BlockHeight_, BlockWidth_> & other, CloneMode clone_mode = CloneMode::Weak)
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
      void convert(const SparseMatrixBCSR<DT2_, IT2_, BlockHeight_, BlockWidth_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] graph The graph to create the matrix from
       *
       * Creates a BCSR blocked matrix based on a given adjacency graph representing the sparsity pattern.
       */
      void convert(const Adjacency::Graph & graph)
      {
        this->move(SparseMatrixBCSR(graph));
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixBCSR & operator= (const SparseLayout<IT_, layout_id> & layout_in)
      {
        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool::release_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool::release_memory(this->_indices.at(i));

        this->_elements.clear();
        this->_indices.clear();
        this->_elements_size.clear();
        this->_indices_size.clear();
        this->_scalar_index.clear();

        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());

        for (auto i : this->_indices)
          MemoryPool::increase_memory(i);

        this->_elements.push_back(MemoryPool::template allocate_memory<DT_>(used_elements<Perspective::pod>()));
        this->_elements_size.push_back(used_elements<Perspective::pod>());

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
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_bcsr, input);
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
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_bcsr, config);
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
          /// \todo read_from_mtx
          /*case FileMode::fm_mtx:
           *   read_from_mtx(filename);
           *   break;*/
          case FileMode::fm_bcsr:
          case FileMode::fm_binary:
            this->template _deserialize<double, std::uint64_t>(FileMode::fm_bcsr, file);
            break;
          default:
            XABORTM("Filemode not supported!");
        }
      }

      /*
       * \brief Read in matrix from MatrixMarket mtx file.
       *
       * \param[in] filename The file that shall be read in.
       */
      /*void read_from_mtx(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);
        read_from_mtx(file);
        file.close();
      }*/

      /*
       * \brief Read in matrix from MatrixMarket mtx stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      /*void read_from_mtx(std::istream& file)
      {
        this->clear();
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);

        itd::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >

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
          _columns() = col;
          _size() = this->rows() * this->columns();
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
        _size() = this->rows() * this->columns();
        _used_elements() = ue;

        DT_ * tval = new DT_[ue];
        IT_ * tcol_ind = new IT_[ue];
        IT_ * trow_ptr = new IT_[rows() + 1];

        IT_ idx(0);
        Index row_idx(0);
        for (auto row : entries)
        {
          trow_ptr[row_idx] = idx;
          for (auto col : row.second )
          {
            tcol_ind[idx] = col.first;
            tval[idx] = col.second;
            ++idx;
          }
          row.second.clear();
          ++row_idx;
        }
        trow_ptr[rows()] = IT_(ue);
        entries.clear();

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(rows() + 1));
        this->_indices_size.push_back(rows() + 1);

        MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
        MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
        MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, rows() + 1);

        delete[] tval;
        delete[] tcol_ind;
        delete[] trow_ptr;
      }*/
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
          buff = new char[LAFEM::FileOutStreamBufferSize];
          file.rdbuf()->pubsetbuf(buff, LAFEM::FileOutStreamBufferSize);
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
          case FileMode::fm_bcsr:
          case FileMode::fm_binary:
            this->template _serialize<double, std::uint64_t>(FileMode::fm_bcsr, file);
            break;
          case FileMode::fm_mtx:
          {
            file << "%%MatrixMarket matrix coordinate real general" << "\n";
            file << this->template rows<Perspective::pod>() << " " << this->template columns<Perspective::pod>() << " " << this->template used_elements<Perspective::pod>() << "\n";

            for (Index row(0) ; row < rows() ; ++row)
            {
              const IT_ end(this->row_ptr()[row + 1]);
              for (IT_ i(this->row_ptr()[row]) ; i < end ; ++i)
              {
                auto block = this->val()[i];
                for (int y(0) ; y < BlockHeight_ ; ++y)
                {
                  for (int x(0) ; x < BlockWidth_ ; ++x)
                  {
                    file << ((int)row * BlockHeight_) + y + 1 << " " << ((int)this->col_ind()[i] * BlockWidth_) + x + 1 << " " << std::scientific << block[y][x] << "\n";
                  }
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
       * \brief Retrieve specific matrix element.
       *
       * \param[in] row The row of the matrix element.
       * \param[in] col The column of the matrix element.
       *
       * \returns Specific matrix element.
       */
      ValueType operator()(Index row, Index col) const
      {
        ASSERT(row < rows());
        ASSERT(col < columns());

        MemoryPool::synchronize();

        for (Index i(this->row_ptr()[row]) ; i < this->row_ptr()[row+1] ; ++i)
        {
          if (this->col_ind()[i] == col)
          {
            return this->val<Perspective::native>()[i];
          }
          if (this->col_ind()[i] > col)
            break; //return zero element
        }

        return ValueType(0.);
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
       * \returns Matrix row count if perspective_ = false, e.g. count every block as one row.
       * \returns Raw matrix row count if perspective_ = true, e.g. row_count * BlockHeight_.
       */
      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        if (perspective_ == Perspective::pod)
          return this->_scalar_index.at(1) * Index(BlockHeight_);
        else
          return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count if perspective_ = false, e.g. count every block as one column.
       * \returns Raw matrix column count if perspective_ = true, e.g. column_count * BlockWidth_.
       */
      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        if (perspective_ == Perspective::pod)
          return this->_scalar_index.at(2) * Index(BlockWidth_);
        else
          return this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count if perspective_ = false, e.g. count every block as one entry.
       * \returns Raw non zero element count if perspective_ = true, e.g. used_elements * BlockHeight_ * BlockWidth_.
       */
      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
      {
        if (perspective_ == Perspective::pod)
          return this->_scalar_index.at(3) * Index(BlockHeight_ * BlockWidth_);
        else
          return this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve column indices array.
       *
       * \returns Column indices array.
       */
      IT_ * col_ind()
      {
        if (this->size() == 0)
          return nullptr;

        return this->_indices.at(0);
      }

      /// \copydoc col_ind()
      /// const version.
      IT_ const * col_ind() const
      {
        if (this->size() == 0)
          return nullptr;

        return this->_indices.at(0);
      }

      /**
       * \brief Retrieve non zero element array.
       *
       * \tparam perspective_ template parameter to choose the return value type
       *
       * \returns Non zero element array if perspective_ = Perspective::native, e.g. treat every block as one block.
       * \returns Raw non zero element array if perspective_ = Perspective::pod, e.g. treat every entry of a block separated.
       */
      template <Perspective perspective_ = Perspective::native>
      auto val() const -> const typename Intern::BCSRPerspectiveHelper<DT_, BlockHeight_, BlockWidth_, perspective_>::Type *
      {
        if (this->size() == 0)
          return nullptr;

        return (const typename Intern::BCSRPerspectiveHelper<DT_, BlockHeight_, BlockWidth_, perspective_>::Type *)(this->_elements.at(0));
      }

      /// \copydoc val()
      /// const version.
      template <Perspective perspective_ = Perspective::native>
      auto val() -> typename Intern::BCSRPerspectiveHelper<DT_, BlockHeight_, BlockWidth_, perspective_>::Type *
      {
        if (this->size() == 0)
          return nullptr;

        return (typename Intern::BCSRPerspectiveHelper<DT_, BlockHeight_, BlockWidth_, perspective_>::Type *)(this->_elements.at(0));
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

      /// \copydoc row_ptr()
      /// const version.
      IT_ const * row_ptr() const
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(1);
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SparseMatrixBCSR";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const SparseMatrixBCSR & x, bool full = false)
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
                const SparseMatrixBCSR & x,
                const SparseMatrixBCSR & y,
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

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);
        Arch::Axpy::value(this->template val<Perspective::pod>(),
            alpha,
            x.template val<Perspective::pod>(),
            y.template val<Perspective::pod>(),
            this->used_elements<Perspective::pod>());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const SparseMatrixBCSR & x, const DT_ alpha)
      {
        XASSERTM(x.rows() == this->rows(), "Row count does not match!");
        XASSERTM(x.columns() == this->columns(), "Column count does not match!");
        XASSERTM(x.used_elements() == this->used_elements(), "Nonzero count does not match!");

        TimeStamp ts_start;
        Statistics::add_flops(this->used_elements<Perspective::pod>());
        Arch::Scale::value(this->template val<Perspective::pod>(), x.template val<Perspective::pod>(), alpha, this->used_elements<Perspective::pod>());
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
        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);
        DT_ result = Arch::Norm2::value(this->template val<Perspective::pod>(),
                                              this->used_elements<Perspective::pod>());
        TimeStamp ts_stop;
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
        XASSERTM(row_norms.size() == this->rows(), "Matrix/Vector dimension mismatch");

        TimeStamp ts_start;
        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::RowNorm::bcsr_norm2(row_norms.template elements<Perspective::pod>(),
          this->template val<Perspective::pod>(),
          col_ind(), row_ptr(), rows(), BlockHeight, BlockWidth);

        TimeStamp ts_stop;
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
        XASSERTM(row_norms.size() == this->rows(), "Matrix/Vector dimension mismatch");

        TimeStamp ts_start;
        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::RowNorm::bcsr_norm2sqr(row_norms.template elements<Perspective::pod>(),
          this->template val<Perspective::pod>(),
          col_ind(), row_ptr(), rows(), BlockHeight, BlockWidth);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Computes the square of the 2-norm for every row, where every row is scaled by a vector
       *
       * \param[out] row_norms
       * For every (scaled) row, this left-vector will contain the square of its 2-norm
       *
       * \param[in] scal
       * The scaling vector
       *
       * This computes
       * \f[
       *    row\_norms_i = \sum_{j=0}^{n-1} scal_j (this_{ij})^2
       * \f]
       * and is used to compute
       * \f[
       *   \mathrm{tr}(B^T \mathrm{diag}(A) B)
       * \f]
       *
       */
      void row_norm2sqr(VectorTypeL& row_norms, const VectorTypeR& scal) const
      {
        XASSERTM(row_norms.size() == this->rows(), "Matrix/Vector dimension mismatch");
        XASSERTM(scal.size() == this->columns(), "Matrix/scalings dimension mismatch");

        TimeStamp ts_start;
        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::RowNorm::bcsr_scaled_norm2sqr(row_norms.template elements<Perspective::pod>(),
           scal.template elements<Perspective::pod>(), this->template val<Perspective::pod>(),
          col_ind(), row_ptr(), rows(), BlockHeight, BlockWidth);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this^\top \f$
       *
       * \return The transposed matrix
       *
       * \note The resulting matrix has transposed block dimensions, too.
       */
      SparseMatrixBCSR<DT_, IT_, BlockWidth_, BlockHeight_> transpose() const
      {
        SparseMatrixBCSR<DT_, IT_, BlockWidth_, BlockHeight_> x_t;
        x_t.transpose(*this);
        return x_t;
      }

      /**
       * \brief Calculate \f$this \leftarrow x^\top \f$
       *
       * \param[in] x The matrix to be transposed.
       */
      void transpose(const SparseMatrixBCSR<DT_, IT_, BlockWidth_, BlockHeight_> & x)
      {
        using XType = SparseMatrixBCSR<DT_, IT_, BlockWidth_, BlockHeight_>;
        if (x.used_elements() == 0)
        {
          this->move(SparseMatrixBCSR(x.rows(), x.columns()));
          return;
        }

        const Index txrows(x.rows());
        const Index txcolumns(x.columns());
        const Index txused_elements(x.used_elements());

        const IT_ * ptxcol_ind(x.col_ind());
        const IT_ * ptxrow_ptr(x.row_ptr());
        const typename XType::ValueType * ptxval(x.val());

        DenseVector<IT_, IT_> tcol_ind(txused_elements);
        DenseVector<DT_, IT_> tval(txused_elements * BlockHeight_ * BlockWidth_);
        DenseVector<IT_, IT_> trow_ptr(txcolumns + 1, IT_(0));

        IT_ * ptcol_ind(tcol_ind.elements());
        ValueType * ptval((ValueType*)tval.elements());
        IT_ * ptrow_ptr(trow_ptr.elements());

        ptrow_ptr[0] = 0;

        for (Index i(0); i < txused_elements; ++i)
        {
          ++ptrow_ptr[ptxcol_ind[i] + 1];
        }

        for (Index i(1); i < txcolumns - 1; ++i)
        {
          ptrow_ptr[i + 1] += ptrow_ptr[i];
        }

        for (Index i(0); i < txrows; ++i)
        {
          for (IT_ k(ptxrow_ptr[i]); k < ptxrow_ptr[i+1]; ++k)
          {
            const IT_ l(ptxcol_ind[k]);
            const IT_ j(ptrow_ptr[l]);
            ptval[j].set_transpose(ptxval[k]);
            ptcol_ind[j] = IT_(i);
            ++ptrow_ptr[l];
          }
        }

        for (Index i(txcolumns); i > 0; --i)
        {
          ptrow_ptr[i] = ptrow_ptr[i - 1];
        }
        ptrow_ptr[0] = 0;

        this->move(SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_> (txcolumns, txrows, tcol_ind, tval, trow_ptr));
      }

      /**
       * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_i\f$
       *
       * \param[in] x The matrix whose rows are to be scaled.
       * \param[in] s The vector to the scale the rows by.
       */
      void scale_rows(const SparseMatrixBCSR & x, const DenseVectorBlocked<DT_,IT_, BlockHeight_> & s)
      {
        XASSERTM(x.rows() == this->rows(), "Row count does not match!");
        XASSERTM(x.columns() == this->columns(), "Column count does not match!");
        XASSERTM(x.used_elements() == this->used_elements(), "Nonzero count does not match!");
        XASSERTM(s.size() == this->rows(), "Vector size does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements<Perspective::pod>());
        Arch::ScaleRows::template bcsr<BlockHeight_, BlockWidth_>(this->val<Perspective::pod>(), x.template val<Perspective::pod>(),
          this->col_ind(), this->row_ptr(), s.template elements<Perspective::pod>(), this->rows(), this->columns(), this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      void scale_cols(const SparseMatrixBCSR & x, const DenseVectorBlocked<DT_,IT_, BlockWidth_> & s)
      {
        XASSERTM(x.rows() == this->rows(), "Row count does not match!");
        XASSERTM(x.columns() == this->columns(), "Column count does not match!");
        XASSERTM(x.used_elements() == this->used_elements(), "Nonzero count does not match!");
        XASSERTM(s.size() == this->columns(), "Vector size does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements<Perspective::pod>());
        Arch::ScaleCols::template bcsr<BlockHeight_, BlockWidth_>(this->val<Perspective::pod>(), x.template val<Perspective::pod>(),
          this->col_ind(), this->row_ptr(), s.template elements<Perspective::pod>(), this->rows(), this->columns(), this->used_elements());

        TimeStamp ts_stop;
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
      void apply(DenseVector<DT_, IT_> & r, const DenseVector< DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->rows<Perspective::pod>(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns<Perspective::pod>(), "Vector size of x does not match!");

        TimeStamp ts_start;
        FEAT_KERNEL_MARKER_START("BCSR-apply");

        if (this->used_elements() == 0)
        {
          r.format();
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);
        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.elements(), DT_(1), x.elements(), DT_(0), r.elements(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \attention r and x must \b not refer to the same vector object!
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<DT_, IT_, BlockHeight_> & r, const DenseVector< DT_, IT_> & x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns<Perspective::pod>(), "Vector size of x does not match!");

        TimeStamp ts_start;

        FEAT_KERNEL_MARKER_START("BCSR-apply");
        if (this->used_elements() == 0)
        {
          r.format();
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), DT_(1), x.elements(), DT_(0), r.template elements<Perspective::pod>(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \attention r and x must \b not refer to the same vector object!
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<DT_, IT_> & r, const DenseVectorBlocked<DT_, IT_, BlockWidth_> & x) const
      {
        XASSERTM(r.size() == this->rows<Perspective::pod>(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        TimeStamp ts_start;

        FEAT_KERNEL_MARKER_START("BCSR-apply");
        if (this->used_elements() == 0)
        {
          r.format();
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.elements(), DT_(1), x.template elements<Perspective::pod>(), DT_(0), r.elements(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
        Statistics::add_time_blas2(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \attention r and x must \b not refer to the same vector object!
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<DT_, IT_, BlockHeight_> & r, const DenseVectorBlocked<DT_, IT_, BlockWidth_> & x) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");

        TimeStamp ts_start;

        FEAT_KERNEL_MARKER_START("BCSR-apply");
        if (this->used_elements() == 0)
        {
          r.format();
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), DT_(1), x.template elements<Perspective::pod>(), DT_(0), r.template elements<Perspective::pod>(), this->template val<Perspective::pod>(),
            this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
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
                 const DenseVector< DT_, IT_> & x,
                 const DenseVector< DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows<Perspective::pod>(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns<Perspective::pod>(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows<Perspective::pod>(), "Vector size of y does not match!");

        TimeStamp ts_start;

        FEAT_KERNEL_MARKER_START("BCSR-apply");
        if (this->used_elements() == 0 || Math::abs(alpha) < Math::eps<DT_>())
        {
          r.copy(y);
          //r.scale(beta);
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops( (this->used_elements<Perspective::pod>() + this->rows<Perspective::pod>()) * 2 );

        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.elements(), alpha, x.elements(), DT_(1), y.elements(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
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
                 DenseVectorBlocked<DT_, IT_, BlockHeight_> & r,
                 const DenseVector< DT_, IT_> & x,
                 const DenseVectorBlocked<DT_, IT_, BlockHeight_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns<Perspective::pod>(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

        TimeStamp ts_start;

        FEAT_KERNEL_MARKER_START("BCSR-apply");
        if (this->used_elements() == 0 || Math::abs(alpha) < Math::eps<DT_>())
        {
          r.copy(y);
          //r.scale(beta);
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops( (this->used_elements<Perspective::pod>() + this->rows<Perspective::pod>()) * 2 );
        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), alpha, x.elements(), DT_(1), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
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
                 const DenseVectorBlocked<DT_, IT_, BlockWidth_> & x,
                 const DenseVector< DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows<Perspective::pod>(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows<Perspective::pod>(), "Vector size of y does not match!");

        TimeStamp ts_start;

        FEAT_KERNEL_MARKER_START("BCSR-apply");
        if (this->used_elements() == 0 || Math::abs(alpha) < Math::eps<DT_>())
        {
          r.copy(y);
          //r.scale(beta);
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops( (this->used_elements<Perspective::pod>() + this->rows<Perspective::pod>()) * 2 );
        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.elements(), alpha, x.template elements<Perspective::pod>(), DT_(1), y.elements(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
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
                 DenseVectorBlocked<DT_, IT_, BlockHeight_> & r,
                 const DenseVectorBlocked<DT_, IT_, BlockWidth_> & x,
                 const DenseVectorBlocked<DT_, IT_, BlockHeight_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows(), "Vector size of y does not match!");

        TimeStamp ts_start;

        FEAT_KERNEL_MARKER_START("BCSR-apply");
        if (this->used_elements() == 0 || Math::abs(alpha) < Math::eps<DT_>())
        {
          r.copy(y);
          //r.scale(beta);
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops( (this->used_elements<Perspective::pod>() + this->rows<Perspective::pod>()) * 2 );
        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), alpha, x.template elements<Perspective::pod>(), DT_(1), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(),
            this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
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
                 DenseVectorBlocked<DT_, IT_, BlockHeight_> & r,
                 const DenseVectorBlocked<DT_, IT_, BlockWidth_> & x,
                 const DenseVector< DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        XASSERTM(r.size() == this->rows(), "Vector size of r does not match!");
        XASSERTM(x.size() == this->columns(), "Vector size of x does not match!");
        XASSERTM(y.size() == this->rows<Perspective::pod>(), "Vector size of y does not match!");

        TimeStamp ts_start;

        FEAT_KERNEL_MARKER_START("BCSR-apply");
        if (this->used_elements() == 0 || Math::abs(alpha) < Math::eps<DT_>())
        {
          r.convert(y);
          //r.scale(beta);
          return;
        }

        XASSERTM(r.template elements<Perspective::pod>() != x.template elements<Perspective::pod>(), "Vector x and r must not share the same memory!");

        Statistics::add_flops( (this->used_elements<Perspective::pod>() + this->rows<Perspective::pod>()) * 2 );
        Arch::Apply::template bcsr<BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), alpha, x.template elements<Perspective::pod>(), DT_(1), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(),
            this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());

        FEAT_KERNEL_MARKER_STOP("BCSR-apply");
        TimeStamp ts_stop;
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
        const LAFEM::SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& d,
        const LAFEM::SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& a,
        const LAFEM::SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& b,
        const DT_ alpha = DT_(1),
        const bool allow_incomplete = false)
      {
        // validate matrix dimensions
        XASSERT(BlockHeight_ == BlockWidth_);
        XASSERT(this->rows() == d.rows());
        XASSERT(d.columns() == a.rows());
        XASSERT(a.columns() == b.rows());
        XASSERT(b.columns() == this->columns());

        // fetch matrix arrays:
        ValueType* data_x = this->val();
        const ValueType* data_d = d.val();
        const ValueType* data_a = a.val();
        const ValueType* data_b = b.val();
        const IT_* row_ptr_x = this->row_ptr();
        const IT_* col_idx_x = this->col_ind();
        const IT_* row_ptr_d = d.row_ptr();
        const IT_* col_idx_d = d.col_ind();
        const IT_* row_ptr_a = a.row_ptr();
        const IT_* col_idx_a = a.col_ind();
        const IT_* row_ptr_b = b.row_ptr();
        const IT_* col_idx_b = b.col_ind();

        // loop over all rows of D and X, resp.
        for(IT_ i(0); i < IT_(this->rows()); ++i)
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
              ValueType omega;
              omega.set_mat_mat_mult(data_d[ik], data_a[kl]);
              omega *= alpha;

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
                  ValueType temp;
                  temp.set_mat_mat_mult(omega, data_b[lj]);
                  data_x[ij] += temp;
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
       * \note
       * This version of the add_double_mat_product operation works on d and b as scalar csr matrices.
       * Their values are treaded as beeing a fully filled block.
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
        const LAFEM::SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& a,
        const LAFEM::SparseMatrixCSR<DT_, IT_>& b,
        const DT_ alpha = DT_(1),
        const bool allow_incomplete = false)
      {
        // validate matrix dimensions
        XASSERT(BlockHeight_ == BlockWidth_);
        XASSERT(this->rows() == d.rows());
        XASSERT(d.columns() == a.rows());
        XASSERT(a.columns() == b.rows());
        XASSERT(b.columns() == this->columns());

        // fetch matrix arrays:
        ValueType* data_x = this->val();
        const DT_* data_d = d.val();
        const ValueType* data_a = a.val();
        const DT_* data_b = b.val();
        const IT_* row_ptr_x = this->row_ptr();
        const IT_* col_idx_x = this->col_ind();
        const IT_* row_ptr_d = d.row_ptr();
        const IT_* col_idx_d = d.col_ind();
        const IT_* row_ptr_a = a.row_ptr();
        const IT_* col_idx_a = a.col_ind();
        const IT_* row_ptr_b = b.row_ptr();
        const IT_* col_idx_b = b.col_ind();

        // loop over all rows of D and X, resp.
        for(IT_ i(0); i < IT_(this->rows()); ++i)
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
              ValueType omega;
              omega.set_mat_mat_mult(ValueType(data_d[ik]), data_a[kl]);
              omega *= alpha;

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
                  ValueType temp;
                  temp.set_mat_mat_mult(omega, ValueType(data_b[lj]));
                  data_x[ij] += temp;
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
      ///@}

      /// \copydoc lump_rows()
      void lump_rows(VectorTypeL& lump) const
      {
        XASSERTM(lump.size() == rows(), "lump vector size does not match matrix row count!");

        Arch::Lumping::bcsr(
          lump.template elements<Perspective::pod>(),
          this->template val<Perspective::pod>(),
          col_ind(), row_ptr(), rows(), BlockHeight, BlockWidth);
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
       * \param[in] diag_indices
       * A vector containing the indices of the diagonal entries
       *
       * \param[out] diag
       * The vector containing the diagonal entry values
       */
      void extract_diag(VectorTypeL & diag, DenseVector<IT_, IT_> & diag_indices) const
      {
        XASSERTM(diag.size() == rows(), "diag size does not match matrix row count!");
        XASSERTM(rows() == columns(), "matrix is not square!");

        for (Index row(0); row < rows(); row++)
        {
          const Index index = diag_indices.elements()[row];
          typename VectorTypeL::ValueType t(0);
          if (index != used_elements())
          {
            ValueType m = this->template val<LAFEM::Perspective::native>()[index];
            for (int i(0) ; i < BlockHeight_ ; ++i)
            {
              t[i] = m[i][i];
            }
          }
          diag(row, t);
        }
      }

      /**
       * \brief extract main diagonal vector from matrix
       *
       * \param[out] diag
       * The vector containing the diagonal entry values
       */
      void extract_diag(VectorTypeL & diag) const
      {
        auto diag_indices = extract_diag_indices();
        extract_diag(diag, diag_indices);
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

      /**
       * \brief extract main diagonal vector from matrix
       *
       * \param[out] diag_indices
       * A vector containing the indices of the diagonal entries
       */
      void extract_diag_indices(DenseVector<IT_, IT_> & diag_indices) const
      {
        XASSERTM(diag_indices.size() == rows(), "diag size does not match matrix row count!");
        XASSERTM(rows() == columns(), "matrix is not square!");

        Arch::Diagonal::csr(diag_indices.elements(), col_ind(), row_ptr(), rows());
      }

      /**
       * \brief extract main diagonal vector from matrix
       *
       * \returns A vector containing the indices of the diagonal entries
       */
      DenseVector<IT_, IT_> extract_diag_indices() const
      {
        DenseVector<IT_, IT_> diag_indices(this->template rows<LAFEM::Perspective::native>());
        extract_diag_indices(diag_indices);
        return diag_indices;
      }

      /// Permutate matrix rows and columns according to the given Permutations
      void permute(Adjacency::Permutation & perm_row, Adjacency::Permutation & perm_col)
      {
        if (perm_row.size() == 0 && perm_col.size() == 0)
          return;

        XASSERTM(perm_row.size() == this->rows(), "Container rows does not match permutation size");
        XASSERTM(perm_col.size() == this->columns(), "Container columns does not match permutation size");

        // http://de.mathworks.com/help/matlab/math/sparse-matrix-operations.html#f6-13070
        IT_ * temp_row_ptr = new IT_[rows() + 1];
        IT_ * temp_col_ind = new IT_[used_elements()];
        ValueType * temp_val = new ValueType[used_elements()];

        Index * perm_pos;
        perm_pos = perm_row.get_perm_pos();

        //permute rows from source to temp_*
        Index new_start(0);
        temp_row_ptr[0] = 0;
        for (Index row(0) ; row < this->rows() ; ++row)
        {
          Index row_size(this->row_ptr()[perm_pos[row] + 1] - this->row_ptr()[perm_pos[row]]);

          //iterate over all elements in single one new and old row
          for (Index i(new_start), j(this->row_ptr()[perm_pos[row]]) ; i < new_start + row_size ; ++i, ++j)
          {
            temp_col_ind[i] = this->col_ind()[j];
            temp_val[i] = this->val()[j];
          }

          new_start += row_size;
          temp_row_ptr[row+1] = (IT_)new_start;
        }

        //use inverse col permutation as lookup table: i -> new location of i
        Adjacency::Permutation perm_col_inv = perm_col.inverse();
        perm_pos = perm_col_inv.get_perm_pos();

        //permute columns from temp_* to source
        ::memcpy(this->row_ptr(), temp_row_ptr, (rows() + 1) * sizeof(IT_));
        ::memcpy(this->val(), temp_val, used_elements() * sizeof(ValueType));
        for (Index i(0) ; i < used_elements() ; ++i)
        {
          this->col_ind()[i] = (IT_)perm_pos[temp_col_ind[i]];
        }

        delete[] temp_row_ptr;
        delete[] temp_col_ind;
        delete[] temp_val;

        //sort columns in every row by column index
        IT_ swap_key;
        ValueType swap_val;
        for (Index row(0) ; row < rows() ; ++row)
        {
          Index offset(this->row_ptr()[row]);
          Index row_size(this->row_ptr()[row+1] - this->row_ptr()[row]);
          for (Index i(1), j ; i < row_size ; ++i)
          {
            swap_key = this->col_ind()[i + offset];
            swap_val = this->val()[i + offset];
            j = i;
            while (j > 0 && this->col_ind()[j - 1 + offset] > swap_key)
            {
              this->col_ind()[j + offset] = this->col_ind()[j - 1 + offset];
              this->val()[j + offset] = this->val()[j - 1 + offset];
              --j;
            }
            this->col_ind()[j + offset] = swap_key;
            this->val()[j + offset] = swap_val;
          }
        }
      }

      /**
       * \brief Adds the trace of a double matrix multiplikation to a vector
       *
       * This function performs the following computation:
       * \f[ v \leftarrow v + \alpha \mathrm{tr}( D \textnormal{diag}(A) \cdot B\f]
       *
       * where \f$ B = \mathrm{this} \f$
       *
       * where
       * - \e v denotes a m vector
       * - \e D denotes a m-by-l matrix
       * - \e A denotes a vector representing a l-by-l diagonal matrix
       * - \e B denotes a l-by-n matrix
       *
       * and the block sizes have to match accordingly:
       *  - \e v has a bm blocks
       *  - \e D has bm-by-bl blocks
       *  - \e A has bl blocks
       *  - \e B has bl-by-bn blocks
       *
       * \note
       * This function currently only supports data in main memory.
       *
       * \tparam[BHD_]
       * The block height of D, which is also the block size for v
       *
       * \param[in] a
       * The vector representing the diagonal matrix A.
       *
       * \param[in] d, b
       * The left and right multiplicant matrices
       *
       * \param[in] alpha
       * The scaling factor for the product.
       *
       * \author Jordi Paul
       *
       */
      template<int BHD_>
      void add_trace_double_mat_mult(
        typename LAFEM::SparseMatrixBCSR< DT_, IT_, BHD_, BlockHeight>::VectorTypeL& v,
        const LAFEM::SparseMatrixBCSR<DT_, IT_, BHD_, BlockHeight>& d,
        const LAFEM::DenseVectorBlocked<DT_, IT_, BlockHeight>& a,
        const DT_ alpha = DT_(1)) const
      {
        typedef LAFEM::SparseMatrixBCSR<DT_, IT_, BHD_, BlockHeight> MatrixD;

        static constexpr IT_ bhd = IT_(BHD_);
        static constexpr IT_ bwd = IT_(BlockHeight);

        static constexpr IT_ bh = IT_(BlockHeight);
        static constexpr IT_ bw = IT_(BlockWidth);

        const auto& b(*this);

        // Check matrix dimensions
        XASSERT(v.size() == d.rows());
        XASSERT(d.columns() == a.size());
        XASSERT(a.size() == b.rows());
        XASSERT(b.columns() == this->columns());

        // Fetch matrix arrays
        // Use POD here to avoid hassle with 1x1 matrices and stuff
        DT_* data_v = v.template elements<LAFEM::Perspective::pod>();
        const DT_* data_d = d.template val<LAFEM::Perspective::pod>();
        const DT_* data_a = a.template elements<LAFEM::Perspective::pod>();
        const DT_* data_b = b.template val<LAFEM::Perspective::pod>();
        const IT_* row_ptr_d = d.row_ptr();
        const IT_* col_idx_d = d.col_ind();
        const IT_* row_ptr_b = b.row_ptr();
        const IT_* col_idx_b = b.col_ind();

        // Loop over all rows i of v
        for(IT_ row(0); row < IT_(v.size()); ++row)
        {
          // For all non-zeros D_ik in row i of D
          for(IT_ pos_d(row_ptr_d[row]); pos_d  < row_ptr_d[row+1]; ++pos_d)
          {
            // Get column index k
            const IT_ col_d(col_idx_d[pos_d]);

            // Pre-compute block factor (alpha * D_ik * A_kk)
            typename MatrixD::ValueType omega(0);
            for(IT_ i(0); i < bhd; ++i)
            {
              for(IT_ j(0); j < bwd; ++j)
              {
                omega[int(i)][int(j)] = alpha * data_d[bhd*bwd*pos_d + bwd*i + j] * data_a[bwd*col_d + j];
              }
            }

            // Loop over all non-zero blocks B_kj in row j of B and perform a "sparse axpy" of B_l onto v_i, i.e.:
            // v_i += alpha * (D_ik * A_kk) * B_k.
            IT_ pos_b(row_ptr_b[col_d]);
            while(pos_b < row_ptr_b[col_d+1])
            {
              if(col_idx_b[pos_b] == row)
              {
                // B_kj contributes to v_i here
                for(IT_ i(0); i < bhd; ++i)
                {
                  for(IT_ j(0); j < bh; ++j)
                  {
                    data_v[bhd*row + i] += omega[int(i)][int(j)]*data_b[bh*bw*pos_b + bw*j + i];
                  }
                }
                break;
              }
              ++pos_b;
            }
          }
        }
      }

      /**
       * \brief SparseMatrixBCSR comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      friend bool operator== (const SparseMatrixBCSR & a, const SparseMatrixBCSR & b)
      {
        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;
        if (a.used_elements() != b.used_elements())
          return false;

        if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && a.get_indices().size() == 0 && b.get_elements().size() == 0 && b.get_indices().size() == 0)
          return true;

        IT_ * col_ind_a;
        IT_ * col_ind_b;
        DT_ * val_a;
        DT_ * val_b;
        IT_ * row_ptr_a;
        IT_ * row_ptr_b;

        bool ret(true);

        col_ind_a = const_cast<IT_*>(a.col_ind());
        val_a = const_cast<DT_*>(a.template val<Perspective::pod>());
        row_ptr_a = const_cast<IT_*>(a.row_ptr());

        col_ind_b = const_cast<IT_*>(b.col_ind());
        val_b = const_cast<DT_*>(b.template val<Perspective::pod>());
        row_ptr_b = const_cast<IT_*>(b.row_ptr());

        for (Index i(0) ; i < a.used_elements() ; ++i)
        {
          if (col_ind_a[i] != col_ind_b[i])
          {
            ret = false;
            break;
          }
        }

        if (ret)
        {
          for (Index i(0) ; i < a.template used_elements<Perspective::pod>() ; ++i)
          {
            if (val_a[i] != val_b[i])
            {
              ret = false;
              break;
            }
          }
        }

        if (ret)
        {
          for (Index i(0) ; i < a.rows() + 1; ++i)
          {
            if (row_ptr_a[i] != row_ptr_b[i])
            {
              ret = false;
              break;
            }
          }
        }

        return ret;
      }

      /**
       * \brief SparseMatrixBCSR streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const SparseMatrixBCSR & b)
      {
        lhs << "[" << "\n";
        for (Index i(0) ; i < b.rows() ; ++i)
        {
          for (int k(0) ; k < BlockHeight_ ; ++k)
          {
            lhs << "[";
            for (Index j(0) ; j < b.columns() ; ++j)
            {
              for (int l(0) ; l < BlockWidth_ ; ++l)
                lhs << "  " << b(i, j).v[k][l];
            }
            lhs << "]" << "\n";
          }
        }
        lhs << "]" << "\n";

        return lhs;
      }

      /// \cond internal
      // Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(this->rows());
      }

      // Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(this->columns());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        const auto * prow_ptr(this->row_ptr());
        const auto trow = row / Index(BlockHeight_);
        return Index(prow_ptr[trow + 1] - prow_ptr[trow]) * Index(BlockWidth_);
      }

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DT_ * const pval_set, IT_ * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        const auto * prow_ptr(this->row_ptr());
        const auto * pcol_ind(this->col_ind());
        const auto * pval(this->val());

        const auto trow = int(row) / BlockHeight_;
        const auto lrow = int(row) - trow * BlockHeight_;

        const Index start((Index(prow_ptr[trow])));
        const Index end((Index(prow_ptr[trow + 1] - prow_ptr[trow])));
        for (Index i(0); i < end; ++i)
        {
          for (Index ti(0); ti < Index(BlockWidth_); ++ti)
          {
            pval_set[(i * Index(BlockWidth_) + ti) * stride] = pval[start + i](lrow, int(ti));
            pcol_set[(i * Index(BlockWidth_) + ti) * stride] = pcol_ind[start + i] * IT_(BlockWidth_) + IT_(ti) + IT_(col_start);
          }
        }
      }

      void set_line_reverse(const Index row, const DT_ * const pval_set, const Index stride = 1)
      {
        const auto * prow_ptr(this->row_ptr());
        auto * pval(this->val());

        const auto trow = int(row) / BlockHeight_;
        const auto lrow = int(row) - trow * BlockHeight_;

        const Index start((Index(prow_ptr[trow])));
        const Index end((Index(prow_ptr[trow + 1] - prow_ptr[trow])));
        for (Index i(0); i < end; ++i)
        {
          for (Index ti(0); ti < Index(BlockWidth_); ++ti)
          {
            pval[start + i](lrow, int(ti)) = pval_set[(i * Index(BlockWidth_) + ti) * stride];
          }
        }
      }
      /// \endcond

      /* ******************************************************************* */
      /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
      /* ******************************************************************* */
    public:
      /** \copydoc Adjactor::get_num_nodes_domain() */
      inline Index get_num_nodes_domain() const
      {
        return rows();
      }

      /** \copydoc Adjactor::get_num_nodes_image() */
      inline Index get_num_nodes_image() const
      {
        return columns();
      }

      /** \copydoc Adjactor::image_begin() */
      inline ImageIterator image_begin(Index domain_node) const
      {
        XASSERTM(domain_node < rows(), "Domain node index out of range");
        return &this->_indices.at(0)[this->_indices.at(1)[domain_node]];
      }

      /** \copydoc Adjactor::image_end() */
      inline ImageIterator image_end(Index domain_node) const
      {
        XASSERTM(domain_node < rows(), "Domain node index out of range");
        return &this->_indices.at(0)[this->_indices.at(1)[domain_node + 1]];
      }
    }; // class SparseMatrixBCSR
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_SPARSE_MATRIX_BCSR_HPP
