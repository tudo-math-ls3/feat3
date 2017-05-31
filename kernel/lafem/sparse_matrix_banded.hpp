#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/apply.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/adjacency/graph.hpp>

#include <fstream>
#include <set>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief sparse banded matrix
     *
     * \tparam Mem_ The \ref FEAT::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its diagonal entries
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: vector of offsets (bottom-left diagonal has offset 0,
     *                                 main diagonal has offset rows - 1 and
     *                                 top-right diagonal has offset row + columns - 1)\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_index[4]: number of offsets \n
     * _scalar_dt[0]: zero element
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
     * - The first diagonal is the one at the bottom-left and gets the offset = 1,
     * - the main diaognal has the offset = rows - 1
     * - and the last offset at the top-right has the offset = rows + columns - 1.
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Christoph Lohmann
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseMatrixBanded : public Container<Mem_, DT_, IT_>
    {
    public:
     /**
       * \brief Scatter-Axpy operation for SparseMatrixBanded
       *
       * \author Christoph Lohmann
       */
      class ScatterAxpy
      {
      public:
        typedef LAFEM::SparseMatrixBanded<Mem::Main, DT_, IT_> MatrixType;
        typedef Mem::Main MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
#ifdef DEBUG
        const IT_ _deadcode;
#endif
        Index _num_rows;
        Index _num_cols;
        Index _num_of_offsets;
        IT_* _offsets;
        IT_* _col_ptr;
        DT_ *_data;

      public:
        explicit ScatterAxpy(MatrixType& matrix) :
#ifdef DEBUG
          _deadcode(~IT_(0)),
#endif
          _num_rows(matrix.rows()),
          _num_cols(matrix.columns()),
          _num_of_offsets(matrix.num_of_offsets()),
          _offsets(matrix.offsets()),
          _col_ptr(nullptr),
          _data(matrix.val())
        {
          // allocate column-pointer array
          _col_ptr = new IT_[matrix.columns()];
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
            // loop over all row entry contributations
            for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
            {
              // fetch row entry weight and pre-multiply by alpha
              DT_ iw = alpha * DT_(row_map.get_weight(i, ic));

              // fetch row index
              Index ix = row_map.get_index(i, ic);

              // build column pointer for this row entry contribution
              for(IT_ k(0); k < _num_of_offsets; ++k)
              {
                if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
                {
                  _col_ptr[_offsets[k] + ix + 1 - _num_rows] = IT_(k * _num_rows + ix);
                }
              }

              // loop over all local column entries
              for(int j(0); j < col_map.get_num_local_dofs(); ++j)
              {
                // loop over all column entry contributions
                for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
                {
                  // fetch trial function dof weight
                  DT_ jw = DT_(col_map.get_weight(j, jc));

                  // fetch column index
                  Index jx = col_map.get_index(j, jc);

                  // ensure that the column pointer is valid for this index
                  ASSERTM(_col_ptr[jx] != _deadcode, "invalid column index");

                  // incorporate data into global matrix
                  _data[_col_ptr[jx]] += (iw * jw) * loc_mat[i][j];

                  // continue with next column contribution
                }
                // continue with next column entry
              }

#ifdef DEBUG
              // reformat column-pointer array
              for(IT_ k(0); k < _num_of_offsets; ++k)
              {
                if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
                {
                  _col_ptr[_offsets[k] + ix + 1 - _num_rows] = _deadcode;
                }
              }
#endif
              // continue with next row contribution
            }
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
        typedef LAFEM::SparseMatrixBanded<Mem::Main, DT_, IT_> MatrixType;
        typedef Mem::Main MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
#ifdef DEBUG
        const IT_ _deadcode;
#endif
        Index _num_rows;
        Index _num_cols;
        Index _num_of_offsets;
        const IT_* _offsets;
        IT_* _col_ptr;
        const DT_ *_data;

      public:
        explicit GatherAxpy(const MatrixType& matrix) :
#ifdef DEBUG
          _deadcode(~IT_(0)),
#endif
          _num_rows(matrix.rows()),
          _num_cols(matrix.columns()),
          _num_of_offsets(matrix.num_of_offsets()),
          _offsets(matrix.offsets()),
          _col_ptr(nullptr),
          _data(matrix.val())
        {
          // allocate column-pointer array
          _col_ptr = new IT_[matrix.columns()];
#ifdef DEBUG
          for(Index i(0); i < _num_cols; ++i)
          {
            _col_ptr[i] = _deadcode;
          }
#endif
        }

        virtual ~GatherAxpy()
        {
          if(_col_ptr != nullptr)
          {
            delete [] _col_ptr;
          }
        }

        template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
        void operator()(LocalMatrix_& loc_mat, const RowMapping_& row_map,
                        const ColMapping_& col_map, DT_ alpha = DT_(1))
        {
          // loop over all local row entries
          for(int i(0); i < row_map.get_num_local_dofs(); ++i)
          {
            // loop over all row entry contributations
            for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
            {
              // fetch row index
              Index ix = row_map.get_index(i, ic);

              // build column pointer for this row entry contribution
              for(IT_ k(0); k < _num_of_offsets; ++k)
              {
                if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
                {
                  _col_ptr[_offsets[k] + ix + 1 - _num_rows] = IT_(k * _num_rows + ix);
                }
              }

              // loop over all local column entries
              for(int j(0); j < col_map.get_num_local_dofs(); ++j)
              {
                // clear  accumulation entry
                DT_ dx(DT_(0));

                // loop over all column entry contributions
                for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
                {
                  // fetch column index
                  Index jx = col_map.get_index(j, jc);

                  // ensure that the column pointer is valid for this index
                  ASSERTM(_col_ptr[jx] != _deadcode, "invalid column index");

                  // update accumulator
                  dx += DT_(col_map.get_weight(j, jc)) * _data[_col_ptr[jx]];

                  // continue with next column contribution
                }

                // update local matrix data
                loc_mat[i][j] += (alpha * DT_(row_map.get_weight(i, ic))) * dx;

                // continue with next column entry
              }

#ifdef DEBUG
              // reformat column-pointer array
              for(IT_ k(0); k < _num_of_offsets; ++k)
              {
                if(_offsets[k] + ix + 1 >= _num_rows && _offsets[k] + ix + 1 < 2 * _num_rows)
                {
                  _col_ptr[_offsets[k] + ix + 1 - _num_rows] = _deadcode;
                }
              }
#endif

              // continue with next row contribution
            }
            // continue with next row entry
          }
        }
      }; // class GatherAxpy

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

      Index & _num_of_offsets()
      {
        return this->_scalar_index.at(4);
      }

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Compatible L-vector type
      typedef DenseVector<MemType, DataType, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<MemType, DataType, IT_> VectorTypeR;
      /// Our used layout type
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_banded;
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
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = class SparseMatrixBanded<Mem2_, DT2_, IT2_>;

      /// this typedef lets you create a matrix container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit SparseMatrixBanded() :
        Container<Mem_, DT_, IT_> (0)
      {
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] layout_in The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixBanded(const SparseLayout<Mem_, IT_, layout_id> & layout_in) :
        Container<Mem_, DT_, IT_> (layout_in._scalar_index.at(0))
      {
        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          MemoryPool<Mem_>::increase_memory(i);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(rows() * num_of_offsets()));
        this->_elements_size.push_back(rows() * num_of_offsets());
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] val_in The vector with non zero elements.
       * \param[in] offsets_in The vector of offsets.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixBanded(const Index rows_in, const Index columns_in,
                                  DenseVector<Mem_, DT_, IT_> & val_in,
                                  DenseVector<Mem_, IT_, IT_> & offsets_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        if (val_in.size() != rows_in * offsets_in.size())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Size of values does not match to number of offsets and row count!");
        }

        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);

        Index tused_elements(0);

        for (Index i(0); i < offsets_in.size(); ++i)
        {
          const Index toffset(Index(offsets_in(i)));

          if (toffset + Index(2) > rows_in + columns_in)
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Offset out of matrix!");
          }

          tused_elements += columns_in + Math::min(rows_in, columns_in + rows_in - toffset - 1) - Math::max(columns_in + rows_in - toffset - 1, columns_in);
        }

        this->_scalar_index.push_back(tused_elements);
        this->_scalar_index.push_back(offsets_in.size());
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(offsets_in.elements());
        this->_indices_size.push_back(offsets_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The graph to create the matrix from
       *
       * Creates a matrix based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit SparseMatrixBanded(const Adjacency::Graph & graph) :
        Container<Mem_, DT_, IT_>(0)
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

        DenseVector<Mem::Main, IT_, IT_> toffsets(Index(moffsets.size()));
        auto * ptoffsets = toffsets.elements();

        IT_ idx(0);
        for (auto off : moffsets)
        {
          ptoffsets[idx] = off;
          ++idx;
        }

        DenseVector<Mem_, IT_, IT_> toffsets_mem;
        toffsets_mem.convert(toffsets);
        DenseVector<Mem_, DT_, IT_> tval_mem(Index(moffsets.size()) * num_rows, DT_(0));

        this->assign(SparseMatrixBanded(num_rows, num_cols, tval_mem, toffsets_mem));
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a banded matrix based on the source filestream.
       */
      explicit SparseMatrixBanded(FileMode mode, String filename) :
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
       * Creates a banded matrix based on the source filestream.
       */
      explicit SparseMatrixBanded(FileMode mode, std::istream& file) :
        Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, file);
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixBanded(std::vector<char> input) :
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
      SparseMatrixBanded(SparseMatrixBanded && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixBanded>(other))
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

/// \cond nodoxy
      InsertWeakClone( SparseMatrixBanded )
/// \endcond

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixBanded<Mem2_, DT2_, IT2_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixBanded & operator= (const SparseLayout<Mem_, IT_, layout_id> & layout_in)
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
        this->_scalar_dt.clear();

        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          MemoryPool<Mem_>::increase_memory(i);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(rows() * num_of_offsets()));
        this->_elements_size.push_back(rows() * num_of_offsets());

        return *this;
      }


      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        switch(mode)
        {
        case FileMode::fm_bm:
          write_out_bm(filename);
          break;
        case FileMode::fm_binary:
          write_out_bm(filename);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
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
          write_out_bm(file);
          break;
        case FileMode::fm_binary:
          write_out_bm(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Write out matrix to banded binary file.
       *
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out_bm(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        write_out_bm(file);
        file.close();
      }

      /**
       * \brief Write out matrix to banded binary file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_bm(std::ostream& file) const
      {
        if (! std::is_same<DT_, double>::value)
          std::cout<<"Warning: You are writing out a banded matrix that is not double precision!"<<std::endl;

        this->template _serialise<double, uint64_t>(FileMode::fm_bm, file);
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

        const Index trows(this->_scalar_index.at(1));

        for (Index i(0); i < this->_scalar_index.at(4); ++i)
        {
          const Index toffset(Index(MemoryPool<Mem_>::get_element(this->_indices.at(0), i)));
          if (row + toffset + Index(1) == col + trows)
          {
            return MemoryPool<Mem_>::get_element(this->_elements.at(0), i * trows + row);
          }
        }
        return zero_element();
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
      const Index & rows() const
      {
        return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      template <Perspective = Perspective::native>
      const Index & columns() const
      {
        return this->_scalar_index.at(2);
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
       * \brief Retrieve number of offsets.
       *
       * \returns Number of offsets.
       */
      const Index & num_of_offsets() const
      {
        return this->_scalar_index.at(4);
      }

      /**
       * \brief Retrieve element array.
       *
       * \returns Non zero element array.
       */
      DT_ * val()
      {
        return this->_elements.at(0);
      }

      DT_ const * val() const
      {
        return this->_elements.at(0);
      }

      /**
       * \brief Retrieve offsets array.
       *
       * \returns offsets array.
       */
      IT_ * offsets()
      {
        return this->_indices.at(0);
      }

      IT_ const * offsets() const
      {
        return this->_indices.at(0);
      }

      /**
       * \brief Retrieve non zero element.
       *
       * \returns Non zero element.
       */
      const DT_ zero_element() const
      {
        return this->_scalar_dt.at(0);
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

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      template <typename Mem2_>
      void copy(const SparseMatrixBanded<Mem2_, DT_, IT_> & x, bool full = false)
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
       */
      void axpy(
                const SparseMatrixBanded & x,
                const SparseMatrixBanded & y,
                const DT_ alpha = DT_(1))
      {
        if (x.rows() != y.rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.num_of_offsets() != y.num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.num_of_offsets() != this->num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.used_elements() != y.used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

        if (Math::abs(alpha) < Math::eps<DT_>())
        {
          this->copy(y);
          //y.scale(beta);
          return;
        }

        Arch::Axpy<Mem_>::dv(this->val(), alpha, x.val(), y.val(), this->rows() * this->num_of_offsets());
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const SparseMatrixBanded & x, const DT_ alpha)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.num_of_offsets() != this->num_of_offsets())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_of_offsets do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

        Arch::Scale<Mem_>::value(this->val(), x.val(), alpha, this->rows() * this->num_of_offsets());
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        return Arch::Norm2<Mem_>::value(this->val(), this->rows() * this->num_of_offsets());
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_>& r, const DenseVector<Mem_, DT_, IT_>& x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        if (r.template elements<Perspective::pod>() == x.template elements<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector x and r must not share the same memory!");

        Arch::Apply<Mem_>::banded(r.elements(),
            DT_(1),
            x.elements(),
            DT_(0),
            r.elements(),
            this->val(),
            this->offsets(),
            this->num_of_offsets(),
            this->rows(),
            this->columns());
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(DenseVector<Mem_,DT_, IT_>& r,
                 const DenseVector<Mem_, DT_, IT_>& x,
                 const DenseVector<Mem_, DT_, IT_>& y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        if (r.template elements<Perspective::pod>() == x.template elements<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector x and r must not share the same memory!");

        if (this->used_elements() == 0 || Math::abs(alpha) < Math::eps<DT_>())
        {
          r.copy(y);
          //r.scale(beta);
          return;
        }

        Arch::Apply<Mem_>::banded(r.elements(),
            alpha,
            x.elements(),
            DT_(1),
            y.elements(),
            this->val(),
            this->offsets(),
            this->num_of_offsets(),
            this->rows(),
            this->columns());
      }
      ///@}

      /// \copydoc extract_diag()
      void extract_diag(VectorTypeL & diag) const
      {
        XASSERTM(diag.size() == rows(), "diag size does not match matrix row count!");
        XASSERTM(rows() == columns(), "matrix is not square!");

        for (Index i(0) ; i < num_of_offsets() ; ++i)
        {
          if (MemoryPool<Mem_>::get_element(offsets(), i) == rows() - 1)
          {
            MemoryPool<Mem_>::copy(diag.elements(), val() + i * rows(), rows());
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
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialise(std::vector<char> input)
      {
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_bm, input);
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialise for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialise()
      {
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_bm);
      }

      /**
       * \brief Read in matrix from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, String filename)
      {
        switch(mode)
        {
        case FileMode::fm_bm:
          read_from_bm(filename);
          break;
        case FileMode::fm_binary:
          read_from_bm(filename);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
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
          read_from_bm(file);
          break;
        case FileMode::fm_binary:
          read_from_bm(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Read in matrix from binary file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_bm(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        read_from_bm(file);
        file.close();
      }

      /**
       * \brief Read in matrix from binary stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_bm(std::istream& file)
      {
        this->template _deserialise<double, uint64_t>(FileMode::fm_bm, file);
      }

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

      /// Returns first row-index of the diagonal matching to the offset i
      Index start_offset(const Index i) const
      {
        return Arch::Intern::ApplyBanded::start_offset(i, offsets(), rows(), columns(), num_of_offsets());
      }

      /// Returns last row-index of the diagonal matching to the offset i
      Index end_offset(const Index i) const
      {
        return Arch::Intern::ApplyBanded::end_offset(i, offsets(), rows(), columns(), num_of_offsets());
      }

      /* ******************************************************************* */
      /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
      /* ******************************************************************* */
    public:
      /** \copydoc Adjactor::get_num_nodes_domain() */
      const inline Index & get_num_nodes_domain() const
      {
        return rows();
      }

      /** \copydoc Adjactor::get_num_nodes_image() */
      const inline Index & get_num_nodes_image() const
      {
        return columns();
      }

      /** \copydoc Adjactor::image_begin() */
      inline ImageIterator image_begin(Index domain_node) const
      {
        XASSERT(domain_node < rows());

        const IT_ * toffsets(offsets());
        const IT_ tnum_of_offsets(num_of_offsets());
        const Index trows(rows());

        for(IT_ k(0); k < tnum_of_offsets; ++k)
        {
          if(toffsets[k] + domain_node + 1 >= trows)
            return ImageIterator(k, domain_node + 1 - trows, toffsets);
        }
      }

      /** \copydoc Adjactor::image_end() */
      inline ImageIterator image_end(Index domain_node) const
      {
        XASSERT(domain_node < rows());

        const IT_ * toffsets(offsets());
        const IT_ tnum_of_offsets(num_of_offsets());
        const Index trows(rows());

        for(IT_ k(0); k < tnum_of_offsets; ++k)
        {
          if(toffsets[k] + domain_node + 1 >= 2 * trows)
          {
            return ImageIterator(k, domain_node + 1 - trows, toffsets);
          }
        }
        return ImageIterator(tnum_of_offsets, domain_node + 1 - trows, toffsets);
      }


      /**
       * \brief SparseMatrixBanded comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_> friend bool operator== (const SparseMatrixBanded & a, const SparseMatrixBanded<Mem2_, DT_, IT_> & b)
      {
        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;
        if (a.num_of_offsets() != b.num_of_offsets())
          return false;
        if (a.used_elements() != b.used_elements())
          return false;
        if (a.zero_element() != b.zero_element())
          return false;

        if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && a.get_indices().size() == 0 && b.get_elements().size() == 0 && b.get_indices().size() == 0)
          return true;

        IT_ * offsets_a;
        IT_ * offsets_b;
        DT_ * val_a;
        DT_ * val_b;

        if(std::is_same<Mem::Main, Mem_>::value)
        {
          offsets_a = (IT_*)a.offsets();
          val_a = (DT_*)a.val();
        }
        else
        {
          offsets_a = new IT_[a.num_of_offsets()];
          MemoryPool<Mem_>::template download<IT_>(offsets_a, a.offsets(), a.num_of_offsets());
          val_a = new DT_[a.num_of_offsets() * a.rows()];
          MemoryPool<Mem_>::template download<DT_>(val_a, a.val(), a.num_of_offsets() * a.rows());
        }
        if(std::is_same<Mem::Main, Mem2_>::value)
        {
          offsets_b = (IT_*)b.offsets();
          val_b = (DT_*)b.val();
        }
        else
        {
          offsets_b = new IT_[b.num_of_offsets()];
          MemoryPool<Mem2_>::template download<IT_>(offsets_b, b.offsets(), b.num_of_offsets());
          val_b = new DT_[b.num_of_offsets() * b.rows()];
          MemoryPool<Mem2_>::template download<DT_>(val_b, b.val(), b.num_of_offsets() * b.rows());
        }

        bool ret(true);

        for (Index i(0); i < a.num_of_offsets(); ++i)
        {
          if (offsets_a[i] != offsets_b[i])
          {
            ret = false;
            break;
          }
        }

        for (Index i(0) ; i < a.num_of_offsets() * a.rows() ; ++i)
        {
          if (val_a[i] != val_b[i])
          {
            ret = false;
            break;
          }
        }

        if(! std::is_same<Mem::Main, Mem_>::value)
        {
          delete[] offsets_a;
          delete[] val_a;
        }
        if(! std::is_same<Mem::Main, Mem2_>::value)
        {
          delete[] offsets_b;
          delete[] val_b;
        }

        return ret;
      }

      /**
       * \brief SparseMatrixBanded streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const SparseMatrixBanded & b)
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
    };

  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_SPARSE_MATRIX_BANDED_HPP
