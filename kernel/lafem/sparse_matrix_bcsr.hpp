#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_BCSR_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_BCSR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>

#include <fstream>

namespace FEAST
{
  namespace LAFEM
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Mem_, typename DT_, typename IT_, int size_>
      struct BCSRVectorHelper
      {
        static_assert(size_ > 1, "invalid block size");
        typedef DenseVectorBlocked<Mem_, DT_, IT_, size_> VectorType;
      };

      template<typename Mem_, typename DT_, typename IT_>
      struct BCSRVectorHelper<Mem_, DT_, IT_, 1>
      {
        typedef DenseVector<Mem_, DT_, IT_> VectorType;
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
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
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
     * _scalar_index[3]: non zero element count (used elements) (multiple of blocksize)\n
     * _scalar_dt[0]: zero element
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
    class SparseMatrixBCSR : public Container<Mem_, DT_, IT_>
    {
      static_assert(BlockHeight_ > 0, "invalid block size");
      static_assert(BlockWidth_ > 0, "invalid block size");

    public:
      /// Our memory architecture type
      typedef Mem_ MemType;
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

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = class SparseMatrixBCSR<Mem2_, DT2_, IT2_, BlockHeight_, BlockWidth_>;

      /// Compatible L-vector type
      typedef typename Intern::BCSRVectorHelper<Mem_, DT_, IT_, BlockHeight_>::VectorType VectorTypeL;
      /// Compatible R-vector type
      typedef typename Intern::BCSRVectorHelper<Mem_, DT_, IT_, BlockWidth_>::VectorType VectorTypeR;

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
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockHeight_, BlockWidth_> MatrixType;
        typedef Mem::Main MemType;
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
        ValueType *_data;

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
              for(IT_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
              {
                _col_ptr[_col_idx[k]] = k;
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

#ifdef DEBUG
                  // ensure that the column pointer is valid for this index
                  ASSERT(_col_ptr[jx] != _deadcode, "invalid column index");
#endif

                  // incorporate data into global matrix
                  _data[_col_ptr[jx]] += (iw * jw) * loc_mat[i][j];

                  // continue with next column contribution
                }
                // continue with next column entry
              }

#ifdef DEBUG
              // reformat column-pointer array
              for(IT_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
              {
                _col_ptr[_col_idx[k]] = _deadcode;
              }
#endif
              // continue with next row contribution
            }
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
        Container<Mem_, DT_, IT_> (0)
      {
        CONTEXT("When creating SparseMatrixBCSR");
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
      explicit SparseMatrixBCSR(const SparseLayout<Mem_, IT_, layout_id> & layout_in) :
        Container<Mem_, DT_, IT_> (layout_in._scalar_index.at(0))
      {
        CONTEXT("When creating SparseMatrixBCSR");
        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          MemoryPool<Mem_>::increase_memory(i);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(used_elements<Perspective::pod>()));
        this->_elements_size.push_back(used_elements<Perspective::pod>());
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The.graph to create the matrix from
       *
       * Creates a CSR blocked matrix based on a given adjacency graph representing the sparsity pattern.
       */
      explicit SparseMatrixBCSR(const Adjacency::Graph & graph) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixCSR");

        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nnze = graph.get_num_indices();

        // Create temporary vectors. Row and column pointer are block wise
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vrow_ptr(num_rows+1);
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vcol_idx(num_nnze);
        // The data array has to account for the block size
        LAFEM::DenseVector<Mem::Main, DT_, IT_> vdata(num_nnze*Index(BlockHeight_*BlockWidth_), DT_(0));

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
        this->assign(SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockHeight_, BlockWidth_>(num_rows, num_cols, vcol_idx, vdata, vrow_ptr));
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a BCSR matrix based on the source file.
       */
      explicit SparseMatrixBCSR(FileMode mode, String filename) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixBCSR");

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
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixBCSR");

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
                                      DenseVector<Mem_, IT_, IT_> & col_ind_in, DenseVector<Mem_, DT_, IT_> & val_in, DenseVector<Mem_, IT_, IT_> & row_ptr_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        CONTEXT("When creating SparseMatrixBCSR");
        ASSERT(val_in.size() % (BlockHeight_ * BlockWidth_) == 0, "Error: " + stringify(val_in.size()) + " not multiple of container blocksize!");
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(val_in.size() / Index(BlockHeight_ * BlockWidth_));
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(col_ind_in.elements());
        this->_indices_size.push_back(col_ind_in.size());
        this->_indices.push_back(row_ptr_in.elements());
        this->_indices_size.push_back(row_ptr_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] std::vector<char> A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixBCSR(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixBCSR");
        deserialise<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixBCSR(SparseMatrixBCSR && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixBCSR>(other))
      {
        CONTEXT("When moving SparseMatrixBCSR");
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
        CONTEXT("When moving SparseMatrixBCSR");

        this->move(std::forward<SparseMatrixBCSR>(other));

        return *this;
      }

      InsertWeakClone( SparseMatrixBCSR );

      /** \brief Shallow copy operation
       *
       * Create a shallow copy of itself.
       *
       */
      SparseMatrixBCSR shared() const
      {
        CONTEXT("When sharing SparseMatrixBCSR");
        SparseMatrixBCSR r;
        r.assign(*this);
        return r;
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixBCSR<Mem2_, DT2_, IT2_, BlockHeight_, BlockWidth_> & other)
      {
        CONTEXT("When converting SparseMatrixBCSR");
        this->assign(other);
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixBCSR & operator= (const SparseLayout<Mem_, IT_, layout_id> & layout_in)
      {
        CONTEXT("When assigning SparseMatrixBCSR");

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

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(used_elements<Perspective::pod>()));
        this->_elements_size.push_back(used_elements<Perspective::pod>());

        return *this;
      }

      /**
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] std::vector<char> A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialise(std::vector<char> input)
      {
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_bcsr, input);
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * \param[in] mode FileMode enum, describing the actual container specialisation.
       * \param[out] std::vector<char> A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAST::LAFEM::Container::_serialise for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialise()
      {
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_bcsr);
      }

      /**
       * \brief Read in matrix from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, String filename)
      {
        CONTEXT("When reading in SparseMatrixBCSR");

        switch(mode)
        {
        /// \todo read_from_mtx
        /*case FileMode::fm_mtx:
          read_from_mtx(filename);
          break;*/
        case FileMode::fm_bcsr:
          read_from_bcsr(filename);
          break;
        case FileMode::fm_binary:
          read_from_bcsr(filename);
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
        CONTEXT("When reading in SparseMatrixBCSR");

        switch(mode)
        {
        /*case FileMode::fm_mtx:
          read_from_mtx(file);
          break;*/
        case FileMode::fm_bcsr:
          read_from_bcsr(file);
          break;
        case FileMode::fm_binary:
          read_from_bcsr(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Read in matrix from MatrixMarket mtx file.
       *
       * \param[in] filename The file that shall be read in.
       */
      /*void read_from_mtx(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        read_from_mtx(file);
        file.close();
      }*/

      /**
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
        this->_scalar_dt.push_back(DT_(0));

        std::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >

        Index ue(0);
        String line;
        std::getline(file, line);
        const bool general((line.find("%%MatrixMarket matrix coordinate real general") != String::npos) ? true : false);
        const bool symmetric((line.find("%%MatrixMarket matrix coordinate real symmetric") != String::npos) ? true : false);

        if (symmetric == false && general == false)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a compatible mtx-file");
        }

        while(!file.eof())
        {
          std::getline(file,line);
          if (file.eof())
            throw InternalError(__func__, __FILE__, __LINE__, "Input-file is empty");

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
       * \brief Read in matrix from binary file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_bcsr(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        read_from_bcsr(file);
        file.close();
      }

      /**
       * \brief Read in matrix from binary stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_bcsr(std::istream& file)
      {
        this->template _deserialise<double, uint64_t>(FileMode::fm_bcsr, file);
      }


      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        CONTEXT("When writing out SparseMatrixBCSR");

        switch(mode)
        {
        case FileMode::fm_bcsr:
          write_out_bcsr(filename);
          break;
        case FileMode::fm_binary:
          write_out_bcsr(filename);
          break;
        case FileMode::fm_mtx:
          write_out_mtx(filename);
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
        CONTEXT("When writing out SparseMatrixBCSR");

        switch(mode)
        {
        case FileMode::fm_bcsr:
          write_out_bcsr(file);
          break;
        case FileMode::fm_binary:
          write_out_bcsr(file);
          break;
        case FileMode::fm_mtx:
          write_out_mtx(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Write out matrix to bcsr binary file.
       *
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out_bcsr(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        write_out_bcsr(file);
        file.close();
      }

      /**
       * \brief Write out matrix to bcsr binary file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_bcsr(std::ostream& file) const
      {
        if (! std::is_same<DT_, double>::value)
          std::cout<<"Warning: You are writing out a bcsr matrix that is not double precision!"<<std::endl;

        this->template _serialise<double, uint64_t>(FileMode::fm_bcsr, file);
      }

      /**
       * \brief Write out matrix to MatrixMarktet mtx file.
       *
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out_mtx(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        write_out_mtx(file);
        file.close();
      }

      /**
       * \brief Write out matrix to MatrixMarktet mtx file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_mtx(std::ostream& file) const
      {
        SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockHeight_, BlockWidth_> temp;
        temp.convert(*this);

        file << "%%MatrixMarket matrix coordinate real general" << std::endl;
        file << temp.template rows<Perspective::pod>() << " " << temp.template columns<Perspective::pod>() << " " << temp.template used_elements<Perspective::pod>() << std::endl;

        for (Index row(0) ; row < rows() ; ++row)
        {
          const IT_ end(temp.row_ptr()[row + 1]);
          for (IT_ i(temp.row_ptr()[row]) ; i < end ; ++i)
          {
            auto block = temp.val()[i];
            for (int y(0) ; y < BlockHeight_ ; ++y)
            {
              for (int x(0) ; x < BlockWidth_ ; ++x)
              {
                file << ((int)row * BlockHeight_) + y + 1 << " " << ((int)temp.col_ind()[i] * BlockWidth_) + x + 1 << " " << std::scientific << block[y][x] << std::endl;
              }
            }
          }
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
      Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> operator()(Index row, Index col) const
      {
        CONTEXT("When retrieving SparseMatrixBCSR element");

        ASSERT(row < rows(), "Error: " + stringify(row) + " exceeds sparse matrix csr row size " + stringify(rows()) + " !");
        ASSERT(col < columns(), "Error: " + stringify(col) + " exceeds sparse matrix csr column size " + stringify(columns()) + " !");

        for (Index i(MemoryPool<Mem_>::get_element(this->_indices.at(1), row)) ; i < MemoryPool<Mem_>::get_element(this->_indices.at(1), row + 1) ; ++i)
        {
          if (MemoryPool<Mem_>::get_element(this->_indices.at(0), i) == col)
          {
            Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> t;
            MemoryPool<Mem_>::download((DT_*)t.v, this->_elements.at(0) + i * Index(BlockHeight_*BlockWidth_), Index(BlockHeight_*BlockWidth_));
            return t;
          }
          if (MemoryPool<Mem_>::get_element(this->_indices.at(0), i) > col)
            break; //return zero element
        }

        Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> ze((DT_(zero_element())));
        return ze;
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
       * \template perspective_ template parameter to choose the return value type
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
        return "SparseMatrixBCSR";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      void copy(const SparseMatrixBCSR & x)
      {
        this->_copy_content(x);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       */
      template <typename Mem2_>
      void copy(const SparseMatrixBCSR<Mem2_, DT_, IT_, BlockHeight_, BlockWidth_> & x)
      {
        this->_copy_content(x);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow y + \alpha x\f$
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
        if (x.rows() != y.rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
        if (x.used_elements() != y.used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

        TimeStamp ts_start;

        // check for special cases
        // r <- x + y
        if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>());
          Arch::Sum<Mem_>::value(this->template val<Perspective::pod>(),
                                 x.template val<Perspective::pod>(),
                                 y.template val<Perspective::pod>(),
                                 this->used_elements<Perspective::pod>());
        }
        // r <- y - x
        else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>());
          Arch::Difference<Mem_>::value(this->template val<Perspective::pod>(),
                                        y.template val<Perspective::pod>(),
                                        x.template val<Perspective::pod>(),
                                        this->used_elements<Perspective::pod>());
        }
        // r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          this->copy(y);
        // r <- y + alpha*x
        else
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);
          Arch::Axpy<Mem_>::dv(this->template val<Perspective::pod>(),
                               alpha,
                               x.template val<Perspective::pod>(),
                               y.template val<Perspective::pod>(),
                               this->used_elements<Perspective::pod>());
        }

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const SparseMatrixBCSR & x, const DT_ alpha)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");

        TimeStamp ts_start;
        Statistics::add_flops(this->used_elements<Perspective::pod>());
        Arch::Scale<Mem_>::value(this->template val<Perspective::pod>(), x.template val<Perspective::pod>(), alpha, this->used_elements<Perspective::pod>());
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
        DT_ result = Arch::Norm2<Mem_>::value(this->template val<Perspective::pod>(),
                                              this->used_elements<Perspective::pod>());
        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));
        return result;
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        if (r.size() != this->rows<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);
        Arch::ProductMatVec<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
          r.elements(), this->template val<Perspective::pod>(), this->col_ind(), this->row_ptr(),
          x.elements(), this->rows(), columns(), used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::ProductMatVec<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
          r.template elements<Perspective::pod>(), this->template val<Perspective::pod>(), this->col_ind(), this->row_ptr(),
          x.elements(), this->rows(), columns(), used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x) const
      {
        if (r.size() != this->rows<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::ProductMatVec<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
          r.elements(), this->template val<Perspective::pod>(), this->col_ind(), this->row_ptr(),
          x.template elements<Perspective::pod>(), this->rows(), columns(), used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r, const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements<Perspective::pod>() * 2);

        Arch::ProductMatVec<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
          r.template elements<Perspective::pod>(), this->template val<Perspective::pod>(), this->col_ind(), this->row_ptr(),
          x.template elements<Perspective::pod>(), this->rows(), columns(), used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        TimeStamp ts_start;

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);

          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.elements(), y.elements(), this->template val<Perspective::pod>(), this->col_ind(), this->row_ptr(),
            x.elements(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*x
        else
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);

          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.elements(), alpha, x.elements(), y.elements(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockHeight_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        TimeStamp ts_start;

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);
          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), x.elements(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          //r.copy(y);
          r.convert(y);
        // r <- y + alpha*x
        else
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);
          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), alpha, x.elements(), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        TimeStamp ts_start;

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);
          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.elements(), y.elements(), this->template val<Perspective::pod>(), this->col_ind(),  this->row_ptr(),
            x.template elements<Perspective::pod>(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          //r.copy(y);
          r.convert(y);
        // r <- y + alpha*x
        else
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);
          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.elements(), alpha, x.template elements<Perspective::pod>(), y.elements(), this->template val<Perspective::pod>(), this->col_ind(),
            this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockHeight_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        TimeStamp ts_start;

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);
          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(), this->col_ind(), this->row_ptr(),
            x.template elements<Perspective::pod>(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*A*x
        else
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);
          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), alpha, x.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(),
            this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockHeight_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockWidth_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        TimeStamp ts_start;

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);
          Arch::Defect<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(), this->col_ind(), this->row_ptr(),
            x.template elements<Perspective::pod>(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.convert(y);
        // r <- y + alpha*A*x
        else
        {
          Statistics::add_flops(this->used_elements<Perspective::pod>() * 3);
          Arch::Axpy<Mem_>::template csrb<DT_, IT_, BlockHeight_, BlockWidth_>(
            r.template elements<Perspective::pod>(), alpha, x.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->template val<Perspective::pod>(),
            this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }
      ///@}


      /// \copydoc extract_diag()
      void extract_diag(VectorTypeL & diag) const
      {
        ASSERT(diag.size() == rows(), "Error: diag size does not match matrix row count!");
        ASSERT(rows() == columns(), "Error: matrix is not square!");

        /// \todo Replace by Arch::kernels
        const Index n = rows();
        for(Index i(0); i < n; ++i)
        {
          auto diag_sub_matrix = (*this)(i, i);
          FEAST::Tiny::Vector<DT_, BlockHeight_> diag_sub_vector;
          for (int j(0) ; j < BlockHeight_ ; ++j)
          {
            diag_sub_vector[j] = diag_sub_matrix[j][j];
          }
          diag(i, diag_sub_vector);
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
       * \brief SparseMatrixBCSR comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_>
      friend bool operator== (const SparseMatrixBCSR & a, const SparseMatrixBCSR<Mem2_, DT_, IT_, BlockHeight_, BlockWidth_> & b)
      {
        CONTEXT("When comparing SparseMatrixBCSRs");

        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;
        if (a.used_elements() != b.used_elements())
          return false;
        if (a.zero_element() != b.zero_element())
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

        if(std::is_same<Mem::Main, Mem_>::value)
        {
          col_ind_a = (IT_*)a.col_ind();
          val_a = (DT_*)a.val();
          row_ptr_a = (IT_*)a.row_ptr();
        }
        else
        {
          col_ind_a = new IT_[a.used_elements()];
          MemoryPool<Mem_>::template download<IT_>(col_ind_a, a.col_ind(), a.used_elements());
          val_a = new DT_[a.template used_elements<Perspective::pod>()];
          MemoryPool<Mem_>::template download<DT_>(val_a, a.template val<Perspective::pod>(), a.template used_elements<Perspective::pod>());
          row_ptr_a = new IT_[a.rows() + 1];
          MemoryPool<Mem_>::template download<IT_>(row_ptr_a, a.row_ptr(), a.rows() + 1);
        }
        if(std::is_same<Mem::Main, Mem2_>::value)
        {
          col_ind_b = (IT_*)b.col_ind();
          val_b = (DT_*)b.val();
          row_ptr_b = (IT_*)b.row_ptr();
        }
        else
        {
          col_ind_b = new IT_[b.used_elements()];
          MemoryPool<Mem2_>::template download<IT_>(col_ind_b, b.col_ind(), b.used_elements());
          val_b = new DT_[b.template used_elements<Perspective::pod>()];
          MemoryPool<Mem2_>::template download<DT_>(val_b, b.template val<Perspective::pod>(), b.template used_elements<Perspective::pod>());
          row_ptr_b = new IT_[b.rows() + 1];
          MemoryPool<Mem2_>::template download<IT_>(row_ptr_b, b.row_ptr(), b.rows() + 1);
        }

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

        if(! std::is_same<Mem::Main, Mem_>::value)
        {
          delete[] col_ind_a;
          delete[] val_a;
          delete[] row_ptr_a;
        }
        if(! std::is_same<Mem::Main, Mem2_>::value)
        {
          delete[] col_ind_b;
          delete[] val_b;
          delete[] row_ptr_b;
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
        lhs << "[" << std::endl;
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
            lhs << "]" << std::endl;
          }
        }
        lhs << "]" << std::endl;

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
      /// \endcond
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_MATRIX_BCSR_HPP
