#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/lafem/sparse_matrix_csr_blocked.hpp>
#include <kernel/lafem/matrix_base.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/scale_row_col.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/adjacency/graph.hpp>

#include <fstream>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief CSR based sparse matrix.
     *
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse row format.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: column index per non zero element \n
     * _indices[1]: row start index (including matrix end index)\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_dt[0]: zero element
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseMatrixCSR : public Container<Mem_, DT_, IT_>, public MatrixBase
    {
      private:
        void _read_from_mtx(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
          _read_from_mtx(file);
          file.close();
        }

        void _read_from_mtx(std::istream& file)
        {
          std::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));

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

            entries[(unsigned int)row].insert(std::pair<IT_, DT_>(col, tval));
            ++ue;
            if (symmetric == true && row != col)
            {
              entries[(unsigned int)col].insert(std::pair<IT_, DT_>(row, tval));
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

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(rows() + 1));
          this->_indices_size.push_back(rows() + 1);

          Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, rows() + 1);

          delete[] tval;
          delete[] tcol_ind;
          delete[] trow_ptr;
        }

        void _read_from_csr(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
          _read_from_csr(file);
          file.close();
        }

        void _read_from_csr(std::istream& file)
        {
          uint64_t rows64;
          uint64_t columns64;
          uint64_t elements64;
          file.read((char *)&rows64, sizeof(uint64_t));
          file.read((char *)&columns64, sizeof(uint64_t));
          file.read((char *)&elements64, sizeof(uint64_t));
          Index trows = (Index)rows64;
          Index tcolumns = (Index)columns64;
          Index elements = (Index)elements64;

          _size() = Index(trows * tcolumns);
          this->_scalar_index.push_back(Index(trows));
          this->_scalar_index.push_back(Index(tcolumns));
          this->_scalar_index.push_back(elements);
          this->_scalar_dt.push_back(DT_(0));

          uint64_t * ccol_ind = new uint64_t[elements];
          file.read((char *)ccol_ind, (long)(elements * sizeof(uint64_t)));
          IT_ * tcol_ind = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(elements);
          for (Index i(0) ; i < elements ; ++i)
            tcol_ind[i] = IT_(ccol_ind[i]);
          delete[] ccol_ind;

          uint64_t * crow_ptr = new uint64_t[std::size_t(trows + 1)];
          file.read((char *)crow_ptr, (long)((trows + 1) * sizeof(uint64_t)));
          IT_ * trow_ptr = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(trows + 1);
          for (Index i(0) ; i < trows + 1 ; ++i)
            trow_ptr[i] = IT_(crow_ptr[i]);
          delete[] crow_ptr;

          double * cval = new double[elements];
          file.read((char *)cval, (long)(elements * sizeof(double)));
          DT_ * tval = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(elements);
          for (Index i(0) ; i < elements ; ++i)
            tval[i] = DT_(cval[i]);
          delete[] cval;

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(elements));
          this->_elements_size.push_back(elements);
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(elements));
          this->_indices_size.push_back(elements);
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_rows() + 1));
          this->_indices_size.push_back(_rows() + 1);

          Util::MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval, elements);
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(0), tcol_ind, elements);
          Util::MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(1), trow_ptr, _rows()  + 1);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tval);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcol_ind);
          Util::MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr);
        }

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
        /// Our datatype
        typedef DT_ DataType;
        /// Our indextype
        typedef IT_ IndexType;
        /// Our memory architecture type
        typedef Mem_ MemType;
        /// Compatible L-vector type
        typedef DenseVector<Mem_, DT_, IT_> VectorTypeL;
        /// Compatible R-vector type
        typedef DenseVector<Mem_, DT_, IT_> VectorTypeR;
        /// Our used layout type
        static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_csr;
        /// ImageIterator typedef for Adjactor interface implementation
        typedef const IT_* ImageIterator;
        /// Our 'base' class type
        template <typename Mem2_, typename DT2_, typename IT2_ = IT_>
        using ContainerType = class SparseMatrixCSR<Mem2_, DT2_, IT2_>;


        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit SparseMatrixCSR() :
          Container<Mem_, DT_, IT_> (0)
        {
          CONTEXT("When creating SparseMatrixCSR");
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] layout The layout to be used.
         *
         * Creates an empty matrix with given layout.
         */
        explicit SparseMatrixCSR(const SparseLayout<Mem_, IT_, layout_id> & layout_in) :
          Container<Mem_, DT_, IT_> (layout_in._scalar_index.at(0))
        {
          CONTEXT("When creating SparseMatrixCSR");
          this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
          this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
          this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
          this->_scalar_dt.push_back(DT_(0));

          for (auto i : this->_indices)
            Util::MemoryPool<Mem_>::instance()->increase_memory(i);

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a CSR matrix based on the source matrix.
         */
      template <typename MT_>
        explicit SparseMatrixCSR(const MT_ & other) :
          Container<Mem_, DT_, IT_>(other.size())
        {
          CONTEXT("When creating SparseMatrixCSR");

          convert(other);
        }

        /**
         * \brief Constructor
         *
         * \param[in] graph The.graph to create the matrix from
         *
         * Creates a CSR matrix based on a given adjacency graph, representing the sparsity pattern.
         */
        explicit SparseMatrixCSR(const Adjacency::Graph & graph) :
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseMatrixCSR");

          Index num_rows = graph.get_num_nodes_domain();
          Index num_cols = graph.get_num_nodes_image();
          Index num_nnze = graph.get_num_indices();

          // create temporary vectors
          LAFEM::DenseVector<Mem::Main, IT_, IT_> vrow_ptr(num_rows+1);
          LAFEM::DenseVector<Mem::Main, IT_, IT_> vcol_idx(num_nnze);
          LAFEM::DenseVector<Mem::Main, DT_, IT_> vdata(num_nnze, DT_(0));

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
          this->assign(SparseMatrixCSR<Mem::Main, DT_, IT_>(num_rows, num_cols, vcol_idx, vdata, vrow_ptr));
        }

        /**
         * \brief Constructor
         *
         * \param[in] mode The used file format.
         * \param[in] filename The source file.
         *
         * Creates a CSR matrix based on the source file.
         */
        explicit SparseMatrixCSR(FileMode mode, String filename) :
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseMatrixCSR");

          switch(mode)
          {
            case FileMode::fm_mtx:
              _read_from_mtx(filename);
              break;
            case FileMode::fm_csr:
              _read_from_csr(filename);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
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
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseMatrixCSR");

          switch(mode)
          {
            case FileMode::fm_mtx:
              _read_from_mtx(file);
              break;
            case FileMode::fm_csr:
              _read_from_csr(file);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
        }

        /**
         * \brief Constructor
         *
         * \param[in] rows The row count of the created matrix.
         * \param[in] columns The column count of the created matrix.
         * \param[in] col_ind Vector with column indices.
         * \param[in] val Vector with non zero elements.
         * \param[in] row_ptr Vector with start indices of all rows into the val/col_ind arrays.
         * Note, that this vector must also contain the end index of the last row and thus has a size of row_count + 1.
         *
         * Creates a matrix with given dimensions and content.
         */
        explicit SparseMatrixCSR(const Index rows_in, const Index columns_in,
            DenseVector<Mem_, IT_, IT_> & col_ind_in, DenseVector<Mem_, DT_, IT_> & val_in, DenseVector<Mem_, IT_, IT_> & row_ptr_in) :
          Container<Mem_, DT_, IT_>(rows_in * columns_in)
        {
          CONTEXT("When creating SparseMatrixCSR");
          this->_scalar_index.push_back(rows_in);
          this->_scalar_index.push_back(columns_in);
          this->_scalar_index.push_back(val_in.size());
          this->_scalar_dt.push_back(DT_(0));

          this->_elements.push_back(val_in.elements());
          this->_elements_size.push_back(val_in.size());
          this->_indices.push_back(col_ind_in.elements());
          this->_indices_size.push_back(col_ind_in.size());
          this->_indices.push_back(row_ptr_in.elements());
          this->_indices_size.push_back(row_ptr_in.size());

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            Util::MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            Util::MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source matrix.
         *
         * Moves a given matrix to this matrix.
         */
        SparseMatrixCSR(SparseMatrixCSR && other) :
          Container<Mem_, DT_, IT_>(std::forward<SparseMatrixCSR>(other))
        {
          CONTEXT("When moving SparseMatrixCSR");
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
          CONTEXT("When moving SparseMatrixCSR");

          this->move(std::forward<SparseMatrixCSR>(other));

          return *this;
        }

        /** \brief Clone operation
         *
         * Create a deep copy of itself.
         *
         * \param[in] clone_indices Should we create a deep copy of the index arrays, too ?
         *
         * \return A deep copy of itself.
         *
         */
        SparseMatrixCSR clone(bool clone_indices = false) const
        {
          CONTEXT("When cloning SparseMatrixiCSR");
          SparseMatrixCSR t;
          t.clone(*this, clone_indices);
          return t;
        }

        using Container<Mem_, DT_, IT_>::clone;

        /**
         * \brief Conversion method
         *
         * \param[in] other The source Matrix.
         *
         * Use source matrix content as content of current matrix
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const SparseMatrixCSR<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting SparseMatrixCSR");
          this->assign(other);
        }

        /**
         * \brief Conversion method
         *
         * \param[in] other The source Matrix.
         *
         * Use source matrix content as content of current matrix
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const SparseMatrixCOO<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting SparseMatrixCSR");

          this->clear();

          this->_scalar_index.push_back(other.size());
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_dt.push_back(other.zero_element());

          SparseMatrixCOO<Mem::Main, DT_, IT_> cother;
          cother.convert(other);

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_rows() + 1));
          this->_indices_size.push_back(_rows() + 1);

          DT_ * tval(nullptr);
          IT_ * tcol_ind(nullptr);
          IT_ * trow_ptr(nullptr);
          if (std::is_same<Mem_, Mem::Main>::value)
          {
            tval = this->_elements.at(0);
            tcol_ind = this->_indices.at(0);
            trow_ptr = this->_indices.at(1);
          }
          else
          {
            tval = new DT_[other.used_elements()];
            tcol_ind = new IT_[other.used_elements()];
            trow_ptr = new IT_[other.rows() + 1];
          }

          IT_ ait(0);
          Index current_row(0);
          trow_ptr[current_row] = IT_(0);
          for (Index it(0) ; it < cother.used_elements() ; ++it)
          {
            Index row(cother.row_indices()[it]);
            IT_ column(cother.column_indices()[it]);

            if (current_row < row)
            {
              for (Index i(current_row + 1) ; i < row ; ++i)
              {
                trow_ptr[i] = ait;
              }
              current_row = row;
              trow_ptr[current_row] = ait;
            }
            tval[ait] = cother.val()[it];
            tcol_ind[ait] = column;
            ++ait;
          }
          for (Index i(current_row + 1) ; i < _rows() ; ++i)
          {
            trow_ptr[i] = ait;
          }
          trow_ptr[_rows()] = ait;

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, _rows() + 1);
            delete[] tval;
            delete[] tcol_ind;
            delete[] trow_ptr;
          }
        }

        /**
         * \brief Conversion method
         *
         * \param[in] other The source Matrix.
         *
         * Use source matrix content as content of current matrix
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const SparseMatrixELL<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting SparseMatrixCSR");

          this->clear();

          this->_scalar_index.push_back(other.size());
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_dt.push_back(other.zero_element());

          SparseMatrixELL<Mem::Main, DT_, IT_> cother;
          cother.convert(other);

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_rows() + 1));
          this->_indices_size.push_back(_rows() + 1);

          DT_ * tval(nullptr);
          IT_ * tcol_ind(nullptr);
          IT_ * trow_ptr(nullptr);
          if (std::is_same<Mem_, Mem::Main>::value)
          {
            tval = this->_elements.at(0);
            tcol_ind = this->_indices.at(0);
            trow_ptr = this->_indices.at(1);
          }
          else
          {
            tval = new DT_[used_elements()];
            tcol_ind = new IT_[used_elements()];
            trow_ptr = new IT_[rows() + 1];
          }

          trow_ptr[0] = IT_(0);
          const Index stride(cother.stride());
          IT_ ue(0);
          for (Index row(0) ; row < cother.rows() ; ++row)
          {
            const IT_ * tAj(cother.Aj());
            const DT_ * tAx(cother.Ax());
            tAj += row;
            tAx += row;

            const Index max(cother.Arl()[row]);
            for(Index n(0); n < max ; n++)
            {
              tval[ue] = *tAx;
              tcol_ind[ue] = *tAj;

              tAj += stride;
              tAx += stride;
              ++ue;
            }
            trow_ptr[row + 1] = ue;
          }

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, _rows() + 1);
            delete[] tval;
            delete[] tcol_ind;
            delete[] trow_ptr;
          }
        }

        /**
         * \brief Convertion method
         *
         * \param[in] other The source Matrix.
         *
         * Use source matrix content as content of current matrix
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const SparseMatrixBanded<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting SparseMatrixCSR");

          this->clear();

          this->_scalar_index.push_back(other.size());
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_dt.push_back(other.zero_element());

          SparseMatrixBanded<Mem::Main, DT_, IT_> cother;
          cother.convert(other);

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_rows() + 1));
          this->_indices_size.push_back(_rows() + 1);

          DT_ * tval(nullptr);
          IT_ * tcol_ind(nullptr);
          IT_ * trow_ptr(nullptr);
          if (std::is_same<Mem_, Mem::Main>::value)
          {
            tval = this->_elements.at(0);
            tcol_ind = this->_indices.at(0);
            trow_ptr = this->_indices.at(1);
          }
          else
          {
            tval = new DT_[other.used_elements()];
            tcol_ind = new IT_[other.used_elements()];
            trow_ptr = new IT_[other.rows() + 1];
          }

          trow_ptr[0] = 0;

          #ifdef START_OFFSET
            #warning Overwriting definition of START_OFFSET
            #undef START_OFFSET
          #endif

          #ifdef END_OFFSET
            #error Overwriting definition of END_OFFSET
            #undef END_OFFSET
          #endif

          #define START_OFFSET(j) ((j == Index(-1)) ? crows : ((j == k) ? 0 : crows - coffsets[j] - 1))
          #define END_OFFSET(j) ((j == Index(-1)) ? crows : ((j == cnum_of_offsets) ? 0 : ccolumns + crows - coffsets[j] - 1))

          const DT_ * cval(cother.val());
          const IT_ * coffsets(cother.offsets());
          const Index cnum_of_offsets(cother.num_of_offsets());
          const Index crows(cother.rows());
          const Index ccolumns(cother.columns());

          // Search first offset of the upper triangular matrix
          Index k(0);
          while (k < cnum_of_offsets && coffsets[k] + 1 < crows)
          {
            ++k;
          }

          IT_ ue(0);
          // iteration over all offsets of the lower triangular matrix
          for (Index i(k + 1); i > 0;)
          {
            --i;

            // iteration over all offsets of the upper triangular matrix
            for (Index j(cnum_of_offsets + 1); j > 0;)
            {
              --j;

              // iteration over all rows which contain the offsets between offset i and offset j
              for (Index l(std::max(START_OFFSET(i), END_OFFSET(j))); l < std::min(START_OFFSET(i-1), END_OFFSET(j-1)); ++l)
              {
                for (Index a(i); a < j; ++a)
                {
                  tval[ue] = cval[a * crows + l];
                  tcol_ind[ue] = IT_(l + coffsets[a] + 1 - crows);
                  ++ue;
                }
                trow_ptr[l + 1] = ue;
              }
            }
          }
          #undef START_OFFSET
          #undef END_OFFSET

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, _rows() + 1);
            delete[] tval;
            delete[] tcol_ind;
            delete[] trow_ptr;
          }
        }

        /**
         * \brief Convertion method
         *
         * \param[in] other The source Matrix.
         *
         * Use source matrix content as content of current matrix
         */
        template <typename Mem2_, typename DT2_, typename IT2_, Index BlockHeight_, Index BlockWidth_>
        void convert(const SparseMatrixCSRBlocked<Mem2_, DT2_, IT2_, BlockHeight_, BlockWidth_> & other)
        {
          CONTEXT("When converting SparseMatrixCSR");

          this->clear();

          this->_scalar_index.push_back(other.raw_rows() * other.raw_columns());
          this->_scalar_index.push_back(other.raw_rows());
          this->_scalar_index.push_back(other.raw_columns());
          this->_scalar_index.push_back(other.raw_used_elements());
          this->_scalar_dt.push_back(other.zero_element());

          SparseMatrixCSRBlocked<Mem::Main, DT_, IT_, BlockHeight_, BlockWidth_> cother;
          cother.convert(other);

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_rows() + 1));
          this->_indices_size.push_back(_rows() + 1);

          DT_ * tval(nullptr);
          IT_ * tcol_ind(nullptr);
          IT_ * trow_ptr(nullptr);
          if (std::is_same<Mem_, Mem::Main>::value)
          {
            tval = this->_elements.at(0);
            tcol_ind = this->_indices.at(0);
            trow_ptr = this->_indices.at(1);
          }
          else
          {
            tval = new DT_[used_elements()];
            tcol_ind = new IT_[used_elements()];
            trow_ptr = new IT_[rows() + 1];
          }

          Index ait(0);
          trow_ptr[0] = IT_(0);
          Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *mval(reinterpret_cast<Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *>(cother.val()));
          for (Index orow(0) ; orow < cother.rows() ; ++orow)
          {
            for (Index row(0) ; row < BlockHeight_ ; ++row)
            {
              for(Index ocol(cother.row_ptr()[orow]) ; ocol < cother.row_ptr()[orow + 1] ; ++ocol)
              {
                Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *tbm(mval + ocol);
                for (Index col(0) ; col < BlockWidth_ ; ++col)
                {
                  tval[ait] = (*tbm)(row,col);
                  tcol_ind[ait] = cother.col_ind()[ocol] * (IT_)BlockWidth_ + (IT_)col;
                  ++ait;
                }
              }
              trow_ptr[orow * BlockHeight_ + row + 1] = (IT_)ait;
            }
          }

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, _rows() + 1);
            delete[] tval;
            delete[] tcol_ind;
            delete[] trow_ptr;
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
          CONTEXT("When converting SparseMatrixCSR");

          typename MT_::template ContainerType<Mem::Main, DT_, IT_> ta;
          ta.convert(a);

          const Index arows(ta.rows());
          const Index acolumns(ta.columns());
          const Index aused_elements(ta.used_elements());

          DenseVector<Mem::Main, DT_, IT_> tval(aused_elements);
          DenseVector<Mem::Main, IT_, IT_> tcol_ind(aused_elements);
          DenseVector<Mem::Main, IT_, IT_> trow_ptr(arows + 1);

          DT_ * pval(tval.elements());
          IT_ * pcol_ind(tcol_ind.elements());
          IT_ * prow_ptr(trow_ptr.elements());

          for (Index i(0); i < arows; ++i)
          {
            prow_ptr[i + 1] = IT_(ta.get_length_of_line(i));
          }

          prow_ptr[0] = IT_(0);

          for (Index i(1); i < arows + 1; ++i)
          {
            prow_ptr[i] += prow_ptr[i - 1];
          }

          for (Index i(0); i < arows; ++i)
          {
            ta.set_line(i, pval + prow_ptr[i], pcol_ind + prow_ptr[i], 0);
          }

          SparseMatrixCSR<Mem::Main, DT_, IT_> ta_csr(arows, acolumns, tcol_ind, tval, trow_ptr);
          SparseMatrixCSR<Mem_, DT_, IT_> a_csr;
          a_csr.convert(ta_csr);

          this->assign(a_csr);
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] layout A sparse matrix layout.
         *
         * Assigns a new matrix layout, discarding all old data
         */
        SparseMatrixCSR & operator= (const SparseLayout<Mem_, IT_, layout_id> & layout_in)
        {
          CONTEXT("When assigning SparseMatrixCSR");

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            Util::MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            Util::MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

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
            Util::MemoryPool<Mem_>::instance()->increase_memory(i);

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());

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
          CONTEXT("When writing out SparseMatrixCSR");

          switch(mode)
          {
            case FileMode::fm_csr:
              write_out_csr(filename);
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
          CONTEXT("When writing out SparseMatrixCSR");

          switch(mode)
          {
            case FileMode::fm_csr:
              write_out_csr(file);
              break;
            case FileMode::fm_mtx:
              write_out_mtx(file);
              break;
            default:
                throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
        }

        /**
         * \brief Write out matrix to csr binary file.
         *
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out_csr(String filename) const
        {
          std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
          write_out_csr(file);
          file.close();
        }

        /**
         * \brief Write out matrix to csr binary file.
         *
         * \param[in] file The stream that shall be written to.
         */
        void write_out_csr(std::ostream& file) const
        {
          if (! std::is_same<DT_, double>::value)
            std::cout<<"Warning: You are writing out an csr matrix with less than double precission!"<<std::endl;

          IT_ * tcol_ind = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(this->_indices_size.at(0));
          Util::MemoryPool<Mem_>::template download<IT_>(tcol_ind, this->_indices.at(0), this->_indices_size.at(0));
          uint64_t * ccol_ind = new uint64_t[this->_indices_size.at(0)];
          for (Index i(0) ; i < this->_indices_size.at(0) ; ++i)
            ccol_ind[i] = tcol_ind[i];
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tcol_ind);

          IT_ * trow_ptr = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(this->_indices_size.at(1));
          Util::MemoryPool<Mem_>::template download<IT_>(trow_ptr, this->_indices.at(1), this->_indices_size.at(1));
          uint64_t * crow_ptr = new uint64_t[this->_indices_size.at(1)];
          for (Index i(0) ; i < this->_indices_size.at(1) ; ++i)
            crow_ptr[i] = trow_ptr[i];
          Util::MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr);

          DT_ * tval = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(this->_elements_size.at(0));
          Util::MemoryPool<Mem_>::template download<DT_>(tval, this->_elements.at(0), this->_elements_size.at(0));
          double * cval = new double[this->_elements_size.at(0)];
          for (Index i(0) ; i < this->_elements_size.at(0) ; ++i)
            cval[i] = (double)tval[i];
          Util::MemoryPool<Mem::Main>::instance()->release_memory(tval);

          uint64_t trows(rows());
          uint64_t tcolumns(columns());
          uint64_t elements(this->_indices_size.at(0));
          file.write((const char *)&trows, sizeof(uint64_t));
          file.write((const char *)&tcolumns, sizeof(uint64_t));
          file.write((const char *)&elements, sizeof(uint64_t));
          file.write((const char *)ccol_ind, (long)(elements * sizeof(uint64_t)));
          file.write((const char *)crow_ptr, (long)((trows + 1) * sizeof(uint64_t)));
          file.write((const char *)cval, (long)(elements * sizeof(double)));

          delete[] ccol_ind;
          delete[] crow_ptr;
          delete[] cval;
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
          SparseMatrixCSR<Mem::Main, DT_, IT_> temp;
          temp.convert(*this);

          file << "%%MatrixMarket matrix coordinate real general" << std::endl;
          file << temp.rows() << " " << temp.columns() << " " << temp.used_elements() << std::endl;

          for (Index row(0) ; row < rows() ; ++row)
          {
            const Index end(temp.row_ptr()[row + 1]);
            for (Index i(temp.row_ptr()[row]) ; i < end ; ++i)
            {
              file << row + 1 << " " << temp.col_ind()[i] + 1 << " " << std::scientific << temp.val()[i] << std::endl;
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
        DT_ operator()(Index row, Index col) const
        {
          CONTEXT("When retrieving SparseMatrixCSR element");

          ASSERT(row < rows(), "Error: " + stringify(row) + " exceeds sparse matrix csr row size " + stringify(rows()) + " !");
          ASSERT(col < columns(), "Error: " + stringify(col) + " exceeds sparse matrix csr column size " + stringify(columns()) + " !");

          for (Index i(Index(Util::MemoryPool<Mem_>::get_element(this->_indices.at(1), row))) ; i < Index(Util::MemoryPool<Mem_>::get_element(this->_indices.at(1), row + 1)) ; ++i)
          {
            if (Index(Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i)) == col)
              return Util::MemoryPool<Mem_>::get_element(this->_elements.at(0), i);
            if (Index(Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i)) > col)
              return zero_element();
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
        const Index & rows() const
        {
          return this->_scalar_index.at(1);
        }

        /**
         * \brief Retrieve matrix column count.
         *
         * \returns Matrix column count.
         */
        const Index & columns() const
        {
          return this->_scalar_index.at(2);
        }

        /**
         * \brief Retrieve non zero element count.
         *
         * \returns Non zero element count.
         */
        const Index & used_elements() const override
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
          return this->_indices.at(0);
        }

        IT_ const * col_ind() const
        {
          return this->_indices.at(0);
        }

        /**
         * \brief Retrieve non zero element array.
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
         * \brief Retrieve row start index array.
         *
         * \returns Row start index array.
         */
        IT_ * row_ptr()
        {
          return this->_indices.at(1);
        }

        IT_ const * row_ptr() const
        {
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
          return "SparseMatrixCSR";
        }

        /**
         * \brief Performs \f$this \leftarrow x\f$.
         *
         * \param[in] x The Matrix to be copied.
         */
        void copy(const SparseMatrixCSR & x)
        {
          this->_copy_content(x);
        }

        /**
         * \brief Performs \f$this \leftarrow x\f$.
         *
         * \param[in] x The Matrix to be copied.
         */
        template <typename Mem2_>
        void copy(const SparseMatrixCSR<Mem2_, DT_, IT_> & x)
        {
          this->_copy_content(x);
        }

        ///@name Linear algebra operations
        ///@{
        /**
         * \brief Calculate \f$this \leftarrow y + \alpha x\f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
         *
         * \param[in] x The first summand matrix to be scaled.
         * \param[in] y The second summand matrix
         * \param[in] alpha A scalar to multiply x with.
         */
        template <typename Algo_>
        void axpy(
          const SparseMatrixCSR & x,
          const SparseMatrixCSR & y,
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

          // check for special cases
          // r <- x + y
          if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
            Arch::Sum<Mem_, Algo_>::value(this->val(), x.val(), y.val(), this->used_elements());
          // r <- y - x
          else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
            Arch::Difference<Mem_, Algo_>::value(this->val(), y.val(), x.val(), this->used_elements());
          // r <- y
          else if(Math::abs(alpha) < Math::eps<DT_>())
            this->copy(y);
          // r <- y + alpha*x
          else
            Arch::Axpy<Mem_, Algo_>::dv(this->val(), alpha, x.val(), y.val(), this->used_elements());
        }

        /**
         * \brief Calculate \f$this \leftarrow \alpha x \f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
         *
         * \param[in] x The matrix to be scaled.
         * \param[in] alpha A scalar to scale x with.
         */
        template <typename Algo_>
        void scale(const SparseMatrixCSR & x, const DT_ alpha)
        {
          if (x.rows() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
          if (x.columns() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");

          Arch::Scale<Mem_, Algo_>::value(this->val(), x.val(), alpha, this->used_elements());
        }

        /**
         * \brief Calculate \f$this^\top \f$
         *
         * \return The transposed matrix
         */
        SparseMatrixCSR transpose()
        {
          SparseMatrixCSR x_t;
          x_t.transpose(*this);
          return x_t;
        }

        /**
         * \brief Calculate \f$this \leftarrow x^\top \f$
         *
         * \param[in] x The matrix to be transposed.
         */
        void transpose(const SparseMatrixCSR & x)
        {
          SparseMatrixCSR<Mem::Main, DT_, IT_> tx;
          tx.convert(x);

          const Index txrows(tx.rows());
          const Index txcolumns(tx.columns());
          const Index txused_elements(tx.used_elements());

          const IT_ * ptxcol_ind(tx.col_ind());
          const IT_ * ptxrow_ptr(tx.row_ptr());
          const DT_ * ptxval(tx.val());

          DenseVector<Mem::Main, IT_, IT_> tcol_ind(txused_elements);
          DenseVector<Mem::Main, DT_, IT_> tval(txused_elements);
          DenseVector<Mem::Main, IT_, IT_> trow_ptr(txcolumns + 1, IT_(0));

          IT_ * ptcol_ind(tcol_ind.elements());
          DT_ * ptval(tval.elements());
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
            for (Index k(ptxrow_ptr[i]); k < ptxrow_ptr[i+1]; ++k)
            {
              const Index l(ptxcol_ind[k]);
              const Index j(ptrow_ptr[l]);
              ptval[j] = ptxval[k];
              ptcol_ind[j] = IT_(i);
              ++ptrow_ptr[l];
            }
          }

          for (Index i(txcolumns); i > 0; --i)
          {
            ptrow_ptr[i] = ptrow_ptr[i - 1];
          }
          ptrow_ptr[0] = 0;

          SparseMatrixCSR<Mem::Main, DT_, IT_> tx_t(txcolumns, txrows, tcol_ind, tval, trow_ptr);

          SparseMatrixCSR<Mem_, DT_, IT_> x_t;
          x_t.convert(tx_t);
          this->assign(x_t);
        }

        /**
         * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_i\f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
         *
         * \param[in] x The matrix whose rows are to be scaled.
         * \param[in] s The vector to the scale the rows by.
         */
        template<typename Algo_>
        void scale_rows(const SparseMatrixCSR & x, const DenseVector<Mem_,DT_,IT_> & s)
        {
          if (x.rows() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
          if (x.columns() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
          if (s.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

          Arch::ScaleRows<Mem_, Algo_>::csr(this->val(), x.val(), this->col_ind(), this->row_ptr(),
            s.elements(), this->rows(), this->columns(), this->used_elements());
        }

        /**
         * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_j\f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
         *
         * \param[in] x The matrix whose columns are to be scaled.
         * \param[in] s The vector to the scale the columns by.
         */
        template<typename Algo_>
        void scale_cols(const SparseMatrixCSR & x, const DenseVector<Mem_,DT_,IT_> & s)
        {
          if (x.rows() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
          if (x.columns() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
          if (s.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

          Arch::ScaleCols<Mem_, Algo_>::csr(this->val(), x.val(), this->col_ind(), this->row_ptr(),
            s.elements(), this->rows(), this->columns(), this->used_elements());
        }

        /**
         * \brief Calculate \f$ r \leftarrow this\cdot x \f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
         *
         * \param[out] r The vector that recieves the result.
         * \param[in] x The vector to be multiplied by this matrix.
         */
        template<typename Algo_>
        void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
        {
          if (r.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
          if (x.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

          Arch::ProductMatVec<Mem_, Algo_>::csr(r.elements(), this->val(), this->col_ind(), this->row_ptr(),
            x.elements(), this->rows(), columns(), used_elements());
        }

        /**
         * \brief Calculate \f$ r \leftarrow y + \alpha this\cdot x \f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
         *
         * \param[out] r The vector that recieves the result.
         * \param[in] x The vector to be multiplied by this matrix.
         * \param[in] y The summand vector.
         * \param[in] alpha A scalar to scale the product with.
         */
        template<typename Algo_>
        void apply(
          DenseVector<Mem_,DT_, IT_> & r,
          const DenseVector<Mem_, DT_, IT_> & x,
          const DenseVector<Mem_, DT_, IT_> & y,
          const DT_ alpha = DT_(1)) const
        {
          if (r.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
          if (x.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
          if (y.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

          // check for special cases
          // r <- y - A*x
          if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          {
            Arch::Defect<Mem_, Algo_>::csr(r.elements(), y.elements(), this->val(), this->col_ind(),
              this->row_ptr(), x.elements(), this->rows(), this->columns(), this->used_elements());
          }
          //r <- y
          else if(Math::abs(alpha) < Math::eps<DT_>())
            r.copy(y);
          // r <- y + alpha*x
          else
          {
            Arch::Axpy<Mem_, Algo_>::csr(r.elements(), alpha, x.elements(), y.elements(),
              this->val(), this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
          }
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
          return prow_ptr[row + 1] - prow_ptr[row];
        }

        /// \cond internal
        /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
        void set_line(const Index row, DT_ * const pval_set, IT_ * const pcol_set,
                             const Index col_start, const Index stride = 1) const
        {
          const IT_ * prow_ptr(this->row_ptr());
          const IT_ * pcol_ind(this->col_ind());
          const DT_ * pval(this->val());

          const Index start(prow_ptr[row]);
          for (Index i(0); i < prow_ptr[row + 1] - prow_ptr[row]; ++i)
          {
            pval_set[i * stride] = pval[start + i];
            pcol_set[i * stride] = pcol_ind[start + i] + IT_(col_start);
          }
        }
        /// \endcond

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
          ASSERT(domain_node < rows(), "Domain node index out of range");
          return &this->_indices.at(0)[this->_indices.at(1)[domain_node]];
        }

        /** \copydoc Adjactor::image_end() */
        inline ImageIterator image_end(Index domain_node) const
        {
          CONTEXT("Graph::image_end()");
          ASSERT(domain_node < rows(), "Domain node index out of range");
          return &this->_indices.at(0)[this->_indices.at(1)[domain_node + 1]];
        }
    };

    /**
     * \brief SparseMatrixCSR comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Mem_, typename Mem2_, typename DT_, typename IT_> bool operator== (const SparseMatrixCSR<Mem_, DT_, IT_> & a, const SparseMatrixCSR<Mem2_, DT_, IT_> & b)
    {
      CONTEXT("When comparing SparseMatrixCSRs");

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
        Util::MemoryPool<Mem_>::instance()->template download<IT_>(col_ind_a, a.col_ind(), a.used_elements());
        val_a = new DT_[a.used_elements()];
        Util::MemoryPool<Mem_>::instance()->template download<DT_>(val_a, a.val(), a.used_elements());
        row_ptr_a = new IT_[a.rows() + 1];
        Util::MemoryPool<Mem_>::instance()->template download<IT_>(row_ptr_a, a.row_ptr(), a.rows() + 1);
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
        Util::MemoryPool<Mem2_>::instance()->template download<IT_>(col_ind_b, b.col_ind(), b.used_elements());
        val_b = new DT_[b.used_elements()];
        Util::MemoryPool<Mem2_>::instance()->template download<DT_>(val_b, b.val(), b.used_elements());
        row_ptr_b = new IT_[b.rows() + 1];
        Util::MemoryPool<Mem2_>::instance()->template download<IT_>(row_ptr_b, b.row_ptr(), b.rows() + 1);
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
     * \brief SparseMatrixCSR streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Mem_, typename DT_, typename IT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixCSR<Mem_, DT_, IT_> & b)
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

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP
