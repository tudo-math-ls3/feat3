#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/algorithm.hpp>
#include <kernel/util/graph.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief CSR based sparse matrix.
     *
     * \tparam Arch_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse row format.
     *
     * \author Dirk Ribbrock
     */
    template <typename Arch_, typename DT_>
    class SparseMatrixCSR : public Container<Arch_, DT_>
    {
      public:
        /// ImageIterator typedef for Adjactor interface implementation
        typedef const Index* ImageIterator;

      private:
        /// Column indices.
        Index * _col_ind;
        /// Non zero values.
        DT_ * _val;
        /// Row start indices (including matrix end index).
        Index * _row_ptr;
        /// Row end indices.
        Index * _row_ptr_end;
        /// Row count.
        Index _rows;
        /// Column count.
        Index _columns;
        /// Our non zero element.
        DT_ _zero_element;
        /// Our non zero element count.
        Index _used_elements;
        /// Our layout related hash
        unsigned long _hash;

        void _read_from_csr(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          _read_from_csr(file);
          file.close();
        }

        void _read_from_csr(std::istream& file)
        {
          uint64_t rows;
          uint64_t columns;
          uint64_t elements;
          file.read((char *)&rows, sizeof(uint64_t));
          file.read((char *)&columns, sizeof(uint64_t));
          file.read((char *)&elements, sizeof(uint64_t));

          this->_size = Index(rows * columns);
          _rows = Index(rows);
          _columns = Index(columns);
          _zero_element = DT_(0);
          _used_elements = elements;

          uint64_t * ccol_ind = new uint64_t[elements];
          file.read((char *)ccol_ind, elements * sizeof(uint64_t));
          _col_ind = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((elements) * sizeof(Index));
          for (Index i(0) ; i < elements ; ++i)
            _col_ind[i] = Index(ccol_ind[i]);
          delete[] ccol_ind;

          uint64_t * crow_ptr = new uint64_t[rows + 1];
          file.read((char *)crow_ptr, (rows + 1) * sizeof(uint64_t));
          _row_ptr = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((rows + 1) * sizeof(Index));
          for (Index i(0) ; i < rows + 1 ; ++i)
            _row_ptr[i] = Index(crow_ptr[i]);
          delete[] crow_ptr;

          uint64_t * crow_ptr_end = new uint64_t[rows];
          file.read((char *)crow_ptr_end, (rows) * sizeof(uint64_t));
          _row_ptr_end = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((rows) * sizeof(Index));
          for (Index i(0) ; i < rows ; ++i)
            _row_ptr_end[i] = Index(crow_ptr_end[i]);
          delete[] crow_ptr_end;

          double * cval = new double[elements];
          file.read((char *)cval, elements * sizeof(double));
          _val = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((elements) * sizeof(DT_));
          for (Index i(0) ; i < elements ; ++i)
            _val[i] = DT_(cval[i]);
          delete[] cval;

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(elements * sizeof(DT_)));
          this->_elements_size.push_back(elements);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(elements * sizeof(Index)));
          this->_indices_size.push_back(elements);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory((_rows + 1)* sizeof(Index)));
          this->_indices_size.push_back(_rows + 1);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory((_rows)* sizeof(Index)));
          this->_indices_size.push_back(_rows);

          MemoryPool<Arch_>::upload(this->get_elements().at(0), _val, elements * sizeof(DT_));
          MemoryPool<Arch_>::upload(this->get_indices().at(0), _col_ind, elements * sizeof(Index));
          MemoryPool<Arch_>::upload(this->get_indices().at(1), _row_ptr, (_rows  + 1) * sizeof(Index));
          MemoryPool<Arch_>::upload(this->get_indices().at(2), _row_ptr_end, (_rows) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(_val);
          MemoryPool<Mem::Main>::instance()->release_memory(_col_ind);
          MemoryPool<Mem::Main>::instance()->release_memory(_row_ptr);
          MemoryPool<Mem::Main>::instance()->release_memory(_row_ptr_end);

          _col_ind = this->_indices.at(0);
          _val = this->_elements.at(0);
          _row_ptr = this->_indices.at(1);
          _row_ptr_end = this->_indices.at(2);

          _hash = MemoryPool<Arch_>::instance()->generate_hash(_row_ptr, (this->_rows + 1) * sizeof(Index));
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Arch_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit SparseMatrixCSR() :
          Container<Arch_, DT_> (0),
          _rows(0),
          _columns(0),
          _zero_element(DT_(0)),
          _used_elements(0),
          _hash(0)
        {
          CONTEXT("When creating SparseMatrixCSR");
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix in COO format.
         *
         * Creates a CSR matrix based on the COO source matrix.
         */
        explicit SparseMatrixCSR(const SparseMatrixCOO<Arch_, DT_> & other) :
          Container<Arch_, DT_>(other.size()),
          _rows(other.rows()),
          _columns(other.columns()),
          _zero_element(other.zero_element()),
          _used_elements(other.used_elements())
        {
          CONTEXT("When creating SparseMatrixCSR");

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(_used_elements * sizeof(DT_)));
          this->_elements_size.push_back(_used_elements);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_used_elements * sizeof(Index)));
          this->_indices_size.push_back(_used_elements);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory((_rows + 1) * sizeof(Index)));
          this->_indices_size.push_back(_rows + 1);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory((_rows) * sizeof(Index)));
          this->_indices_size.push_back(_rows);

          _col_ind = this->_indices.at(0);
          _val = this->_elements.at(0);
          _row_ptr = this->_indices.at(1);
          _row_ptr_end = this->_indices.at(2);

          Index ait(0);
          Index current_row(0);
          _row_ptr[current_row] = 0;
          for (Index it(0) ; it < other.used_elements() ; ++it)
          {
            Index row(other.row()[it]);
            Index column(other.column()[it]);

            if (current_row < row)
            {
              _row_ptr_end[current_row] = ait;
              for (unsigned long i(current_row + 1) ; i < row ; ++i)
              {
                _row_ptr[i] = ait;
                _row_ptr_end[i] = ait;
              }
              current_row = row;
              _row_ptr[current_row] = ait;
            }
            _val[ait] = other.val()[it];
            _col_ind[ait] = column;
            ++ait;
          }
          _row_ptr_end[current_row] = ait;
          for (unsigned long i(current_row + 1) ; i < _rows ; ++i)
          {
            _row_ptr[i] = ait;
            _row_ptr_end[i] = ait;
          }
          _row_ptr[_rows] = ait;

          _hash = MemoryPool<Arch_>::instance()->generate_hash(_row_ptr, (this->_rows + 1) * sizeof(Index));
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix in COO format.
         *
         * Creates a CSR matrix based on the COO source matrix from another memory architecture.
         */
        template <typename Arch2_>
        explicit SparseMatrixCSR(const SparseMatrixCOO<Arch2_, DT_> & other) :
          Container<Arch_, DT_>(other.size()),
          _rows(other.rows()),
          _columns(other.columns()),
          _zero_element(other.zero_element()),
          _used_elements(other.used_elements())
        {
          CONTEXT("When creating SparseMatrixCSR");

          SparseMatrixCOO<Mem::Main, DT_> xother(other);
          SparseMatrixCSR<Mem::Main, DT_> tother(xother);

          this->_size = tother.size();
          this->_rows = tother.rows();
          this->_columns = tother.columns();
          this->_used_elements = tother.used_elements();
          this->_zero_element = tother.zero_element();
          this->_hash = tother.hash();

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(tother.used_elements() * sizeof(DT_)));
          this->_elements_size.push_back(this->_used_elements);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(tother.used_elements() * sizeof(Index)));
          this->_indices_size.push_back(this->_used_elements);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory((_rows + 1) * sizeof(Index)));
          this->_indices_size.push_back(_rows + 1);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_rows * sizeof(Index)));
          this->_indices_size.push_back(_rows);

          _col_ind = this->_indices.at(0);
          _val = this->_elements.at(0);
          _row_ptr = this->_indices.at(1);
          _row_ptr_end = this->_indices.at(2);

          Index src_size(tother.get_elements_size().at(0) * sizeof(DT_));
          Index dest_size(tother.get_elements_size().at(0) * sizeof(DT_));
          void * temp(::malloc(src_size));
          MemoryPool<Mem::Main>::download(temp, tother.get_elements().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_elements().at(0), temp, dest_size);
          ::free(temp);

          src_size = (tother.get_indices_size().at(0) * sizeof(Index));
          dest_size = (tother.get_indices_size().at(0) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Mem::Main>::download(temp, tother.get_indices().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(0), temp, dest_size);
          ::free(temp);

          src_size = (tother.get_indices_size().at(1) * sizeof(Index));
          dest_size = (tother.get_indices_size().at(1) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Mem::Main>::download(temp, tother.get_indices().at(1), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(1), temp, dest_size);
          ::free(temp);

          src_size = (tother.get_indices_size().at(2) * sizeof(Index));
          dest_size = (tother.get_indices_size().at(2) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Mem::Main>::download(temp, tother.get_indices().at(2), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(2), temp, dest_size);
          ::free(temp);
        }

        /**
         * \brief Constructor
         *
         * \param[in] filename The source file in binary CSR format.
         *
         * Creates a CSR matrix based on the source file.
         */
        explicit SparseMatrixCSR(String filename) :
          Container<Arch_, DT_>(0)
        {
          CONTEXT("When creating SparseMatrixCSR");

          _read_from_csr(filename);
        }

        /**
         * \brief Constructor
         *
         * \param[in] file The source file in binary CSR format.
         *
         * Creates a CSR matrix based on the source file.
         */
        explicit SparseMatrixCSR(std::istream& file) :
          Container<Arch_, DT_>(0)
        {
          CONTEXT("When creating SparseMatrixCSR");

          _read_from_csr(file);
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
         * \param[in] row_ptr_end Vector with end indices of all rows.
         *
         * Creates a matrix with given dimensions and content.
         * During creation, the input data are copied. Thus the matrix is independent of later input vector modifications.
         */
        explicit SparseMatrixCSR(Index rows, Index columns, const DenseVector<Arch_, Index> & col_ind, const DenseVector<Arch_, DT_> & val, const DenseVector<Arch_, Index> & row_ptr, const DenseVector<Arch_, Index> & row_ptr_end) :
          Container<Arch_, DT_>(rows * columns),
          _rows(rows),
          _columns(columns),
          _zero_element(DT_(0)),
          // \TODO use real element count - be aware of non used elements
          _used_elements(val.size())
        {
          CONTEXT("When creating SparseMatrixCSR");

          this->_elements.push_back(val.get_elements().at(0));
          this->_elements_size.push_back(val.size());
          this->_indices.push_back(col_ind.get_elements().at(0));
          this->_indices_size.push_back(col_ind.size());
          this->_indices.push_back(row_ptr.get_elements().at(0));
          this->_indices_size.push_back(row_ptr.size());
          this->_indices.push_back(row_ptr_end.get_elements().at(0));
          this->_indices_size.push_back(row_ptr_end.size());

          this->_val = this->_elements.at(0);
          this->_col_ind = this->_indices.at(0);
          this->_row_ptr = this->_indices.at(1);
          this->_row_ptr_end = this->_indices.at(2);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

          _hash = MemoryPool<Arch_>::instance()->generate_hash(_row_ptr, (this->_rows + 1) * sizeof(Index));
        }

        /**
         * \brief Constructor
         *
         * \param[in] graph The Graph, the matrix will be created from.
         *
         * Creates a matrix from a given graph.
         */
        explicit SparseMatrixCSR(const Graph & graph) :
          Container<Arch_, DT_>(graph.get_num_nodes_domain() * graph.get_num_nodes_image()),
          _zero_element(DT_(0)),
          _used_elements(graph.get_num_indices())
        {
          CONTEXT("When creating SparseMatrixCSR");

          this->_rows = graph.get_num_nodes_domain();
          this->_columns = graph.get_num_nodes_image();

          const Index* dom_ptr = graph.get_domain_ptr();
          const Index* dom_end = graph.get_domain_end();
          const Index* img_idx = graph.get_image_idx();
          if(dom_end == nullptr)
          {
            dom_end = &dom_ptr[1];
          }

          DenseVector<Arch_, Index> col_ind(_used_elements);
          DenseVector<Arch_, DT_> val(_used_elements);
          DenseVector<Arch_, Index> row_ptr(_rows + 1);
          DenseVector<Arch_, Index> row_ptr_end(_rows);

          Index* prow_ptr = row_ptr.elements();
          for(Index i(0); i <= _rows; ++i)
          {
            prow_ptr[i] = dom_ptr[i];
          }

          Index* prow_end = row_ptr_end.elements();
          for(Index i(0); i < _rows; ++i)
          {
            prow_end[i] = dom_end[i];
          }

          Index* pcol_ind = col_ind.elements();
          for(Index i(0); i < _used_elements; ++i)
          {
            pcol_ind[i] = img_idx[i];
          }

          this->_elements.push_back(val.get_elements().at(0));
          this->_elements_size.push_back(val.size());
          this->_indices.push_back(col_ind.get_elements().at(0));
          this->_indices_size.push_back(col_ind.size());
          this->_indices.push_back(row_ptr.get_elements().at(0));
          this->_indices_size.push_back(row_ptr.size());
          this->_indices.push_back(row_ptr_end.get_elements().at(0));
          this->_indices_size.push_back(row_ptr_end.size());

          this->_val = this->_elements.at(0);
          this->_col_ind = this->_indices.at(0);
          this->_row_ptr = this->_indices.at(1);
          this->_row_ptr_end = this->_indices.at(2);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

          _hash = MemoryPool<Arch_>::instance()->generate_hash(_row_ptr, (this->_rows + 1) * sizeof(Index));
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a shallow copy of a given matrix.
         */
        SparseMatrixCSR(const SparseMatrixCSR<Arch_, DT_> & other) :
          Container<Arch_, DT_>(other),
          _rows(other._rows),
          _columns(other._columns),
          _zero_element(other._zero_element),
          _used_elements(other._used_elements),
          _hash(other.hash())
        {
          CONTEXT("When copying SparseMatrixCSR");

          this->_val = this->_elements.at(0);
          this->_col_ind = this->_indices.at(0);
          this->_row_ptr = this->_indices.at(1);
          this->_row_ptr_end = this->_indices.at(2);
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a copy of a given matrix from another memory architecture.
         */
        template <typename Arch2_, typename DT2_>
        SparseMatrixCSR(const SparseMatrixCSR<Arch2_, DT2_> & other) :
          Container<Arch_, DT_>(other),
          _rows(other.rows()),
          _columns(other.columns()),
          _zero_element(other.zero_element()),
          _used_elements(other.used_elements()),
          _hash(other.hash())
        {
          CONTEXT("When copying SparseMatrixCSR");

          this->_val = this->_elements.at(0);
          this->_col_ind = this->_indices.at(0);
          this->_row_ptr = this->_indices.at(1);
          this->_row_ptr_end = this->_indices.at(2);
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        SparseMatrixCSR<Arch_, DT_> clone()
        {
          CONTEXT("When cloning SparseMatrixCSR");

          DenseVector<Arch_, Index> col_ind(this->_used_elements, this->_col_ind);
          DenseVector<Arch_, DT_> val(this->_used_elements, this->_val);
          DenseVector<Arch_, Index> row_ptr(this->_rows + 1, this->_row_ptr);
          DenseVector<Arch_, Index> row_ptr_end(this->_rows, this->_row_ptr_end);

          SparseMatrixCSR<Arch_, DT_> t(this->_rows, this->_columns, col_ind.clone(), val.clone(),
              row_ptr.clone(), row_ptr_end.clone());

          return t;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns another matrix to the target matrix.
         */
        SparseMatrixCSR<Arch_, DT_> & operator= (const SparseMatrixCSR<Arch_, DT_> & other)
        {
          CONTEXT("When assigning SparseMatrixCSR");

          if (this == &other)
            return *this;

          this->_size = other.size();
          this->_rows = other.rows();
          this->_columns = other.columns();
          this->_used_elements = other.used_elements();
          this->_zero_element = other._zero_element;
          this->_hash = other.hash();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();

          std::vector<DT_ *> new_elements = other.get_elements();
          std::vector<Index *> new_indices = other.get_indices();

          this->_elements.assign(new_elements.begin(), new_elements.end());
          this->_indices.assign(new_indices.begin(), new_indices.end());
          this->_elements_size.assign(other.get_elements_size().begin(), other.get_elements_size().end());
          this->_indices_size.assign(other.get_indices_size().begin(), other.get_indices_size().end());

          _col_ind = this->_indices.at(0);
          _val = this->_elements.at(0);
          _row_ptr = this->_indices.at(1);
          _row_ptr_end = this->_indices.at(2);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

          return *this;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns a matrix from another memory architecture to the target matrix.
         */
        template <typename Arch2_, typename DT2_>
        SparseMatrixCSR<Arch_, DT_> & operator= (const SparseMatrixCSR<Arch2_, DT2_> & other)
        {
          CONTEXT("When assigning SparseMatrixCSR");

          this->_size = other.size();
          this->_rows = other.rows();
          this->_columns = other.columns();
          this->_used_elements = other.used_elements();
          this->_zero_element = other.zero_element();
          this->_hash = other.hash();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(other.used_elements() * sizeof(DT_)));
          this->_elements_size.push_back(this->_used_elements);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(other.used_elements() * sizeof(Index)));
          this->_indices_size.push_back(this->_used_elements);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory((_rows + 1) * sizeof(Index)));
          this->_indices_size.push_back(_rows + 1);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_rows * sizeof(Index)));
          this->_indices_size.push_back(_rows);

          _col_ind = this->_indices.at(0);
          _val = this->_elements.at(0);
          _row_ptr = this->_indices.at(1);
          _row_ptr_end = this->_indices.at(2);

          Index src_size(other.get_elements_size().at(0) * sizeof(DT2_));
          Index dest_size(other.get_elements_size().at(0) * sizeof(DT_));
          void * temp(::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_elements().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_elements().at(0), temp, dest_size);
          ::free(temp);

          src_size = (other.get_indices_size().at(0) * sizeof(Index));
          dest_size = (other.get_indices_size().at(0) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_indices().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(0), temp, dest_size);
          ::free(temp);

          src_size = (other.get_indices_size().at(1) * sizeof(Index));
          dest_size = (other.get_indices_size().at(1) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_indices().at(1), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(1), temp, dest_size);
          ::free(temp);

          src_size = (other.get_indices_size().at(2) * sizeof(Index));
          dest_size = (other.get_indices_size().at(2) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_indices().at(2), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(2), temp, dest_size);
          ::free(temp);

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
            case fm_csr:
              write_out_csr(filename);
              break;
            case fm_m:
              write_out_m(filename);
              break;
            default:
                throw InternalError("Filemode not supported!");
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
            case fm_csr:
              write_out_csr(file);
              break;
            case fm_m:
              write_out_m(file);
              break;
            default:
                throw InternalError("Filemode not supported!");
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
            throw InternalError("Unable to open Matrix file " + filename);
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
          if (typeid(DT_) != typeid(double))
            std::cout<<"Warning: You are writing out an csr matrix with less than double precission!"<<std::endl;

          Index * col_ind = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory(this->_indices_size.at(0) * sizeof(Index));
          MemoryPool<Arch_>::download(col_ind, _col_ind, this->_indices_size.at(0) * sizeof(Index));
          uint64_t * ccol_ind = new uint64_t[this->_indices_size.at(0)];
          for (Index i(0) ; i < this->_indices_size.at(0) ; ++i)
            ccol_ind[i] = col_ind[i];
          MemoryPool<Mem::Main>::instance()->release_memory(col_ind);

          Index * row_ptr = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory(this->_indices_size.at(1) * sizeof(Index));
          MemoryPool<Arch_>::download(row_ptr, _row_ptr, this->_indices_size.at(1) * sizeof(Index));
          uint64_t * crow_ptr = new uint64_t[this->_indices_size.at(1)];
          for (Index i(0) ; i < this->_indices_size.at(1) ; ++i)
            crow_ptr[i] = row_ptr[i];
          MemoryPool<Mem::Main>::instance()->release_memory(row_ptr);

          Index * row_ptr_end = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory(this->_indices_size.at(2) * sizeof(Index));
          MemoryPool<Arch_>::download(row_ptr_end, _row_ptr_end, this->_indices_size.at(2) * sizeof(Index));
          uint64_t * crow_ptr_end = new uint64_t[this->_indices_size.at(2)];
          for (Index i(0) ; i < this->_indices_size.at(2) ; ++i)
            crow_ptr_end[i] = row_ptr_end[i];
          MemoryPool<Mem::Main>::instance()->release_memory(row_ptr_end);

          DT_ * val = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory(this->_elements_size.at(0) * sizeof(DT_));
          MemoryPool<Arch_>::download(val, _val, this->_elements_size.at(0) * sizeof(DT_));
          double * cval = new double[this->_elements_size.at(0)];
          for (Index i(0) ; i < this->_elements_size.at(0) ; ++i)
            cval[i] = val[i];
          MemoryPool<Mem::Main>::instance()->release_memory(val);

          uint64_t rows(_rows);
          uint64_t columns(_columns);
          uint64_t elements(this->_indices_size.at(0));
          file.write((const char *)&rows, sizeof(uint64_t));
          file.write((const char *)&columns, sizeof(uint64_t));
          file.write((const char *)&elements, sizeof(uint64_t));
          file.write((const char *)ccol_ind, elements * sizeof(uint64_t));
          file.write((const char *)crow_ptr, (rows + 1) * sizeof(uint64_t));
          file.write((const char *)crow_ptr_end, rows * sizeof(uint64_t));
          file.write((const char *)cval, elements * sizeof(double));

          delete[] ccol_ind;
          delete[] crow_ptr;
          delete[] crow_ptr_end;
          delete[] cval;
        }

        /**
         * \brief Write out matrix to matlab m file.
         *
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out_m(String filename) const
        {
          std::ofstream file(filename.c_str(), std::ofstream::out);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          write_out_m(file);
          file.close();
        }

        /**
         * \brief Write out matrix to matlab m file.
         *
         * \param[in] file The stream that shall be written to.
         */
        void write_out_m(std::ostream& file) const
        {
          SparseMatrixCSR<Mem::Main, DT_> temp(*this);

          file << "data = [" << std::endl;
          for (Index row(0) ; row < _rows ; ++row)
          {
            const Index end(temp.row_ptr_end()[row]);
            for (Index i(temp.row_ptr()[row]) ; i < end ; ++i)
            {
              if (temp.val()[i] != DT_(0))
              {
                file << row + 1 << " " << temp.col_ind()[i] + 1 << " " << std::scientific << (double)temp.val()[i] << ";" << std::endl;
              }
            }
          }
          file << "];" << std::endl;
          file << "mat=sparse(data(:,1),data(:,2),data(:,3));";
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

          ASSERT(row < this->_rows, "Error: " + stringify(row) + " exceeds sparse matrix csr row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + " exceeds sparse matrix csr column size " + stringify(this->_columns) + " !");

          if (typeid(Arch_) == typeid(Mem::Main))
          {
            for (unsigned long i(_row_ptr[row]) ; i < _row_ptr_end[row] ; ++i)
            {
              if (_col_ind[i] == col)
                return _val[i];
              if (_col_ind[i] > col)
                return _zero_element;
            }
            return _zero_element;
          }
          else
          {
            SparseMatrixCSR<Mem::Main, DT_> temp(*this);
            return temp(row, col);
          }
        }

        /**
         * \brief Retrieve matrix row count.
         *
         * \returns Matrix row count.
         */
        const Index & rows() const
        {
          return this->_rows;
        }

        /**
         * \brief Retrieve matrix column count.
         *
         * \returns Matrix column count.
         */
        const Index & columns() const
        {
          return this->_columns;
        }

        /**
         * \brief Retrieve non zero element count.
         *
         * \returns Non zero element count.
         */
        const Index & used_elements() const
        {
          return this->_used_elements;
        }

        /**
         * \brief Retrieve column indices array.
         *
         * \returns Column indices array.
         */
        Index * col_ind() const
        {
          return _col_ind;
        }

        /**
         * \brief Retrieve non zero element array.
         *
         * \returns Non zero element array.
         */
        DT_ * val() const
        {
          return _val;
        }

        /**
         * \brief Retrieve row start index array.
         *
         * \returns Row start index array.
         */
        Index * row_ptr() const
        {
          return _row_ptr;
        }

        /**
         * \brief Retrieve row end index array.
         *
         * \returns Row end index array.
         */
        Index * row_ptr_end() const
        {
          return _row_ptr_end;
        }

        /**
         * \brief Retrieve non zero element.
         *
         * \returns Non zero element.
         */
        const DT_ zero_element() const
        {
          return _zero_element;
        }

        /**
         * \brief Retrieve layout hash.
         *
         * \returns Hash value of matrix layout.
         */
        unsigned long hash() const
        {
          return _hash;
        }

        /* ******************************************************************* */
        /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
        /* ******************************************************************* */
      public:
        /** \copydoc Adjactor::get_num_nodes_domain() */
        inline Index get_num_nodes_domain() const
        {
          return _rows;
        }

        /** \copydoc Adjactor::get_num_nodes_image() */
        inline Index get_num_nodes_image() const
        {
          return _columns;
        }

        /** \copydoc Adjactor::image_begin() */
        inline ImageIterator image_begin(Index domain_node) const
        {
          ASSERT(domain_node < _rows, "Domain node index out of range");
          return &_col_ind[_row_ptr[domain_node]];
        }

        /** \copydoc Adjactor::image_end() */
        inline ImageIterator image_end(Index domain_node) const
        {
          CONTEXT("Graph::image_end()");
          ASSERT(domain_node < _rows, "Domain node index out of range");
          return &_col_ind[_row_ptr_end[domain_node]];
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String type_name()
        {
          return "SparseMatrixCSR";
        }
    };

    /**
     * \brief SparseMatrixCSR comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Arch_, typename Arch2_, typename DT_> bool operator== (const SparseMatrixCSR<Arch_, DT_> & a, const SparseMatrixCSR<Arch2_, DT_> & b)
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
      if (a.hash() != b.hash())
        return false;

      for (Index i(0) ; i < a.rows() ; ++i)
      {
        for (Index j(0) ; j < a.columns() ; ++j)
        {
          if (a(i, j) != b(i, j))
            return false;
        }
      }

      return true;
    }

    /**
     * \brief SparseMatrixCSR streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Arch_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixCSR<Arch_, DT_> & b)
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
