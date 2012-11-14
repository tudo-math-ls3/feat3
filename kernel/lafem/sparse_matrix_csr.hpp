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

      public:
        /// Our datatype
        typedef DT_ data_type;
        /// Our memory architecture type
        typedef Arch_ mem_type;

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
          _used_elements(0)
        {
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
          _zero_element(DT_(0)),
          _used_elements(other.used_elements())
        {
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
          for (typename std::map<unsigned long, DT_>::const_iterator it(other.elements().begin()) ; it != other.elements().end() ; ++it)
          {
            Index row(it->first / _columns);
            Index column(it->first % _columns);

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
            _val[ait] = it->second;
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
          SparseMatrixCSR<Mem::Main, DT_> tother(other);

          this->_size = tother.size();
          this->_rows = tother.rows();
          this->_columns = tother.columns();
          this->_used_elements = tother.used_elements();
          this->_zero_element = tother.zero_element();

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
          MemoryPool<Arch2_>::download(temp, tother.get_elements().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_elements().at(0), temp, dest_size);
          ::free(temp);

          src_size = (tother.get_indices_size().at(0) * sizeof(Index));
          dest_size = (tother.get_indices_size().at(0) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, tother.get_indices().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(0), temp, dest_size);
          ::free(temp);

          src_size = (tother.get_indices_size().at(1) * sizeof(Index));
          dest_size = (tother.get_indices_size().at(1) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, tother.get_indices().at(1), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(1), temp, dest_size);
          ::free(temp);

          src_size = (tother.get_indices_size().at(2) * sizeof(Index));
          dest_size = (tother.get_indices_size().at(2) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, tother.get_indices().at(2), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(2), temp, dest_size);
          ::free(temp);
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
         */
        explicit SparseMatrixCSR(Index rows, Index columns, const DenseVector<Arch_, Index> & col_ind, const DenseVector<Arch_, DT_> & val, const DenseVector<Arch_, Index> & row_ptr, const DenseVector<Arch_, Index> & row_ptr_end) :
          Container<Arch_, DT_>(rows * columns),
          _rows(rows),
          _columns(columns),
          _zero_element(DT_(0)),
          // \TODO use real element count - be aware of non used elements
          _used_elements(val.size())
        {
          DenseVector<Arch_, Index> tcol_ind(col_ind.size());
          copy(tcol_ind, col_ind);
          DenseVector<Arch_, DT_> tval(val.size());
          copy(tval, val);
          DenseVector<Arch_, Index> trow_ptr(row_ptr.size());
          copy(trow_ptr, row_ptr);
          DenseVector<Arch_, Index> trow_ptr_end(row_ptr_end.size());
          copy(trow_ptr_end, row_ptr_end);

          this->_elements.push_back(tval.get_elements().at(0));
          this->_elements_size.push_back(tval.size());
          this->_indices.push_back(tcol_ind.get_elements().at(0));
          this->_indices_size.push_back(tcol_ind.size());
          this->_indices.push_back(trow_ptr.get_elements().at(0));
          this->_indices_size.push_back(trow_ptr.size());
          this->_indices.push_back(trow_ptr_end.get_elements().at(0));
          this->_indices_size.push_back(trow_ptr_end.size());

          this->_val = this->_elements.at(0);
          this->_col_ind = this->_indices.at(0);
          this->_row_ptr = this->_indices.at(1);
          this->_row_ptr_end = this->_indices.at(2);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));
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
          _used_elements(other._used_elements)
        {
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
          _used_elements(other.used_elements())
        {
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
          if (this == &other)
            return *this;

          this->_size = other.size();
          this->_rows = other.rows();
          this->_columns = other.columns();
          this->_used_elements = other.used_elements();
          this->_zero_element = other._zero_element;

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
          this->_size = other.size();
          this->_rows = other.rows();
          this->_columns = other.columns();
          this->_used_elements = other.used_elements();
          this->_zero_element = other.zero_element();

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
         * \brief Retrieve specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         *
         * \returns Specific matrix element.
         */
        DT_ operator()(Index row, Index col) const
        {
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
    };

    /**
     * \brief SparseMatrixCSR comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Arch_, typename Arch2_, typename DT_> bool operator== (const SparseMatrixCSR<Arch_, DT_> & a, const SparseMatrixCSR<Arch2_, DT_> & b)
    {
      if (a.rows() != b.rows())
        return false;
      if (a.columns() != b.columns())
        return false;
      if (a.used_elements() != b.used_elements())
        return false;
      if (a.zero_element() != b.zero_element())
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
