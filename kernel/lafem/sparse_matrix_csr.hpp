#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/matrix_base.hpp>
#include <kernel/lafem/algorithm.hpp>
#include <kernel/adjacency/graph.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    //forward declarations
    template <typename Mem_, typename DT_>
    class SparseMatrixELL;

    template <typename Mem_, typename DT_>
    class SparseMatrixCOO;

    /**
     * \brief CSR based sparse matrix.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse row format.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: column index per non zero element \n
     * _indices[1]: row start index (including matrix end index)\n
     * _indices[2]: row end index \n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_index[4]: layout related hash value \n
     * _scalar_dt[0]: zero element
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_>
    class SparseMatrixCSR : public Container<Mem_, DT_>, public MatrixBase
    {
      public:
        /// ImageIterator typedef for Adjactor interface implementation
        typedef const Index* ImageIterator;

      private:
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
          uint64_t elements64;
          file.read((char *)&rows, sizeof(uint64_t));
          file.read((char *)&columns, sizeof(uint64_t));
          file.read((char *)&elements64, sizeof(uint64_t));
          Index elements = (Index)elements64;

          this->_scalar_index.at(0) = Index(rows * columns);
          this->_scalar_index.push_back(Index(rows));
          this->_scalar_index.push_back(Index(columns));
          this->_scalar_index.push_back(elements);
          this->_scalar_dt.push_back(DT_(0));

          uint64_t * ccol_ind = new uint64_t[elements];
          file.read((char *)ccol_ind, (long)(elements * sizeof(uint64_t)));
          Index * tcol_ind = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>((elements) * sizeof(Index));
          for (Index i(0) ; i < elements ; ++i)
            tcol_ind[i] = Index(ccol_ind[i]);
          delete[] ccol_ind;

          uint64_t * crow_ptr = new uint64_t[std::size_t(rows + 1)];
          file.read((char *)crow_ptr, (long)((rows + 1) * sizeof(uint64_t)));
          Index * trow_ptr = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(Index(rows + 1) * sizeof(Index));
          for (Index i(0) ; i < rows + 1 ; ++i)
            trow_ptr[i] = Index(crow_ptr[i]);
          delete[] crow_ptr;

          uint64_t * crow_ptr_end = new uint64_t[std::size_t(rows)];
          file.read((char *)crow_ptr_end, (long)((rows) * sizeof(uint64_t)));
          Index * trow_ptr_end = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(Index(rows) * sizeof(Index));
          for (Index i(0) ; i < rows ; ++i)
            trow_ptr_end[i] = Index(crow_ptr_end[i]);
          delete[] crow_ptr_end;

          double * cval = new double[elements];
          file.read((char *)cval, (long)(elements * sizeof(double)));
          DT_ * tval = (DT_*)MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>((elements) * sizeof(DT_));
          for (Index i(0) ; i < elements ; ++i)
            tval[i] = DT_(cval[i]);
          delete[] cval;

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(elements * sizeof(DT_)));
          this->_elements_size.push_back(elements);
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>(elements * sizeof(Index)));
          this->_indices_size.push_back(elements);
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1) + 1)* sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1) + 1);
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1))* sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval, elements * sizeof(DT_));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(0), tcol_ind, elements * sizeof(Index));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(1), trow_ptr, (this->_scalar_index.at(1)  + 1) * sizeof(Index));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(2), trow_ptr_end, (this->_scalar_index.at(1)) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(tval);
          MemoryPool<Mem::Main>::instance()->release_memory(tcol_ind);
          MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr);
          MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr_end);

          this->_scalar_index.push_back(MemoryPool<Mem_>::instance()->generate_hash(this->_indices.at(1), (this->_scalar_index.at(1) + 1) * sizeof(Index)));
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Mem_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit SparseMatrixCSR() :
          Container<Mem_, DT_> (0)
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
         * \param[in] other The source matrix in ELL format.
         *
         * Creates a CSR matrix based on the ELL source matrix.
         */
        template <typename Arch2_>
        explicit SparseMatrixCSR(const SparseMatrixELL<Arch2_, DT_> & other_orig) :
          Container<Mem_, DT_>(other_orig.size())
        {
          CONTEXT("When creating SparseMatrixCSR");
          this->_scalar_index.push_back(other_orig.rows());
          this->_scalar_index.push_back(other_orig.columns());
          this->_scalar_index.push_back(other_orig.used_elements());
          this->_scalar_dt.push_back(other_orig.zero_element());

          SparseMatrixELL<Mem::Main, DT_> ccother(other_orig);
          SparseMatrixCOO<Mem::Main, DT_> cother(ccother);

          DT_ *tval = (DT_*)MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3) * sizeof(DT_));
          Index * tcol_ind = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(3) * sizeof(Index));
          Index * trow_ptr = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1) + 1) * sizeof(Index));
          Index * trow_ptr_end = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1)) * sizeof(Index));

          Index ait(0);
          Index current_row(0);
          trow_ptr[current_row] = 0;
          for (Index it(0) ; it < cother.used_elements() ; ++it)
          {
            Index row(cother.row()[it]);
            Index column(cother.column()[it]);

            if (current_row < row)
            {
              trow_ptr_end[current_row] = ait;
              for (unsigned long i(current_row + 1) ; i < row ; ++i)
              {
                trow_ptr[i] = ait;
                trow_ptr_end[i] = ait;
              }
              current_row = row;
              trow_ptr[current_row] = ait;
            }
            tval[ait] = cother.val()[it];
            tcol_ind[ait] = column;
            ++ait;
          }
          trow_ptr_end[current_row] = ait;
          for (unsigned long i(current_row + 1) ; i < this->_scalar_index.at(1) ; ++i)
          {
            trow_ptr[i] = ait;
            trow_ptr_end[i] = ait;
          }
          trow_ptr[this->_scalar_index.at(1)] = ait;

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1) + 1)* sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1) + 1);
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1))* sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval, this->_scalar_index.at(3) * sizeof(DT_));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(0), tcol_ind, this->_scalar_index.at(3) * sizeof(Index));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(1), trow_ptr, (this->_scalar_index.at(1)  + 1) * sizeof(Index));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(2), trow_ptr_end, (this->_scalar_index.at(1)) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(tval);
          MemoryPool<Mem::Main>::instance()->release_memory(tcol_ind);
          MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr);
          MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr_end);

          this->_scalar_index.push_back(MemoryPool<Mem_>::instance()->generate_hash(this->_indices.at(1), (this->_scalar_index.at(1) + 1) * sizeof(Index)));
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix in COO format.
         *
         * Creates a CSR matrix based on the COO source matrix.
         */
        template <typename Arch2_>
        explicit SparseMatrixCSR(const SparseMatrixCOO<Arch2_, DT_> & other) :
          Container<Mem_, DT_>(other.size())
        {
          CONTEXT("When creating SparseMatrixCSR");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_dt.push_back(other.zero_element());

          SparseMatrixCOO<Mem::Main, DT_> cother(other);

          DT_ * tval = (DT_*)MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3) * sizeof(DT_));
          Index * tcol_ind = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(3) * sizeof(Index));
          Index * trow_ptr = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1) + 1) * sizeof(Index));
          Index * trow_ptr_end = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1)) * sizeof(Index));

          Index ait(0);
          Index current_row(0);
          trow_ptr[current_row] = 0;
          for (Index it(0) ; it < cother.used_elements() ; ++it)
          {
            Index row(cother.row()[it]);
            Index column(cother.column()[it]);

            if (current_row < row)
            {
              trow_ptr_end[current_row] = ait;
              for (unsigned long i(current_row + 1) ; i < row ; ++i)
              {
                trow_ptr[i] = ait;
                trow_ptr_end[i] = ait;
              }
              current_row = row;
              trow_ptr[current_row] = ait;
            }
            tval[ait] = cother.val()[it];
            tcol_ind[ait] = column;
            ++ait;
          }
          trow_ptr_end[current_row] = ait;
          for (unsigned long i(current_row + 1) ; i < this->_scalar_index.at(1) ; ++i)
          {
            trow_ptr[i] = ait;
            trow_ptr_end[i] = ait;
          }
          trow_ptr[this->_scalar_index.at(1)] = ait;

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1) + 1)* sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1) + 1);
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->template allocate_memory<Index>((this->_scalar_index.at(1))* sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval, this->_scalar_index.at(3) * sizeof(DT_));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(0), tcol_ind, this->_scalar_index.at(3) * sizeof(Index));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(1), trow_ptr, (this->_scalar_index.at(1)  + 1) * sizeof(Index));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(2), trow_ptr_end, (this->_scalar_index.at(1)) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(tval);
          MemoryPool<Mem::Main>::instance()->release_memory(tcol_ind);
          MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr);
          MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr_end);

          this->_scalar_index.push_back(MemoryPool<Mem_>::instance()->generate_hash(this->_indices.at(1), (this->_scalar_index.at(1) + 1) * sizeof(Index)));
        }

        /**
         * \brief Constructor
         *
         * \param[in] filename The source file in binary CSR format.
         *
         * Creates a CSR matrix based on the source file.
         */
        explicit SparseMatrixCSR(String filename) :
          Container<Mem_, DT_>(0)
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
          Container<Mem_, DT_>(0)
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
        explicit SparseMatrixCSR(Index rows, Index columns, const DenseVector<Mem_, Index> & col_ind, const DenseVector<Mem_, DT_> & val, const DenseVector<Mem_, Index> & row_ptr, const DenseVector<Mem_, Index> & row_ptr_end) :
          Container<Mem_, DT_>(rows * columns)
        {
          CONTEXT("When creating SparseMatrixCSR");
          this->_scalar_index.push_back(rows);
          this->_scalar_index.push_back(columns);
          // \TODO use real element count - be aware of non used elements
          this->_scalar_index.push_back(val.size());
          this->_scalar_dt.push_back(DT_(0));

          this->_elements.push_back(val.get_elements().at(0));
          this->_elements_size.push_back(val.size());
          this->_indices.push_back(col_ind.get_elements().at(0));
          this->_indices_size.push_back(col_ind.size());
          this->_indices.push_back(row_ptr.get_elements().at(0));
          this->_indices_size.push_back(row_ptr.size());
          this->_indices.push_back(row_ptr_end.get_elements().at(0));
          this->_indices_size.push_back(row_ptr_end.size());

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));

          this->_scalar_index.push_back(MemoryPool<Mem_>::instance()->generate_hash(this->_indices.at(1), (this->_scalar_index.at(1) + 1) * sizeof(Index)));
        }

        /**
         * \brief Constructor
         *
         * \param[in] graph The Graph, the matrix will be created from.
         *
         * Creates a matrix from a given graph.
         */
        explicit SparseMatrixCSR(const Adjacency::Graph & graph) :
          Container<Mem_, DT_>(graph.get_num_nodes_domain() * graph.get_num_nodes_image())
        {
          CONTEXT("When creating SparseMatrixCSR");
          this->_scalar_index.push_back(graph.get_num_nodes_domain());
          this->_scalar_index.push_back(graph.get_num_nodes_image());
          this->_scalar_index.push_back(graph.get_num_indices());
          this->_scalar_dt.push_back(DT_(0));

          const Index* dom_ptr = graph.get_domain_ptr();
          const Index* dom_end = graph.get_domain_end();
          const Index* img_idx = graph.get_image_idx();
          if(dom_end == nullptr)
          {
            dom_end = &dom_ptr[1];
          }

          DenseVector<Mem_, Index> col_ind(this->_scalar_index.at(3));
          DenseVector<Mem_, DT_> val(this->_scalar_index.at(3));
          DenseVector<Mem_, Index> row_ptr(this->_scalar_index.at(1) + 1);
          DenseVector<Mem_, Index> row_ptr_end(this->_scalar_index.at(1));

          Index* prow_ptr = row_ptr.elements();
          for(Index i(0); i <= this->_scalar_index.at(1); ++i)
          {
            prow_ptr[i] = dom_ptr[i];
          }

          Index* prow_end = row_ptr_end.elements();
          for(Index i(0); i < this->_scalar_index.at(1); ++i)
          {
            prow_end[i] = dom_end[i];
          }

          Index* pcol_ind = col_ind.elements();
          for(Index i(0); i < this->_scalar_index.at(3); ++i)
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

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));

          this->_scalar_index.push_back(MemoryPool<Mem_>::instance()->generate_hash(this->_indices.at(1), (this->_scalar_index.at(1) + 1) * sizeof(Index)));
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a shallow copy of a given matrix.
         */
        SparseMatrixCSR(const SparseMatrixCSR<Mem_, DT_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseMatrixCSR");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(other.hash());
          this->_scalar_dt.push_back(other.zero_element());
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
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseMatrixCSR");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(other.hash());
          this->_scalar_dt.push_back(other.zero_element());
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        SparseMatrixCSR<Mem_, DT_> clone()
        {
          CONTEXT("When cloning SparseMatrixCSR");

          SparseMatrixCSR<Mem_, DT_> t;
          ((Container<Mem_, DT_>&)t).clone(*this);
          return t;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns another matrix to the target matrix.
         */
        SparseMatrixCSR<Mem_, DT_> & operator= (const SparseMatrixCSR<Mem_, DT_> & other)
        {
          CONTEXT("When assigning SparseMatrixCSR");

          this->assign(other);

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
        SparseMatrixCSR<Mem_, DT_> & operator= (const SparseMatrixCSR<Arch2_, DT2_> & other)
        {
          CONTEXT("When assigning SparseMatrixCSR");

          this->assign(other);

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

          Index * col_ind = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_indices_size.at(0) * sizeof(Index));
          MemoryPool<Mem_>::template download<Index>(col_ind, this->_indices.at(0), this->_indices_size.at(0) * sizeof(Index));
          uint64_t * ccol_ind = new uint64_t[this->_indices_size.at(0)];
          for (Index i(0) ; i < this->_indices_size.at(0) ; ++i)
            ccol_ind[i] = col_ind[i];
          MemoryPool<Mem::Main>::instance()->release_memory(col_ind);

          Index * row_ptr = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_indices_size.at(1) * sizeof(Index));
          MemoryPool<Mem_>::template download<Index>(row_ptr, this->_indices.at(1), this->_indices_size.at(1) * sizeof(Index));
          uint64_t * crow_ptr = new uint64_t[this->_indices_size.at(1)];
          for (Index i(0) ; i < this->_indices_size.at(1) ; ++i)
            crow_ptr[i] = row_ptr[i];
          MemoryPool<Mem::Main>::instance()->release_memory(row_ptr);

          Index * row_ptr_end = (Index*)MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_indices_size.at(2) * sizeof(Index));
          MemoryPool<Mem_>::template download<Index>(row_ptr_end, this->_indices.at(2), this->_indices_size.at(2) * sizeof(Index));
          uint64_t * crow_ptr_end = new uint64_t[this->_indices_size.at(2)];
          for (Index i(0) ; i < this->_indices_size.at(2) ; ++i)
            crow_ptr_end[i] = row_ptr_end[i];
          MemoryPool<Mem::Main>::instance()->release_memory(row_ptr_end);

          DT_ * val = (DT_*)MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(this->_elements_size.at(0) * sizeof(DT_));
          MemoryPool<Mem_>::template download<DT_>(val, this->_elements.at(0), this->_elements_size.at(0) * sizeof(DT_));
          double * cval = new double[this->_elements_size.at(0)];
          for (Index i(0) ; i < this->_elements_size.at(0) ; ++i)
            cval[i] = Type::Traits<DT_>::to_double(val[i]);
          MemoryPool<Mem::Main>::instance()->release_memory(val);

          uint64_t rows(this->_scalar_index.at(1));
          uint64_t columns(this->_scalar_index.at(2));
          uint64_t elements(this->_indices_size.at(0));
          file.write((const char *)&rows, sizeof(uint64_t));
          file.write((const char *)&columns, sizeof(uint64_t));
          file.write((const char *)&elements, sizeof(uint64_t));
          file.write((const char *)ccol_ind, (long)(elements * sizeof(uint64_t)));
          file.write((const char *)crow_ptr, (long)((rows + 1) * sizeof(uint64_t)));
          file.write((const char *)crow_ptr_end, (long)(rows * sizeof(uint64_t)));
          file.write((const char *)cval, (long)(elements * sizeof(double)));

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
          for (Index row(0) ; row < this->_scalar_index.at(1) ; ++row)
          {
            const Index end(temp.row_ptr_end()[row]);
            for (Index i(temp.row_ptr()[row]) ; i < end ; ++i)
            {
              if (temp.val()[i] != DT_(0))
              {
                file << row + 1 << " " << temp.col_ind()[i] + 1 << " " << std::scientific << Type::Traits<DT_>::to_double(temp.val()[i]) << ";" << std::endl;
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

          ASSERT(row < this->_scalar_index.at(1), "Error: " + stringify(row) + " exceeds sparse matrix csr row size " + stringify(this->_scalar_index.at(1)) + " !");
          ASSERT(col < this->_scalar_index.at(2), "Error: " + stringify(col) + " exceeds sparse matrix csr column size " + stringify(this->_scalar_index.at(2)) + " !");

          if (typeid(Mem_) == typeid(Mem::Main))
          {
            for (unsigned long i(this->_indices.at(1)[row]) ; i < this->_indices.at(2)[row] ; ++i)
            {
              if (this->_indices.at(0)[i] == col)
                return this->_elements.at(0)[i];
              if (this->_indices.at(0)[i] > col)
                return this->_scalar_dt.at(0);
            }
            return this->_scalar_dt.at(0);
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
        Index used_elements() const
        {
          return this->_scalar_index.at(3);
        }

        /**
         * \brief Retrieve column indices array.
         *
         * \returns Column indices array.
         */
        Index * col_ind() const
        {
          return this->_indices.at(0);
        }

        /**
         * \brief Retrieve non zero element array.
         *
         * \returns Non zero element array.
         */
        DT_ * val() const
        {
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve row start index array.
         *
         * \returns Row start index array.
         */
        Index * row_ptr() const
        {
          return this->_indices.at(1);
        }

        /**
         * \brief Retrieve row end index array.
         *
         * \returns Row end index array.
         */
        Index * row_ptr_end() const
        {
          return this->_indices.at(2);
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
         * \brief Retrieve layout hash.
         *
         * \returns Hash value of matrix layout.
         */
        unsigned long hash() const
        {
          return this->_scalar_index.at(4);
        }

        /* ******************************************************************* */
        /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
        /* ******************************************************************* */
      public:
        /** \copydoc Adjactor::get_num_nodes_domain() */
        inline Index get_num_nodes_domain() const
        {
          return this->_scalar_index.at(1);
        }

        /** \copydoc Adjactor::get_num_nodes_image() */
        inline Index get_num_nodes_image() const
        {
          return this->_scalar_index.at(2);
        }

        /** \copydoc Adjactor::image_begin() */
        inline ImageIterator image_begin(Index domain_node) const
        {
          ASSERT(domain_node < this->_scalar_index.at(1), "Domain node index out of range");
          return &this->_indices.at(0)[this->_indices.at(1)[domain_node]];
        }

        /** \copydoc Adjactor::image_end() */
        inline ImageIterator image_end(Index domain_node) const
        {
          CONTEXT("Graph::image_end()");
          ASSERT(domain_node < this->_scalar_index.at(1), "Domain node index out of range");
          return &this->_indices.at(0)[this->_indices.at(2)[domain_node]];
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
    template <typename Mem_, typename Arch2_, typename DT_> bool operator== (const SparseMatrixCSR<Mem_, DT_> & a, const SparseMatrixCSR<Arch2_, DT_> & b)
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

      for (Index i(0) ; i < a.used_elements() ; ++i)
      {
        if (MemoryPool<Mem_>::get_element(a.col_ind(), i) != MemoryPool<Arch2_>::get_element(b.col_ind(), i))
          return false;
        if (MemoryPool<Mem_>::get_element(a.val(), i) != MemoryPool<Arch2_>::get_element(b.val(), i))
          return false;
      }
      for (Index i(0) ; i < a.rows() + 1; ++i)
      {
        if (MemoryPool<Mem_>::get_element(a.row_ptr(), i) != MemoryPool<Arch2_>::get_element(b.row_ptr(), i))
          return false;
      }
      for (Index i(0) ; i < a.rows(); ++i)
      {
        if (MemoryPool<Mem_>::get_element(a.row_ptr_end(), i) != MemoryPool<Arch2_>::get_element(b.row_ptr_end(), i))
          return false;
      }

      return true;
    }

    /**
     * \brief SparseMatrixCSR streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Mem_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixCSR<Mem_, DT_> & b)
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
