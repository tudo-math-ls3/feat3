#pragma once
#ifndef KERNEL_LAFEM_DENSE_MATRIX_HPP
#define KERNEL_LAFEM_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/matrix_base.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Dense data matrix class template.
     *
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a matrix of continuous data in memory. \n\n
     * Data survey: \n
     * _elements[0]: raw number values \n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class DenseMatrix : public Container<Mem_, DT_, IT_>, public MatrixBase
    {
      private:
        Index & _rows()
        {
          return _rows();
        }

        Index & _columns()
        {
          return _columns();
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our indextype
        typedef IT_ IndexType;
        /// Our memory architecture type
        typedef Mem_ MemType;
        /// Our 'base' class type
        template <typename Mem2_, typename DT2_, typename IT2_ = IT_>
        using ContainerType = class DenseMatrix<Mem2_, DT2_, IT2_>;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit DenseMatrix() :
          Container<Mem_, DT_, IT_> (0)
        {
          CONTEXT("When creating DenseMatrix");

          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
        }

        /**
         * \brief Constructor
         *
         * \param[in] rows The row count of the created matrix.
         * \param[in] columns The column count of the created matrix.
         *
         * Creates a matrix with given dimensions.
         */
        explicit DenseMatrix(Index rows_in, Index columns_in) :
          Container<Mem_, DT_, IT_>(rows_in * columns_in)
        {
          CONTEXT("When creating DenseMatrix");

          this->_scalar_index.at(0) = rows_in * columns_in;
          this->_scalar_index.push_back(rows_in);
          this->_scalar_index.push_back(columns_in);

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(0)));
          this->_elements_size.push_back(this->_scalar_index.at(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] rows The row count of the created matrix.
         * \param[in] columns The column count of the created matrix.
         * \param[in] value The value, each element will be set to.
         *
         * Creates a matrix with given dimensions and value.
         */
        explicit DenseMatrix(Index rows_in, Index columns_in, DT_ value) :
          Container<Mem_, DT_, IT_>(rows_in * columns_in)
        {
          CONTEXT("When creating DenseMatrix");

          this->_scalar_index.at(0) = rows_in * columns_in;
          this->_scalar_index.push_back(rows_in);
          this->_scalar_index.push_back(columns_in);
          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(0)));
          this->_elements_size.push_back(this->_scalar_index.at(0));
          MemoryPool<Mem_>::instance()->set_memory(this->_elements.at(0), value, this->_scalar_index.at(0));
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source matrix.
         *
         * Moves a given matrix to this matrix.
         */
        DenseMatrix(DenseMatrix && other) :
          Container<Mem_, DT_, IT_>(std::forward<DenseMatrix>(other))
        {
          CONTEXT("When moving DenseMatrix");
        }

        /**
         * \brief Move operator
         *
         * \param[in] other The source matrix.
         *
         * Moves another matrix to the target matrix.
         */
        DenseMatrix & operator= (DenseMatrix && other)
        {
          CONTEXT("When moving DenseMatrix");

          this->move(std::forward<DenseMatrix>(other));

          return *this;
        }

        /** \brief Clone operation
         *
         * Create a deep copy of itself.
         *
         * \return A deep copy of itself.
         *
         */
        DenseMatrix clone() const
        {
          DenseMatrix t;
          t.clone(*this);
          return t;
        }

        using Container<Mem_, DT_, IT_>::clone;

        /**
         * \brief Conversion method
         *
         * \param[in] other The source matrix.
         *
         * Use source matrix content as content of current matrix
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const DenseMatrix<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting DenseMatrix");
          this->assign(other);
        }

        /**
         * \brief Get a pointer to the data array.
         *
         * \returns Pointer to the data array.
         */
        DT_ * elements()
        {
          return this->_elements.at(0);
        }

        DT_ const * elements() const
        {
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         *
         * \returns Specific matrix element.
         */
        const DT_ operator()(Index row, Index col) const
        {
          CONTEXT("When retrieving DenseMatrix element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds dense matrix row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds dense matrix column size " + stringify(this->columns()) + " !");
          return MemoryPool<Mem_>::get_element(this->_elements.at(0), row * this->columns() + col);
        }

        /**
         * \brief Set specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         * \param[in] value The value to be set.
         */
        void operator()(Index row, Index col, DT_ value)
        {
          CONTEXT("When setting DenseMatrix element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds dense matrix row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds dense matrix column size " + stringify(this->columns()) + " !");
          MemoryPool<Mem_>::set_memory(this->_elements.at(0) + row * this->columns() + col, value);
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
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String name()
        {
          return "DenseMatrix";
        }

        /**
         * \brief Performs \f$this \leftarrow x\f$.
         *
         * \param[in] x The Matrix to be copied.
         */
        void copy(const DenseMatrix & x)
        {
          this->_copy_content(x);
        }

        /**
         * \brief Performs \f$this \leftarrow x\f$.
         *
         * \param[in] x The Matrix to be copied.
         */
        template <typename Mem2_>
        void copy(const DenseMatrix<Mem2_, DT_, IT_> & x)
        {
          this->_copy_content(x);
        }

        Index get_length_of_line(const Index /*row*/) const
        {
          return this->columns();
        }
      ///@}
    };

    /**
     * \brief DenseMatrix comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Mem_, typename Mem2_, typename DT_, typename IT_> bool operator== (const DenseMatrix<Mem_, DT_, IT_> & a, const DenseMatrix<Mem2_, DT_, IT_> & b)
    {
      CONTEXT("When comparing DenseMatrices");

      if (a.size() != b.size())
        return false;
      if (a.get_elements().size() != b.get_elements().size())
        return false;
      if (a.get_indices().size() != b.get_indices().size())
        return false;
      if (a.rows() != b.rows())
        return false;
      if (a.columns() != b.columns())
        return false;

      if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && b.get_elements().size() == 0)
        return true;

      bool ret(true);

      DT_ * ta;
      DT_ * tb;

      if(std::is_same<Mem::Main, Mem_>::value)
        ta = (DT_*)a.elements();
      else
      {
        ta = new DT_[a.size()];
        MemoryPool<Mem_>::instance()->template download<DT_>(ta, a.elements(), a.size());
      }
      if(std::is_same<Mem::Main, Mem2_>::value)
        tb = (DT_*)b.elements();
      else
      {
        tb = new DT_[b.size()];
        MemoryPool<Mem2_>::instance()->template download<DT_>(tb, b.elements(), b.size());
      }

      for (Index i(0) ; i < a.size() ; ++i)
        if (ta[i] != tb[i])
        {
          ret = false;
          break;
        }

      if(! std::is_same<Mem::Main, Mem_>::value)
        delete[] ta;
      if(! std::is_same<Mem::Main, Mem2_>::value)
        delete[] tb;

      return ret;
    }

    /**
     * \brief DenseMatrix streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Mem_, typename DT_, typename IT_>
      std::ostream &
      operator<< (std::ostream & lhs, const DenseMatrix<Mem_, DT_, IT_> & b)
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

#endif // KERNEL_LAFEM_DENSE_MATRIX_HPP
