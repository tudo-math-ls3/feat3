#pragma once
#ifndef KERNEL_LAFEM_DENSE_MATRIX_HPP
#define KERNEL_LAFEM_DENSE_MATRIX_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/product_matmat.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/dense_vector.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Dense data matrix class template.
     *
     * \tparam Mem_ The \ref FEAT::Mem "memory architecture" to be used.
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
    class DenseMatrix : public Container<Mem_, DT_, IT_>
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
      /// Compatible L-vector type
      typedef DenseVector<Mem_, DT_, IT_> VectorTypeL;
      /// Compatible R-vector type
      typedef DenseVector<Mem_, DT_, IT_> VectorTypeR;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = class DenseMatrix<Mem2_, DT2_, IT2_>;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit DenseMatrix() :
        Container<Mem_, DT_, IT_> (0)
      {
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
        this->_scalar_index.at(0) = rows_in * columns_in;
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(this->_scalar_index.at(0)));
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
        this->_scalar_index.at(0) = rows_in * columns_in;
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(this->_scalar_index.at(0)));
        this->_elements_size.push_back(this->_scalar_index.at(0));
        MemoryPool<Mem_>::set_memory(this->_elements.at(0), value, this->_scalar_index.at(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] std::vector<char> A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit DenseMatrix(std::vector<char> input) :
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
      DenseMatrix(DenseMatrix && other) :
        Container<Mem_, DT_, IT_>(std::forward<DenseMatrix>(other))
      {
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
        this->move(std::forward<DenseMatrix>(other));

        return *this;
      }

/// \cond nodoxy
      InsertDeepClone ( DenseMatrix )
/// \endcond

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
        ASSERT(row < this->rows());
        ASSERT(col < this->columns());
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
        ASSERT(row < this->rows());
        ASSERT(col < this->columns());
        MemoryPool<Mem_>::set_memory(this->_elements.at(0) + row * this->columns() + col, value);
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
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_dm, input);
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * \param[in] mode FileMode enum, describing the actual container specialisation.
       * \param[out] std::vector<char> A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialise for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialise()
      {
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_dm);
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
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const DenseMatrix & x, bool full = false)
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
      void copy(const DenseMatrix<Mem2_, DT_, IT_> & x, bool full = false)
      {
        this->_copy_content(x, full);
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

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index /*row*/) const
      {
        return this->columns();
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow \alpha x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const DenseMatrix & x, const DT_ alpha)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");

        Arch::Scale<Mem_>::value(this->elements(), x.elements(), alpha, this->used_elements());
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        return Arch::Norm2<Mem_>::value(this->elements(), this->used_elements());
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        if (r.template elements<Perspective::pod>() == x.template elements<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector x and r must not share the same memory!");

        Arch::ProductMatVec<Mem_>::dense(r.elements(), this->elements(),
                                         x.elements(), this->rows(), this->columns());
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
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        if (r.template elements<Perspective::pod>() == x.template elements<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector x and r must not share the same memory!");

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Arch::Defect<Mem_>::dense(r.elements(), y.elements(), this->elements(),
                                    x.elements(), this->rows(), this->columns());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*x
        else
        {
          Arch::Axpy<Mem_>::dense(r.elements(), alpha, y.elements(), this->elements(),
                                  x.elements(), this->rows(), this->columns());
        }
      }

      void multiply(DenseMatrix & x, DenseMatrix & y)
      {
        XASSERTM(x.columns() == y.rows(), "dimension mismatch!");
        XASSERTM(this->rows() == x.rows(), "dimension mismatch!");
        XASSERTM(this->columns() == y.columns(), "dimension mismatch!");

        Arch::ProductMatMat<Mem_>::dense(this->elements(), x.elements(),
                                         y.elements(), this->rows(), this->columns(), x.columns());

      }
      ///@}

      /// \cond internal

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DT_ * const pval_set, IT_ * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        const DT_ * pval(this->elements());

        for (Index i(0); i < columns(); ++i)
        {
          pval_set[i * stride] = pval[columns() * row + i];
          pcol_set[i * stride] = IT_(i) + IT_(col_start);
        }
      }
      /// \endcond

      /**
       * \brief DenseMatrix comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_> friend bool operator== (const DenseMatrix & a, const DenseMatrix<Mem2_, DT_, IT_> & b)
      {
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
          MemoryPool<Mem_>::template download<DT_>(ta, a.elements(), a.size());
        }
        if(std::is_same<Mem::Main, Mem2_>::value)
          tb = (DT_*)b.elements();
        else
        {
          tb = new DT_[b.size()];
          MemoryPool<Mem2_>::template download<DT_>(tb, b.elements(), b.size());
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
      friend std::ostream & operator<< (std::ostream & lhs, const DenseMatrix & b)
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

#endif // KERNEL_LAFEM_DENSE_MATRIX_HPP
