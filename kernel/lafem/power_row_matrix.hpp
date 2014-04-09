#pragma once
#ifndef KERNEL_LAFEM_POWER_ROW_MATRIX_HPP
#define KERNEL_LAFEM_POWER_ROW_MATRIX_HPP 1

// includes, FEAST
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Power-Row-Matrix meta class template
     *
     * This class template implements a horizontal composition of \e n sub-matrices of the same class.
     * This can be interpreted as a dense 1-by-m matrix of other matrices.
     *
     * \tparam SubType_
     * The type of the sub-matrix.
     *
     * \tparam blocks_
     * The bumber of sub-matrix blocks.
     *
     * \author Peter Zajac
     */
    template<
      typename SubType_,
      Index blocks_>
    class PowerRowMatrix :
      protected PowerRowMatrix<SubType_, blocks_-1>
    {
      // declare this class template as a friend for recursive inheritance
      template<typename, Index>
      friend class PowerRowMatrix;

      /// base-class typedef
      typedef PowerRowMatrix<SubType_, blocks_-1> BaseClass;

    public:
      /// sub-matrix type
      typedef SubType_ SubMatrixType;
      /// sub-matrix memory type
      typedef typename SubMatrixType::MemType MemType;
      /// sub-matrix data type
      typedef typename SubMatrixType::DataType DataType;
      /// sub-matrix index type
      typedef typename SubMatrixType::IndexType IndexType;
      /// sub-matrix layout type
      static constexpr SparseLayoutType LayoutType = SubMatrixType::LayoutType;
      /// Compatible L-vector type
      typedef typename SubMatrixType::VectorTypeL VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, blocks_> VectorTypeR;

      /// dummy enum
      enum
      {
        /// number of row blocks (vertical size)
        num_row_blocks = 1,
        /// number of column blocks (horizontal size)
        num_col_blocks = blocks_
      };

    protected:
      /// the last sub-matrix
      SubMatrixType _sub_matrix;

      /// base-class constructor; this one is protected for a reason
      explicit PowerRowMatrix(BaseClass&& other_base, SubMatrixType&& last_sub) :
        BaseClass(std::move(other_base)),
        _sub_matrix(std::move(last_sub))
      {
      }

    public:
      /// default ctor
      PowerRowMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerRowMatrix(const SparseLayout<MemType, IndexType, LayoutType>& layout) :
        BaseClass(layout),
        _sub_matrix(layout)
      {
      }

      /// move ctor
      PowerRowMatrix(PowerRowMatrix&& other) :
        BaseClass(static_cast<BaseClass&&>(other)),
        _sub_matrix(std::move(other._sub_matrix))
      {
      }

      /// move-assign operator
      PowerRowMatrix& operator=(PowerRowMatrix&& other)
      {
        base().operator=(static_cast<BaseClass&&>(other));
        _sub_matrix = std::move(other._sub_matrix);
        return *this;
      }

      /// deleted copy-ctor
      PowerRowMatrix(const PowerRowMatrix&) = delete;
      /// deleted copy-assign operator
      PowerRowMatrix& operator=(const PowerRowMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerRowMatrix()
      {
      }

      /**
       * \brief Creates and returns a deep copy of this matrix.
       */
      PowerRowMatrix clone() const
      {
        return PowerRowMatrix(base().clone(), _sub_matrix.clone());
      }

      /**
       * \brief Returns a sub-matrix block.
       *
       * \tparam i_
       * The row index of the sub-matrix block that is to be returned.
       *
       * \tparam j_
       * The column index of the sub-matrix block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-matrix at position <em>(i_,j_)</em>.
       */
      template<Index i_, Index j_>
      SubMatrixType& at()
      {
        static_assert(i_ == 0, "invalid sub-matrix index");
        static_assert(j_ < blocks_, "invalid sub-matrix index");
        return static_cast<PowerRowMatrix<SubType_, j_+1>&>(*this)._sub_matrix;
      }

      /** \copydoc at() */
      template<Index i_, Index j_>
      const SubMatrixType& at() const
      {
        static_assert(i_ == 0, "invalid sub-matrix index");
        static_assert(j_ < blocks_, "invalid sub-matrix index");
        return static_cast<const PowerRowMatrix<SubType_, j_+1>&>(*this)._sub_matrix;
      }

      /// \cond internal
      SubMatrixType& last()
      {
        return _sub_matrix;
      }

      const SubMatrixType& last() const
      {
        return _sub_matrix;
      }

      PowerRowMatrix<SubType_, blocks_-1>& base()
      {
        return static_cast<BaseClass&>(*this);
      }

      const PowerRowMatrix<SubType_, blocks_-1>& base() const
      {
        return static_cast<const BaseClass&>(*this);
      }

      Index row_blocks() const
      {
        return Index(num_row_blocks);
      }

      Index col_blocks() const
      {
        return Index(num_col_blocks);
      }
      /// \endcond

      /// Returns the total number of rows in this matrix.
      Index rows() const
      {
        return last().rows();
      }

      /// Returns the total number of columns in this matrix.
      Index columns() const
      {
        return base().columns() + last().columns();
      }

      /// Returns the total number of non-zeros in this matrix.
      Index used_elements() const
      {
        return base().used_elements() + last().used_elements();
      }

      /**
       * \brief Clears this matrix.
       *
       * \param[in] value
       * The value to which the matrix' entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        base().format(value);
        last().format(value);
      }

      /**
       * \brief Applies this matrix onto a vector.
       *
       * This function performs
       *  \f[r \leftarrow this\cdot x \f]
       *
       * \param[out] r
       * The vector the receives the result.
       *
       * \param[in] x
       * The multiplicant vector.
       */
      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x)
      {
        base().template apply<Algo_>(r, x.base());
        last().template apply<Algo_>(r, x.last(), r, DataType(1));
      }

      /**
       * \brief Applies this matrix onto a vector.
       *
       * This function performs
       *  \f[r \leftarrow y + \alpha\cdot this\cdot x \f]
       *
       * \param[out] r
       * The vector the receives the result.
       *
       * \param[in] x
       * The multiplicant vector.
       *
       * \param[in] y
       * The summand vector
       * \param[in] alpha A scalar to scale the product with.
       */
      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1))
      {
        base().template apply<Algo_>(r, x.base(), y, alpha);
        last().template apply<Algo_>(r, x.last(), r, alpha);
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return base().create_vector_l();
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(base().create_vector_r(), last().create_vector_r());
      }
    };

    /// \cond internal
    template<typename SubType_>
    class PowerRowMatrix<SubType_, 1>
    {
      template<typename, Index>
      friend class PowerRowMatrix;

    public:
      typedef SubType_ SubMatrixType;
      typedef typename SubMatrixType::MemType MemType;
      typedef typename SubMatrixType::DataType DataType;
      typedef typename SubMatrixType::IndexType IndexType;
      /// sub-matrix layout type
      static constexpr SparseLayoutType LayoutType = SubMatrixType::LayoutType;
      /// Compatible L-vector type
      typedef typename SubMatrixType::VectorTypeL VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, 1> VectorTypeR;

      enum
      {
        num_row_blocks = 1,
        num_col_blocks = 1
      };

    protected:
      SubMatrixType _sub_matrix;

      /// base-class constructor; this one is protected for a reason
      explicit PowerRowMatrix(SubMatrixType&& last_sub) :
        _sub_matrix(std::move(last_sub))
      {
      }

    public:
      /// default ctor
      PowerRowMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerRowMatrix(const SparseLayout<MemType, IndexType, LayoutType>& layout) :
        _sub_matrix(layout)
      {
      }

      /// move ctor
      PowerRowMatrix(PowerRowMatrix&& other) :
        _sub_matrix(std::move(other._sub_matrix))
      {
      }

      /// move-assign operator
      PowerRowMatrix& operator=(PowerRowMatrix&& other)
      {
        _sub_matrix = std::move(other._sub_matrix);
        return *this;
      }

      /// deleted copy-ctor
      PowerRowMatrix(const PowerRowMatrix&) = delete;
      /// deleted copy-assign operator
      PowerRowMatrix& operator=(const PowerRowMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerRowMatrix()
      {
      }

      PowerRowMatrix clone() const
      {
        return PowerRowMatrix(_sub_matrix.clone());
      }

      template<Index i, Index j>
      SubMatrixType& at()
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _sub_matrix;
      }

      template<Index i, Index j>
      const SubMatrixType& at() const
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _sub_matrix;
      }

      SubMatrixType& last()
      {
        return _sub_matrix;
      }

      const SubMatrixType& last() const
      {
        return _sub_matrix;
      }

      Index row_blocks() const
      {
        return Index(1);
      }

      Index col_blocks() const
      {
        return Index(1);
      }

      Index rows() const
      {
        return last().rows();
      }

      Index columns() const
      {
        return last().columns();
      }

      Index used_elements() const
      {
        return last().used_elements();
      }

      void format(DataType value = DataType(0))
      {
        last().format(value);
      }

      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x)
      {
        last().template apply<Algo_>(r, x.last());
      }

      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1))
      {
        last().template apply<Algo_>(r, x.last(), y, alpha);
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return last().create_vector_l();
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(last().create_vector_r());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_ROW_MATRIX_HPP
