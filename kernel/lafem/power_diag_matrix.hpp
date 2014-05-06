#pragma once
#ifndef KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP
#define KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP 1

// includes, FEAST
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Power-Diag-Matrix meta class template
     *
     * This class template implements a diagonal composition of \e n sub-matrices of the same class.
     * This can be interpreted as a diagonal m-by-m matrix of other matrices.
     *
     * \tparam SubType_
     * The type of the sub-matrix.
     *
     * \tparam blocks_
     * The number of sub-matrix blocks.
     *
     * \author Peter Zajac
     */
    template<
      typename SubType_,
      Index blocks_>
    class PowerDiagMatrix :
      protected PowerDiagMatrix<SubType_, blocks_-1>
    {
      // declare this class template as a friend for recursive inheritance
      template<typename, Index>
      friend class PowerDiagMatrix;

      /// base-class typedef
      typedef PowerDiagMatrix<SubType_, blocks_-1> BaseClass;

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
      static constexpr SparseLayoutId layout_id = SubMatrixType::layout_id;
      /// Compatible L-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeL, blocks_> VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, blocks_> VectorTypeR;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_, typename IT2_ = IndexType>
      using ContainerType = class PowerDiagMatrix<typename SubType_::template ContainerType<Mem2_, DT2_, IT2_>, blocks_>;

      /// dummy enum
      enum
      {
        /// number of row blocks (vertical size)
        num_row_blocks = blocks_,
        /// number of column blocks (horizontal size)
        num_col_blocks = blocks_
      };

    protected:
      /// the last sub-matrix
      SubMatrixType _sub_matrix;

      /// base-class constructor; this one is protected for a reason
      explicit PowerDiagMatrix(BaseClass&& other_base, SubMatrixType&& last_sub) :
        BaseClass(std::move(other_base)),
        _sub_matrix(std::move(last_sub))
      {
      }

    public:
      /// default ctor
      PowerDiagMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerDiagMatrix(const SparseLayout<MemType, IndexType, layout_id>& layout) :
        BaseClass(layout),
        _sub_matrix(layout)
      {
      }

      /// move ctor
      PowerDiagMatrix(PowerDiagMatrix&& other) :
        BaseClass(static_cast<BaseClass&&>(other)),
        _sub_matrix(std::move(other._sub_matrix))
      {
      }

      /// move-assign operator
      PowerDiagMatrix& operator=(PowerDiagMatrix&& other)
      {
        base().operator=(static_cast<BaseClass&&>(other));
        _sub_matrix = std::move(other._sub_matrix);
        return *this;
      }

      /// deleted copy-ctor
      PowerDiagMatrix(const PowerDiagMatrix&) = delete;
      /// deleted copy-assign operator
      PowerDiagMatrix& operator=(const PowerDiagMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerDiagMatrix()
      {
      }

      /**
       * \brief Creates and returns a deep copy of this matrix.
       */
      PowerDiagMatrix clone() const
      {
        return PowerDiagMatrix(base().clone(), _sub_matrix.clone());
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
        static_assert(i_ == j_, "invalid sub-matrix index");
        static_assert(i_ < blocks_, "invalid sub-matrix index");
        return static_cast<PowerDiagMatrix<SubType_, i_+1>&>(*this)._sub_matrix;
      }

      /** \copydoc at() */
      template<Index i_, Index j_>
      const SubMatrixType& at() const
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        static_assert(i_ < blocks_, "invalid sub-matrix index");
        return static_cast<const PowerDiagMatrix<SubType_, i_+1>&>(*this)._sub_matrix;
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

      PowerDiagMatrix<SubType_, blocks_-1>& base()
      {
        return static_cast<BaseClass&>(*this);
      }

      const PowerDiagMatrix<SubType_, blocks_-1>& base() const
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
        return base().rows() + last().rows();
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
      void clear(DataType value = DataType(0))
      {
        base().clear(value);
        last().clear(value);
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
        base().template apply<Algo_>(r.base(), x.base());
        last().template apply<Algo_>(r.last(), x.last());
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
        base().template apply<Algo_>(r.base(), x.base(), y.base(), alpha);
        last().template apply<Algo_>(r.last(), x.last(), y.last(), alpha);
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(base().create_vector_l(), last().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(base().create_vector_r(), last().create_vector_r());
      }

      Index get_length_of_line(const Index row) const
      {
        const Index brows(this->base().rows());

        if (row < brows)
        {
          return this->base().get_length_of_line(row);
        }
        else
        {
          return this->last().get_length_of_line(row - brows);
        }
      }

      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set,
                     const Index col_start, const Index stride = 1) const
      {
        const Index brows(this->base().rows());
        const Index bcolumns(this->base().columns());

        if (row < brows)
        {
          this->base().set_line(row, pval_set, pcol_set, col_start, stride);
        }
        else
        {
          this->last().set_line(row - brows, pval_set, pcol_set, col_start + bcolumns, stride);
        }
      }

    };

    /// \cond internal
    template<typename SubType_>
    class PowerDiagMatrix<SubType_, 1>
    {
      template<typename, Index>
      friend class PowerDiagMatrix;

    public:
      typedef SubType_ SubMatrixType;
      typedef typename SubMatrixType::MemType MemType;
      typedef typename SubMatrixType::DataType DataType;
      typedef typename SubMatrixType::IndexType IndexType;
      /// sub-matrix layout type
      static constexpr SparseLayoutId layout_id = SubMatrixType::layout_id;
      /// Compatible L-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeL, 1> VectorTypeL;
      /// Compatible R-vector type
      typedef PowerVector<typename SubMatrixType::VectorTypeR, 1> VectorTypeR;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_, typename IT2_ = IndexType>
      using ContainerType = class PowerDiagMatrix<typename SubType_::template ContainerType<Mem2_, DT2_, IT2_>, 1>;

      enum
      {
        num_row_blocks = 1,
        num_col_blocks = 1
      };

    protected:
      SubMatrixType _sub_matrix;

      /// base-class constructor; this one is protected for a reason
      explicit PowerDiagMatrix(SubMatrixType&& last_sub) :
        _sub_matrix(std::move(last_sub))
      {
      }

    public:
      /// default ctor
      PowerDiagMatrix()
      {
      }

      /// sub-matrix layout ctor
      explicit PowerDiagMatrix(const SparseLayout<MemType, IndexType, layout_id>& layout) :
        _sub_matrix(layout)
      {
      }

      /// move ctor
      PowerDiagMatrix(PowerDiagMatrix&& other) :
        _sub_matrix(std::move(other._sub_matrix))
      {
      }

      /// move-assign operator
      PowerDiagMatrix& operator=(PowerDiagMatrix&& other)
      {
        _sub_matrix = std::move(other._sub_matrix);
        return *this;
      }

      /// deleted copy-ctor
      PowerDiagMatrix(const PowerDiagMatrix&) = delete;
      /// deleted copy-assign operator
      PowerDiagMatrix& operator=(const PowerDiagMatrix&) = delete;

      /// virtual destructor
      virtual ~PowerDiagMatrix()
      {
      }

      PowerDiagMatrix clone() const
      {
        return PowerDiagMatrix(_sub_matrix.clone());
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

      void clear(DataType value = DataType(0))
      {
        last().clear(value);
      }

      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x)
      {
        last().template apply<Algo_>(r.last(), x.last());
      }

      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1))
      {
        last().template apply<Algo_>(r.last(), x.last(), y.last(), alpha);
      }
      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(last().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(last().create_vector_r());
      }

      Index get_length_of_line(const Index row) const
      {
        return this->last().get_length_of_line(row);
      }

      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        this->last().set_line(row, pval_set, pcol_set, col_start, stride);
      }

    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_DIAG_MATRIX_HPP
