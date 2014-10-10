#pragma once
#ifndef KERNEL_LAFEM_TUPLE_DIAG_MATRIX_HPP
#define KERNEL_LAFEM_TUPLE_DIAG_MATRIX_HPP 1

// includes, FEAST
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/meta_element.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Tuple-Diag-Matrix meta class template
     *
     * This class template implements a diagonal composition of \e n sub-matrices of arbitrary classes.
     * This can be interpreted as a diagonal m-by-m matrix of other matrices.
     *
     * \tparam First_, ...Rest_
     * A sequence of (meta) matrix classes which are to be composed.
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class TupleDiagMatrix
    {
      // declare this class template as a friend for recursive inheritance
      template<typename, typename...>
      friend class TupleDiagMatrix;

      /// rest-class typedef
      typedef TupleDiagMatrix<Rest_...> RestClass;

    public:
      /// sub-matrix memory type
      typedef typename First_::MemType MemType;
      /// sub-matrix data type
      typedef typename First_::DataType DataType;
      /// sub-matrix index type
      typedef typename First_::IndexType IndexType;

#ifdef FEAST_COMPILER_MICROSOFT
      template<typename... RestVecL_>
      using VecLHelper = typename RestClass::template VecLHelper<RestVecL_..., typename First_::VectorTypeL>;
      typedef VecLHelper<> VectorTypeL;
      template<typename... RestVecR_>
      using VecRHelper = typename RestClass::template VecRHelper<RestVecR_..., typename First_::VectorTypeR>;
      typedef VecRHelper<> VectorTypeR;

      /// helper class template for ContainerType template
      template <typename Mem2_, typename DT2_, typename IT2_, typename... RestCont_>
      using ContHelper = typename RestClass::template
        ContHelper<Mem2_, DT2_, IT2_, RestCont_..., typename First_::template ContainerType<Mem2_, DT2_, IT2_>>;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = ContHelper<Mem2_, DT2_, IT2_>;
#else
      /// Compatible L-vector type
      typedef TupleVector<typename First_::VectorTypeL, typename Rest_::VectorTypeL...> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename First_::VectorTypeR, typename Rest_::VectorTypeR...> VectorTypeR;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleDiagMatrix<
        typename First_::template ContainerType<Mem2_, DT2_, IT2_>,
        typename Rest_::template ContainerType<Mem2_, DT2_, IT2_>...>;
#endif

      /// number of row blocks (vertical size)
      static constexpr Index num_row_blocks = RestClass::num_row_blocks + Index(1);
      /// number of column blocks (horizontal size)
      static constexpr Index num_col_blocks = RestClass::num_col_blocks + Index(1);

    protected:
      /// the first sub-matrix
      First_ _first;
      /// the remaining part
      RestClass _rest;

      /// base-class constructor; this one is protected for a reason
      explicit TupleDiagMatrix(First_&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

      /// Returns a list of all sub-matrix type names
      static String sub_name_list()
      {
        return First_::name() + "," + RestClass::sub_name_list();
      }

    public:
      /// default ctor
      TupleDiagMatrix()
      {
      }

      /// Sub-Matrix emplacement constructor
      explicit TupleDiagMatrix(First_&& the_first, Rest_&&... the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest...))
      {
      }

      /// move ctor
      TupleDiagMatrix(TupleDiagMatrix&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      TupleDiagMatrix& operator=(TupleDiagMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleDiagMatrix(const TupleDiagMatrix&) = delete;
      /// deleted copy-assign operator
      TupleDiagMatrix& operator=(const TupleDiagMatrix&) = delete;

      /// virtual destructor
      virtual ~TupleDiagMatrix()
      {
      }

      /**
       * \brief Creates and returns a deep copy of this matrix.
       */
      TupleDiagMatrix clone() const
      {
        return TupleDiagMatrix(_first.clone(), _rest.clone());
      }

      /**
       * \brief Creates and returns the transpose of this matrix.
       */
      TupleDiagMatrix transpose() const
      {
        return TupleDiagMatrix(_first.transpose(), _rest.transpose());
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
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<Index i_, Index j_>
      const typename TupleElement<i_, First_, Rest_...>::Type& at() const
      {
        static_assert(i_ == j_, "invalid sub-matrix index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /// \cond internal
      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      RestClass& rest()
      {
        return _rest;
      }

      const RestClass& rest() const
      {
        return _rest;
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
        return first().rows() + rest().rows();
      }

      /// Returns the total number of columns in this matrix.
      Index columns() const
      {
        return first().columns() + rest().columns();
      }

      /// Returns the total number of non-zeros in this matrix.
      Index used_elements() const
      {
        return first().used_elements() + rest().used_elements();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("TupleDiagMatrix<") + sub_name_list() + ">";
      }

      /**
       * \brief Clears this matrix.
       *
       * \param[in] value
       * The value to which the matrix' entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
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
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().template apply<Algo_>(r.first(), x.first());
        rest().template apply<Algo_>(r.rest(), x.rest());
      }

      template<typename Algo_>
      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, rest().rows(), first().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        first().template apply<Algo_>(r_first, x_first);
        rest().template apply<Algo_>(r_rest, x_rest);
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
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().template apply<Algo_>(r.first(), x.first(), y.first(), alpha);
        rest().template apply<Algo_>(r.rest(), x.rest(), y.rest(), alpha);
      }

      template<typename Algo_>
      void apply(DenseVector<MemType, DataType, IndexType>& r, const DenseVector<MemType, DataType, IndexType>& x,
                 const DenseVector<MemType, DataType, IndexType>& y, DataType alpha = DataType(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        DenseVector<MemType, DataType, IndexType> r_first(r, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> r_rest(r, rest().rows(), first().rows());

        DenseVector<MemType, DataType, IndexType> x_first(x, first().columns(), 0);
        DenseVector<MemType, DataType, IndexType> x_rest(x, rest().columns(), first().columns());

        DenseVector<MemType, DataType, IndexType> y_first(y, first().rows(), 0);
        DenseVector<MemType, DataType, IndexType> y_rest(y, rest().rows(), first().rows());

        first().template apply<Algo_>(r_first, x_first, y_first, alpha);
        rest().template apply<Algo_>(r_rest, x_rest, y_rest, alpha);
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(first().create_vector_l(), rest().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(first().create_vector_r(), rest().create_vector_r());
      }

      template<typename Algo_>
      void scale_rows(const TupleDiagMatrix& a, const VectorTypeL& w)
      {
        first().template scale_rows<Algo_>(a.first(), w.first());
        rest().template scale_rows<Algo_>(a.rest(), w.rest());
      }

      template<typename Algo_>
      void scale_cols(const TupleDiagMatrix& a, const VectorTypeL& w)
      {
        first().template scale_cols<Algo_>(a.first(), w.first());
        rest().template scale_cols<Algo_>(a.rest(), w.rest());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        const Index brows(this->first().rows());

        if (row < brows)
        {
          return this->first().get_length_of_line(row);
        }
        else
        {
          return this->rest().get_length_of_line(row - brows);
        }
      }

      /// \cond internal
      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set,
                     const Index col_start, const Index stride = 1) const
      {
        const Index brows(this->first().rows());
        const Index bcolumns(this->first().columns());

        if (row < brows)
        {
          this->first().set_line(row, pval_set, pcol_set, col_start, stride);
        }
        else
        {
          this->rest().set_line(row - brows, pval_set, pcol_set, col_start + bcolumns, stride);
        }
      }
      /// \endcond

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename First2_, typename... Rest2_>
      void convert(const TupleDiagMatrix<First2_, Rest2_...>& other)
      {
        CONTEXT("When converting TupleDiagMatrix");

        this->first().convert(other.first());
        this->rest().convert(other.rest());
      }
    };

    /// \cond internal
    template<typename First_>
    class TupleDiagMatrix<First_>
    {
      template<typename,typename...>
      friend class TupleDiagMatrix;

    public:
      typedef typename First_::MemType MemType;
      typedef typename First_::DataType DataType;
      typedef typename First_::IndexType IndexType;
      /// Compatible L-vector type
      typedef TupleVector<typename First_::VectorTypeL> VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename First_::VectorTypeR> VectorTypeR;

#ifdef FEAST_COMPILER_MICROSOFT
      template<typename... RestVecL_>
      using VecLHelper = TupleVector<RestVecL_..., typename First_::VectorTypeL>;
      template<typename... RestVecR_>
      using VecRHelper = TupleVector<RestVecR_..., typename First_::VectorTypeR>;

      template <typename Mem2_, typename DT2_, typename IT2_, typename... RestCont_>
      using ContHelper = class TupleDiagMatrix<RestCont_..., typename First_::template ContainerType<Mem2_, DT2_, IT2_> >;
      template <typename Mem2_, typename DT2_, typename IT2_, typename... Dummy_>
      using ContainerType = class TupleDiagMatrix<typename First_::template ContainerType<Mem2_, DT2_, IT2_>, Dummy_...>;
#else
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class TupleDiagMatrix<typename First_::template ContainerType<Mem2_, DT2_, IT2_> >;
#endif

      static constexpr Index num_row_blocks = 1;
      static constexpr Index num_col_blocks = 1;

    protected:
      First_ _first;

      static String sub_name_list()
      {
        return First_::name();
      }

    public:
      /// default ctor
      TupleDiagMatrix()
      {
      }

      explicit TupleDiagMatrix(First_&& the_first) :
        _first(std::move(the_first))
      {
      }

      /// move ctor
      TupleDiagMatrix(TupleDiagMatrix&& other) :
        _first(std::move(other._first))
      {
      }

      /// move-assign operator
      TupleDiagMatrix& operator=(TupleDiagMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleDiagMatrix(const TupleDiagMatrix&) = delete;
      /// deleted copy-assign operator
      TupleDiagMatrix& operator=(const TupleDiagMatrix&) = delete;

      /// virtual destructor
      virtual ~TupleDiagMatrix()
      {
      }

      TupleDiagMatrix clone() const
      {
        return TupleDiagMatrix(_first.clone());
      }

      TupleDiagMatrix transpose() const
      {
        return TupleDiagMatrix(_first.transpose());
      }

      template<Index i, Index j>
      typename TupleElement<i, First_>::Type& at()
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _first;
      }

      template<Index i, Index j>
      const typename TupleElement<i, First_>::Type& at() const
      {
        static_assert(i == 0, "invalid sub-matrix index");
        static_assert(j == 0, "invalid sub-matrix index");
        return _first;
      }

      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
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
        return first().rows();
      }

      Index columns() const
      {
        return first().columns();
      }

      Index used_elements() const
      {
        return first().used_elements();
      }

      void format(DataType value = DataType(0))
      {
        first().format(value);
      }

      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().template apply<Algo_>(r.first(), x.first());
      }

      template<typename Algo_>
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().template apply<Algo_>(r.first(), x.first(), y.first(), alpha);
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(first().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(first().create_vector_r());
      }

      template<typename Algo_>
      void scale_rows(const TupleDiagMatrix& a, const VectorTypeL& w)
      {
        first().template scale_rows<Algo_>(a.first(), w.first());
      }

      template<typename Algo_>
      void scale_cols(const TupleDiagMatrix& a, const VectorTypeL& w)
      {
        first().template scale_cols<Algo_>(a.first(), w.first());
      }

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        return this->first().get_length_of_line(row);
      }

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        this->first().set_line(row, pval_set, pcol_set, col_start, stride);
      }

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
#ifdef FEAST_COMPILER_MICROSOFT
      template <typename First2_, typename... Rest2_>
      void convert(const TupleDiagMatrix<First2_, Rest2_...>& other)
      {
        static_assert(sizeof...(Rest2_) == std::size_t(0), "invalid TupleDiagMatrix size");
#else
      template <typename First2_>
      void convert(const TupleDiagMatrix<First2_>& other)
      {
#endif
        CONTEXT("When converting TupleDiagMatrix");

        this->first().convert(other.first());
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_DIAG_MATRIX_HPP
