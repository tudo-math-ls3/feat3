#pragma once
#ifndef KERNEL_LAFEM_TUPLE_VECTOR_HPP
#define KERNEL_LAFEM_TUPLE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/lafem/meta_element.hpp>

// includes, system
#include <iostream>
#include <type_traits>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Variadic TupleVector class template
     *
     * This class template implements a composition of sub-vectors of arbitrary classes.
     *
     * \tparam First_, ...Rest_
     * A sequence of (meta) vector classes which are to be composed.
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class TupleVector
    {
      template<typename,typename...>
      friend class TupleVector;

      typedef TupleVector<Rest_...> RestClass;

    public:
      /// number of vector blocks
      static constexpr Index num_blocks = TupleVector<Rest_...>::num_blocks + 1;


      /// sub-vector mem-type
      typedef typename First_::MemType MemType;
      /// sub-vector data-type
      typedef typename First_::DataType DataType;
      /// sub-vector index-type
      typedef typename First_::IndexType IndexType;

#ifdef FEAST_COMPILER_MICROSOFT
      /// helper class template for ContainerType template
      template <typename Mem2_, typename DT2_, typename IT2_, typename... RestCont_>
      using ContHelper = typename RestClass::template
        ContHelper<Mem2_, DT2_, IT2_, RestCont_..., typename First_::template ContainerType<Mem2_, DT2_, IT2_>>;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = ContHelper<Mem2_, DT2_, IT2_>;
#else
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleVector<
        typename First_::template ContainerType<Mem2_, DT2_, IT2_>,
        typename Rest_::template ContainerType<Mem2_, DT2_, IT2_>...>;
#endif

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value,
                    "sub-vectors have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value,
                    "sub-vectors have different data-types");
      static_assert(std::is_same<IndexType, typename RestClass::IndexType>::value,
                    "sub-vectors have different index-types");

    protected:
      /// the first sub-vector
      First_ _first;
      /// all remaining sub-vectors
      RestClass _rest;

      /// Returns a list of all sub-vector type names
      static String sub_name_list()
      {
        return First_::name() + "," + RestClass::sub_name_list();
      }

    public:
      /// default CTOR
      TupleVector()
      {
      }

      /// rest-class emplacement ctor; for internal use only
      explicit TupleVector(First_&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

      /// Sub-Vector emplacement constructor
      explicit TupleVector(First_&& the_first, Rest_&&... the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest...))
      {
      }

      /// move ctor
      TupleVector(TupleVector&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      TupleVector& operator=(TupleVector&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleVector(const TupleVector&) = delete;
      /// deleted copy-assign operator
      TupleVector& operator=(const TupleVector&) = delete;

      /**
       * \brief Creates and returns a copy of this vector
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      TupleVector clone(CloneMode clone_mode) const
      {
        return TupleVector(_first.clone(clone_mode), _rest.clone(clone_mode));
      }

      /**
       * \brief Creates and returns a copy of this vector
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be seperate.
       *
       */
      TupleVector clone() const
      {
        return TupleVector(_first.clone(), _rest.clone());
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      void clone(const TupleVector& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
        _rest.clone(other._rest, clone_mode);
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be seperate.
       *
       */
      void clone(const TupleVector& other)
      {
        _first.clone(other._first);
        _rest.clone(other._rest);
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

      TupleVector<Rest_...>& rest()
      {
        return _rest;
      }

      const TupleVector<Rest_...>& rest() const
      {
        return _rest;
      }
      /// \endcond

      /**
       * \brief Returns a sub-vector block.
       *
       * \tparam i_
       * The index of the sub-vector block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-vector at position \p i_.
       */
      template<Index i_>
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        static_assert(i_ < Index(num_blocks), "invalid sub-vector index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<Index i_>
      typename TupleElement<i_, First_, Rest_...>::Type const& at() const
      {
        static_assert(i_ < Index(num_blocks), "invalid sub-vector index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /// Returns the total size of this tuple-vector.
      Index size() const
      {
        return _first.size() + _rest.size();
      }

      /// Returns the number of blocks in this tuple-vector.
      Index blocks() const
      {
        return Index(num_blocks);
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("TupleVector<") + sub_name_list() + ">";
      }

      /// Clears the vector
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
      }

      /// Clears the vector
      void clear()
      {
        first().clear();
        rest().clear();
      }

      //template<typename First2_, typename... Rest2_>
      void copy(const TupleVector/*<First2_, Rest2_...>*/& x)
      {
        first().copy(x.first());
        rest().copy(x.rest());
      }

      void axpy(const TupleVector& x, const TupleVector& y, DataType alpha = DataType(1))
      {
        first().axpy(x.first(), y.first(), alpha);
        rest().axpy(x.rest(), y.rest(), alpha);
      }

      void component_product(const TupleVector & x, const TupleVector & y)
      {
        first().component_product(x.first(), y.first());
        rest().component_product(x.rest(), y.rest());
      }

      void component_invert(const TupleVector& x, DataType alpha = DataType(1))
      {
        first().component_invert(x.first(), alpha);
        rest().component_invert(x.rest(), alpha);
      }

      void scale(const TupleVector& x, DataType alpha)
      {
        first().scale(x.first(), alpha);
        rest().scale(x.rest(), alpha);
      }

      DataType dot(const TupleVector& x) const
      {
        return first().dot(x.first()) + rest().dot(x.rest());
      }

      /**
       * \copydoc LAFEM::DenseVector::triple_dot()
       **/
      DataType triple_dot(const TupleVector& x, const TupleVector& y) const
      {
        return first().triple_dot(x.first(), y.first())
          + rest().triple_dot(x.rest(), y.rest());
      }

      DataType norm2sqr() const
      {
        return first().norm2sqr() + rest().norm2sqr();
      }

      DataType norm2() const
      {
        return Math::sqrt(norm2sqr());
      }

      const typename First_::DataType operator()(Index index) const
      {
        CONTEXT("When retrieving TupleVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds tuple vector size " + stringify(size()) + " !");
        if (index < first().size())
        {
          return first()(index);
        }
        else
        {
          return rest()(index - first().size());
        }
      }

      void operator()(Index index, typename First_::DataType value)
      {
        CONTEXT("When retrieving TupleVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds tuple vector size " + stringify(size()) + " !");
        if (index < first().size())
        {
          first()(index, value);
        }
        else
        {
          rest()(index - first().size(), value);
        }
      }

      /// \cond internal
      /// Writes the vector-entries in an allocated array
      void set_vec(DataType * const pval_set) const
      {
        this->first().set_vec(pval_set);
        this->rest().set_vec(pval_set + this->first().size());
      }

      /// Writes data of an array in the vector
      void set_vec_inv(const DataType * const pval_set)
      {
        this->first().set_vec_inv(pval_set);
        this->rest().set_vec_inv(pval_set + this->first().size());
      }
      /// \endcond

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Vector.
       *
       * Use source vector content as content of current vector
       */
      template <typename First2_, typename... Rest2_>
      void convert(const TupleVector<First2_, Rest2_...>& other)
      {
        CONTEXT("When converting TupleVector");

        this->first().convert(other.first());
        this->rest().convert(other.rest());
      }
    }; // class TupleVector<...>

    /// \cond internal
    template<typename First_>
    class TupleVector<First_>
    {
      template<typename,typename...>
      friend class TupleVector;

    public:
      static constexpr Index num_blocks = 1;

      typedef typename First_::MemType MemType;
      typedef typename First_::DataType DataType;
      typedef typename First_::IndexType IndexType;

#ifdef FEAST_COMPILER_MICROSOFT
      template <typename Mem2_, typename DT2_, typename IT2_, typename... RestCont_>
      using ContHelper = class TupleVector<RestCont_..., typename First_::template ContainerType<Mem2_, DT2_, IT2_> >;
      template <typename Mem2_, typename DT2_, typename IT2_, typename... Dummy_>
      using ContainerType = class TupleVector<typename First_::template ContainerType<Mem2_, DT2_, IT2_>, Dummy_...>;
#else
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class TupleVector<typename First_::template ContainerType<Mem2_, DT2_, IT2_> >;
#endif

    protected:
      First_ _first;

      static String sub_name_list()
      {
        return First_::name();
      }

    public:
      TupleVector()
      {
      }

      explicit TupleVector(First_&& the_first) :
        _first(std::move(the_first))
      {
      }

      /// move-ctor
      TupleVector(TupleVector&& other) :
        _first(std::move(other._first))
      {
      }

      /// move-assign operator
      TupleVector& operator=(TupleVector&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleVector(const TupleVector&) = delete;
      /// deleted copy-assign operator
      TupleVector& operator=(const TupleVector&) = delete;

      /**
       * \brief Creates and returns a copy of this vector
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      TupleVector clone(CloneMode clone_mode) const
      {
        return TupleVector(_first.clone(clone_mode));
      }

      /**
       * \brief Creates and returns a copy of this vector
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be seperate.
       */
      TupleVector clone() const
      {
        return TupleVector(_first.clone());
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      void clone(const TupleVector& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be seperate.
       *
       */
      void clone(const TupleVector& other)
      {
        _first.clone(other._first);
      }

      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      Index size() const
      {
        return _first.size();
      }

      /// Returns the number of blocks in this tuple-vector.
      Index blocks() const
      {
        return Index(num_blocks);
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("TupleVector<") + sub_name_list() + ">";
      }

      template<Index i_>
      typename TupleElement<i_, First_>::Type& at()
      {
        static_assert(i_ == 0, "invalid sub-vector index");
        return first();
      }

      template<Index i_>
      typename TupleElement<i_, First_>::Type const& at() const
      {
        static_assert(i_ == 0, "invalid sub-vector index");
        return first();
      }

      void format(DataType value = DataType(0))
      {
        _first.format(value);
      }

      void clear()
      {
        _first.clear();
      }

      //template<typename First2_>
      void copy(const TupleVector/*<First2_>*/& x)
      {
        first().copy(x.first());
      }

      void axpy(const TupleVector& x, const TupleVector& y, DataType alpha = DataType(1))
      {
        first().axpy(x.first(), y.first(), alpha);
      }

      void component_product(const TupleVector & x, const TupleVector & y)
      {
        first().component_product(x.first(), y.first());
      }

      void component_invert(const TupleVector& x, DataType alpha = DataType(1))
      {
        first().component_invert(x.first(), alpha);
      }

      void scale(const TupleVector& x, DataType alpha)
      {
        first().scale(x.first(), alpha);
      }

      DataType dot(const TupleVector& x) const
      {
        return first().dot(x.first());
      }

      DataType triple_dot(const TupleVector& x, const TupleVector& y) const
      {
        return first().triple_dot(x.first(), y.first());
      }

      DataType triple_dot_i(const TupleVector& x, const TupleVector& y) const
      {
        return first().triple_dot_i(x.first(), y.first());
      }

      DataType norm2sqr() const
      {
        return first().norm2sqr();
      }

      DataType norm2() const
      {
        return first().norm2();
      }

      const typename First_::DataType operator()(Index index) const
      {
        CONTEXT("When retrieving TupleVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds tuple vector size " + stringify(size()) + " !");
        return first()(index);
      }

      void operator()(Index index, typename First_::DataType value)
      {
        CONTEXT("When retrieving TupleVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds tuple vector size " + stringify(size()) + " !");
        first()(index, value);
      }

      void set_vec(DataType * const pval_set) const
      {
        this->first().set_vec(pval_set);
      }

      void set_vec_inv(const DataType * const pval_set)
      {
        this->first().set_vec_inv(pval_set);
      }

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Vector.
       *
       * Use source vector content as content of current vector
       */
#ifdef FEAST_COMPILER_MICROSOFT
      template <typename First2_, typename... Rest2_>
      void convert(const TupleVector<First2_, Rest2_...>& other)
      {
        static_assert(sizeof...(Rest2_) == std::size_t(0), "invalud TupleVector size");
#else
      template <typename First2_>
      void convert(const TupleVector<First2_>& other)
      {
#endif
        CONTEXT("When converting TupleVector");

        this->first().convert(other.first());
      }
    };
    /// \endcond

    /// \cond internal
    template <typename First_>
    inline void dump_tuple_vector(std::ostream & os, const TupleVector<First_>& x)
    {
      os << x.first();
    }

    template <typename First_, typename... Rest_>
    inline void dump_tuple_vector(std::ostream & os, const TupleVector<First_, Rest_...>& x)
    {
      os << x.first() << ",";
      dump_tuple_vector<Rest_...>(os, x.rest());
    }

    template <typename First_, typename... Rest_>
    inline std::ostream& operator<< (std::ostream & os, const TupleVector<First_, Rest_...>& x)
    {
      os << "[";
      dump_tuple_vector(os, x);
      os << "]";
      return os;
    }
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_VECTOR_HPP
