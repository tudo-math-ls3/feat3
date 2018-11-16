#pragma once
#ifndef KERNEL_LAFEM_TUPLE_MIRROR_HPP
#define KERNEL_LAFEM_TUPLE_MIRROR_HPP 1

#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/meta_element.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief TupleVector meta-mirror class template
     *
     * This class template implements a vector-mirror which is capable of mirroring
     * TupleVector objects.
     *
     * \tparam First_, Rest_...
     * A sequence of (meta) vector mirrors which are to be composed.
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class TupleMirror
    {
      template<typename,typename...>
      friend class TupleMirror;

      typedef TupleMirror<Rest_...> RestClass;

    public:
      /// number of mirror blocks
      static constexpr int num_blocks = RestClass::num_blocks + 1;

      /// sub-mirror mem-type
      typedef typename First_::MemType MemType;
      /// sub-mirror data-type
      typedef typename First_::DataType DataType;
      /// sub-mirror index-type
      typedef typename First_::IndexType IndexType;

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value, "sub-mirrors have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value, "sub-mirrors have different data-types");
      static_assert(std::is_same<IndexType, typename RestClass::IndexType>::value, "sub-mirrors have different index-types");

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using MirrorType = TupleMirror<
        typename First_::template MirrorType<Mem2_, DT2_, IT2_>,
        typename Rest_::template MirrorType<Mem2_, DT2_, IT2_>...>;

      /// this typedef lets you create a mirror with new Memory, Data and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using MirrorTypeByMDI = MirrorType<Mem2_, DataType2_, IndexType2_>;

    protected:
      /// the first sub-mirror
      First_ _first;
      /// all remaining sub-mirrors
      RestClass _rest;

      /// data-move ctor; this one is protected for a reason
      TupleMirror(First_&& the_first, RestClass&& the_rest) :
        _first(std::forward<First_>(the_first)),
        _rest(std::forward<RestClass>(the_rest))
      {
      }

    public:
      /// default ctor
      TupleMirror()
      {
      }

      /// sub-mirror emplacement ctor
      explicit TupleMirror(First_&& the_first, Rest_&&... the_rest) :
        _first(std::forward<First_>(the_first)),
        _rest(std::forward<Rest_>(the_rest)...)
      {
      }

      /// move-ctor
      TupleMirror(TupleMirror&& other) :
        _first(std::forward<First_>(other._first)),
        _rest(std::forward<RestClass>(other._rest))
      {
      }

      /// move-assign operator
      TupleMirror& operator=(TupleMirror&& other)
      {
        if(this != &other)
        {
          _first = std::forward<First_>(other._first);
          _rest = std::forward<RestClass>(other._rest);
        }
        return *this;
      }

      /**
       * \brief Clone operation
       *
       * \returns A deep copy of this tuple mirror.
       */
      TupleMirror clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        return TupleMirror(_first.clone(clone_mode), _rest.clone(clone_mode));
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other
       * The source mirror.
       */
      template<typename... SubMirror2_>
      void convert(const TupleMirror<SubMirror2_...>& other)
      {
        this->_first.convert(other._first);
        this->_rest.convert(other._rest);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
      }

      /**
       * \brief Checks whether the mirror is empty.
       *
       * \returns \c true, if there are no indices in the mirror, otherwise \c false.
       */
      bool empty() const
      {
        return _first.empty() && _rest.empty();
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
      /// \endcond

      /**
       * \brief Computes the required buffer size for a TupleVector.
       *
       * \tparam[in] vector
       * The vector whose buffer size is to be computed.
       */
      template<typename Tv_, typename... Tw_>
      Index buffer_size(const TupleVector<Tv_, Tw_...>& vector) const
      {
        return _first.buffer_size(vector.first()) + _rest.buffer_size(vector.rest());
      }

      /**
       * \brief Creates a new buffer vector for a TupleVector.
       *
       * \tparam[in] vector
       * The vector for which the buffer is to be created.
       */
      template<typename Tv_, typename... Tw_>
      DenseVector<MemType, DataType, IndexType> create_buffer(const TupleVector<Tv_, Tw_...>& vector) const
      {
        return DenseVector<MemType, DataType, IndexType>(buffer_size(vector), Pinning::disabled);
      }

      /**
       * \brief Returns a sub-mirror block.
       *
       * \tparam i_
       * The index of the sub-mirror block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-mirror at position \p i_.
       */
      template<int i_>
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<int i_>
      typename TupleElement<i_, First_, Rest_...>::Type const& at() const
      {
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc VectorMirror::gather() */
      template<typename Tv_, typename... Tw_>
      void gather(
        LAFEM::DenseVector<MemType, DataType, IndexType>& buffer,
        const LAFEM::TupleVector<Tv_, Tw_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.gather(buffer, vector.first(), buffer_offset);
        _rest.gather(buffer, vector.rest(), buffer_offset + _first.buffer_size(vector.first()));
      }

      /** \copydoc VectorMirror::scatter_axpy() */
      template<typename Tv_, typename... Tw_>
      void scatter_axpy(
        LAFEM::TupleVector<Tv_, Tw_...>& vector,
        const LAFEM::DenseVector<MemType, DataType, IndexType>& buffer,
        const DataType alpha = DataType(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.scatter_axpy(vector.first(), buffer, alpha, buffer_offset);
        _rest.scatter_axpy(vector.rest(), buffer, alpha, buffer_offset + _first.buffer_size(vector.first()));
      }
    }; // class TupleMirror<...>

    /// \cond internal
    template<typename First_>
    class TupleMirror<First_>
    {
      template<typename,typename...>
      friend class TupleMirror;

    public:
      /// number of mirror blocks
      static constexpr int num_blocks = 1;

      typedef typename First_::MemType MemType;
      typedef typename First_::DataType DataType;
      typedef typename First_::IndexType IndexType;

      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using MirrorType = class TupleMirror<typename First_::template MirrorType<Mem2_, DT2_, IT2_> >;

      /// this typedef lets you create a mirror with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using MirrorTypeByMDI = MirrorType<Mem2_, DataType2_, IndexType2_>;

    protected:
      First_ _first;

    public:
      TupleMirror()
      {
      }

      explicit TupleMirror(First_&& the_first) :
        _first(std::forward<First_>(the_first))
      {
      }

      TupleMirror(TupleMirror&& other) :
        _first(std::forward<First_>(other._first))
      {
      }

      TupleMirror& operator=(TupleMirror&& other)
      {
        if(this != &other)
        {
          _first = std::forward<First_>(other._first);
        }
        return *this;
      }

      TupleMirror clone() const
      {
        return TupleMirror(_first.clone());
      }

      template<typename SubMirror2_>
      void convert(const TupleMirror<SubMirror2_>& other)
      {
        this->_first.convert(other._first);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes();
      }
      bool empty() const
      {
        return _first.empty();
      }

      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      template<typename Tv_>
      Index buffer_size(const TupleVector<Tv_>& vector) const
      {
        return _first.buffer_size(vector.first());
      }

      template<typename Tv_>
      DenseVector<MemType, DataType, IndexType> create_buffer(const TupleVector<Tv_>& vector) const
      {
        return DenseVector<MemType, DataType, IndexType>(buffer_size(vector), Pinning::disabled);
      }

      template<int i_>
      First_& at()
      {
        static_assert(i_ == 0, "invalid sub-mirror index");
        return _first;
      }

      template<int i_>
      const First_& at() const
      {
        static_assert(i_ == 0, "invalid sub-mirror index");
        return _first;
      }

      template<typename Tv_>
      void gather(
        LAFEM::DenseVector<MemType, DataType, IndexType>& buffer,
        const LAFEM::TupleVector<Tv_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.gather(buffer, vector.first(), buffer_offset);
      }

      template<typename Tv_>
      void scatter_axpy(
        LAFEM::TupleVector<Tv_>& vector,
        const LAFEM::DenseVector<MemType, DataType, IndexType>& buffer,
        const DataType alpha = DataType(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.scatter_axpy(vector.first(), buffer, alpha, buffer_offset);
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_TUPLE_MIRROR_HPP
