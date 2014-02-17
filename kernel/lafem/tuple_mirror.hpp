#pragma once
#ifndef KERNEL_LAFEM_TUPLE_MIRROR_HPP
#define KERNEL_LAFEM_TUPLE_MIRROR_HPP 1

#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    template<
      Index i_,
      typename First_,
      typename... Rest_>
    struct TupleMirrorElement;

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
      /// dummy enum
      enum
      {
        /// number of mirror blocks
        num_blocks = RestClass::num_blocks + 1
      };

      /// sub-mirror mem-type
      typedef typename First_::MemType MemType;
      /// sub-mirror data-type
      typedef typename First_::DataType DataType;

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value, "sub-mirrors have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value, "sub-mirrors have different data-types");

    protected:
      /// the first sub-mirror
      First_ _first;
      /// all remaining sub-mirrors
      RestClass _rest;

      /// data-move ctor; this one is protected for a reason
      TupleMirror(First_&& first, RestClass&& rest) :
        _first(std::move(first)),
        _rest(std::move(rest))
      {
      }

    public:
      /// default ctor
      TupleMirror()
      {
      }

      /// sub-mirror emplacement ctor
      explicit TupleMirror(First_&& first, Rest_&&... rest) :
        _first(std::move(first)),
        _rest(std::move(rest...))
      {
      }

      /// move-ctor
      TupleMirror(TupleMirror&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      TupleMirror& operator=(TupleMirror&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
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

      /// Returns the total size of the mirror.
      Index size() const
      {
        return _first.size() + _rest.size();
      }

      template<Index i_>
      typename TupleMirrorElement<i_, First_, Rest_...>::Type& at()
      {
        return TupleMirrorElement<i_, First_, Rest_...>::get(*this);
      }

      template<Index i_>
      typename TupleMirrorElement<i_, First_, Rest_...>::Type const& at() const
      {
        return TupleMirrorElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc VectorMirror::gather_prim() */
      template<typename Tx_, typename Ty_, typename... Tv_>
      void gather_prim(
        LAFEM::DenseVector<MemType, Tx_>& buffer,
        const LAFEM::TupleVector<Ty_,Tv_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.gather_prim(buffer, vector.first(), buffer_offset);
        _rest.gather_prim(buffer, vector.rest(), buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::scatter_prim() */
      template<typename Tx_, typename Ty_, typename... Tv_>
      void scatter_prim(
        LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        _first.scatter_prim(vector.first(), buffer, buffer_offset);
        _rest.scatter_prim(vector.rest(), buffer, buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::gather_dual() */
      template<typename Tx_, typename Ty_, typename... Tv_>
      void gather_dual(
        LAFEM::DenseVector<MemType, Tx_>& buffer,
        const LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.gather_dual(buffer, vector.first(), buffer_offset);
        _rest.gather_dual(buffer, vector.rest(), buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::scatter_dual() */
      template<typename Tx_, typename Ty_, typename... Tv_>
      void scatter_dual(
        LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        _first.scatter_dual(vector.first(), buffer, buffer_offset);
        _rest.scatter_dual(vector.rest(), buffer, buffer_offset + _first.size());
      }
    }; // class TupleMirror<...>

    /// \cond internal
    template<typename First_>
    class TupleMirror<First_>
    {
      template<typename,typename...>
      friend class TupleMirror;

    public:
      /// dummy enum
      enum
      {
        /// number of mirror blocks
        num_blocks = 1
      };

      typedef typename First_::MemType MemType;
      typedef typename First_::DataType DataType;

    protected:
      First_ _first;

    public:
      TupleMirror()
      {
      }

      explicit TupleMirror(First_&& first) :
        _first(std::move(first))
      {
      }

      TupleMirror(TupleMirror&& other) :
        _first(std::move(other._first))
      {
      }

      TupleMirror& operator=(TupleMirror&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
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
      /// \endcond

      Index size() const
      {
        return _first.size();
      }

#ifdef FEAST_COMPILER_MICROSOFT
      // The MSVC compiler has a bug, which does not allow us to specify 'Tv_' as a single
      // template parameter rather than a parameter pack.
      // In consequence, we need to allow a whole pack here, and catch invalid sizes by
      // a static_assert within the function body...
      template<typename Tx_, typename... Tv_>
      void gather_prim(
        LAFEM::DenseVector<MemType, Tx_>& buffer,
        const LAFEM::TupleVector<Tv_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.gather_prim(buffer, vector.first(), buffer_offset);
      }

      template<typename Tx_, typename... Tv_>
      void scatter_prim(
        LAFEM::TupleVector<Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.scatter_prim(vector.first(), buffer, buffer_offset);
      }

      template<typename Tx_, typename... Tv_>
      void gather_dual(
        LAFEM::DenseVector<MemType, Tx_>& buffer,
        const LAFEM::TupleVector<Tv_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.gather_dual(buffer, vector.first(), buffer_offset);
      }

      template<typename Tx_, typename... Tv_>
      void scatter_dual(
        LAFEM::TupleVector<Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.scatter_dual(vector.first(), buffer, buffer_offset);
      }
#else // all other compilers
      template<typename Tx_, typename Tv_>
      void gather_prim(
        LAFEM::DenseVector<MemType, Tx_>& buffer,
        const LAFEM::TupleVector<Tv_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.gather_prim(buffer, vector.first(), buffer_offset);
      }

      template<typename Tx_, typename Tv_>
      void scatter_prim(
        LAFEM::TupleVector<Tv_>& vector,
        const LAFEM::DenseVector<MemType, Tx_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        _first.scatter_prim(vector.first(), buffer, buffer_offset);
      }

      template<typename Tx_, typename Tv_>
      void gather_dual(
        LAFEM::DenseVector<MemType, Tx_>& buffer,
        const LAFEM::TupleVector<Tv_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.gather_dual(buffer, vector.first(), buffer_offset);
      }

      template<typename Tx_, typename Tv_>
      void scatter_dual(
        LAFEM::TupleVector<Tv_>& vector,
        const LAFEM::DenseVector<MemType, Tx_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        _first.scatter_dual(vector.first(), buffer, buffer_offset);
      }
#endif // FEAST_COMPILER_MICROSOFT
    };
    /// \endcond

    template<
      Index i_,
      typename First_,
      typename... Rest_>
    struct TupleMirrorElement
    {
      typedef typename TupleMirrorElement<i_-1, Rest_...>::Type Type;

      static Type& get(TupleMirror<First_, Rest_...>& meta)
      {
        return TupleMirrorElement<i_-1, Rest_...>::get(meta.rest());
      }

      static const Type& get(const TupleMirror<First_, Rest_...>& meta)
      {
        return TupleMirrorElement<i_-1, Rest_...>::get(meta.rest());
      }
    };

    template<
      typename First_,
      typename... Rest_>
    struct TupleMirrorElement<0, First_, Rest_...>
    {
      typedef First_ Type;

      static Type& get(TupleMirror<First_, Rest_...>& meta)
      {
        return meta.first();
      }

      static const Type& get(const TupleMirror<First_, Rest_...>& meta)
      {
        return meta.first();
      }
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_MIRROR_HPP
