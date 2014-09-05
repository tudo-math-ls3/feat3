#pragma once
#ifndef KERNEL_LAFEM_TUPLE_MIRROR_HPP
#define KERNEL_LAFEM_TUPLE_MIRROR_HPP 1

#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/meta_element.hpp>

namespace FEAST
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
      /// sub-mirror index-type
      typedef typename First_::IndexType IndexType;

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value, "sub-mirrors have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value, "sub-mirrors have different data-types");
      static_assert(std::is_same<IndexType, typename RestClass::IndexType>::value, "sub-mirrors have different index-types");

    protected:
      /// the first sub-mirror
      First_ _first;
      /// all remaining sub-mirrors
      RestClass _rest;

      /// data-move ctor; this one is protected for a reason
      TupleMirror(First_&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

    public:
      /// default ctor
      TupleMirror()
      {
      }

      /// sub-mirror emplacement ctor
      explicit TupleMirror(First_&& the_first, Rest_&&... the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest...))
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

      /**
       * \brief Clone operation
       *
       * \returns A deep copy of this tuple mirror.
       */
      TupleMirror clone() const
      {
        return TupleMirror(_first.clone(), _rest.clone());
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
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      template<Index i_>
      typename TupleElement<i_, First_, Rest_...>::Type const& at() const
      {
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc VectorMirror::gather_prim() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Ty_, typename... Tv_>
      void gather_prim(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Ty_,Tv_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.template gather_prim<Algo_>(buffer, vector.first(), buffer_offset);
        _rest.template gather_prim<Algo_>(buffer, vector.rest(), buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::gather_axpy_prim() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Ty_, typename... Tv_>
      void gather_axpy_prim(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Ty_,Tv_...>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.template gather_axpy_prim<Algo_>(buffer, vector.first(), alpha, buffer_offset);
        _rest.template gather_axpy_prim<Algo_>(buffer, vector.rest(), alpha, buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::scatter_prim() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Ty_, typename... Tv_>
      void scatter_prim(
        LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        _first.template scatter_prim<Algo_>(vector.first(), buffer, buffer_offset);
        _rest.template scatter_prim<Algo_>(vector.rest(), buffer, buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::scatter_axpy_prim() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Ty_, typename... Tv_>
      void scatter_axpy_prim(
        LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.template scatter_axpy_prim<Algo_>(vector.first(), buffer, alpha, buffer_offset);
        _rest.template scatter_axpy_prim<Algo_>(vector.rest(), buffer, alpha, buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::gather_dual() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Ty_, typename... Tv_>
      void gather_dual(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.template gather_dual<Algo_>(buffer, vector.first(), buffer_offset);
        _rest.template gather_dual<Algo_>(buffer, vector.rest(), buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::gather_axpy_dual() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Ty_, typename... Tv_>
      void gather_axpy_dual(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.template gather_axpy_dual<Algo_>(buffer, vector.first(), alpha, buffer_offset);
        _rest.template gather_axpy_dual<Algo_>(buffer, vector.rest(), alpha, buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::scatter_dual() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Ty_, typename... Tv_>
      void scatter_dual(
        LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        _first.template scatter_dual<Algo_>(vector.first(), buffer, buffer_offset);
        _rest.template scatter_dual<Algo_>(vector.rest(), buffer, buffer_offset + _first.size());
      }

      /** \copydoc VectorMirror::scatter_axpy_dual() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Ty_, typename... Tv_>
      void scatter_axpy_dual(
        LAFEM::TupleVector<Ty_, Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.template scatter_axpy_dual<Algo_>(vector.first(), buffer, alpha, buffer_offset);
        _rest.template scatter_axpy_dual<Algo_>(vector.rest(), buffer, alpha, buffer_offset + _first.size());
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
      typedef typename First_::IndexType IndexType;

    protected:
      First_ _first;

    public:
      TupleMirror()
      {
      }

      explicit TupleMirror(First_&& the_first) :
        _first(std::move(the_first))
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

      TupleMirror clone() const
      {
        return TupleMirror(_first.clone());
      }

#ifdef FEAST_COMPILER_MICROSOFT
      template<typename... SubMirror2_>
      void convert(const TupleMirror<SubMirror2_...>& other)
      {
        static_assert(sizeof...(SubMirror2_) == std::size_t(1), "invalid TupleVector size");
        this->_first.convert(other._first);
      }
#else
      template<typename SubMirror2_>
      void convert(const TupleMirror<SubMirror2_>& other)
      {
        this->_first.convert(other._first);
      }
#endif

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

      template<Index i_>
      First_& at()
      {
        static_assert(i_ == 0, "invalid sub-mirror index");
        return _first;
      }

      template<Index i_>
      const First_& at() const
      {
        static_assert(i_ == 0, "invalid sub-mirror index");
        return _first;
      }

#ifdef FEAST_COMPILER_MICROSOFT
      // The MSVC compiler has a bug, which does not allow us to specify 'Tv_' as a single
      // template parameter rather than a parameter pack.
      // In consequence, we need to allow a whole pack here, and catch invalid sizes by
      // a static_assert within the function body...
      template<typename Algo_, typename Tx_, typename Ix_, typename... Tv_>
      void gather_prim(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Tv_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.template gather_prim<Algo_>(buffer, vector.first(), buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename... Tv_>
      void gather_axpy_prim(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Tv_...>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.template gather_axpy_prim<Algo_>(buffer, vector.first(), alpha, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename... Tv_>
      void scatter_prim(
        LAFEM::TupleVector<Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.template scatter_prim<Algo_>(vector.first(), buffer, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename... Tv_>
      void scatter_axpy_prim(
        LAFEM::TupleVector<Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.template scatter_axpy_prim<Algo_>(vector.first(), buffer, alpha, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename... Tv_>
      void gather_dual(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Tv_...>& vector,
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.template gather_dual<Algo_>(buffer, vector.first(), buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename... Tv_>
      void gather_axpy_dual(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Tv_...>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.template gather_axpy_dual<Algo_>(buffer, vector.first(), alpha, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename... Tv_>
      void scatter_dual(
        LAFEM::TupleVector<Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.template scatter_dual<Algo_>(vector.first(), buffer, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename... Tv_>
      void scatter_axpy_dual(
        LAFEM::TupleVector<Tv_...>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        static_assert(sizeof...(Tv_) == std::size_t(1), "invalid TupleVector size");
        _first.template scatter_axpy_dual<Algo_>(vector.first(), buffer, alpha, buffer_offset);
      }
#else // all other compilers
      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void gather_prim(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Tv_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.template gather_prim<Algo_>(buffer, vector.first(), buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void gather_axpy_prim(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Tv_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.template gather_axpy_prim<Algo_>(buffer, vector.first(), alpha, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void scatter_prim(
        LAFEM::TupleVector<Tv_>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        _first.template scatter_prim<Algo_>(vector.first(), buffer, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void scatter_axpy_prim(
        LAFEM::TupleVector<Tv_>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.template scatter_axpy_prim<Algo_>(vector.first(), buffer, alpha, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void gather_dual(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Tv_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        _first.template gather_dual<Algo_>(buffer, vector.first(), buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void gather_axpy_dual(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::TupleVector<Tv_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.template gather_axpy_dual<Algo_>(buffer, vector.first(), alpha, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void scatter_dual(
        LAFEM::TupleVector<Tv_>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        _first.template scatter_dual<Algo_>(vector.first(), buffer, buffer_offset);
      }

      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void scatter_axpy_dual(
        LAFEM::TupleVector<Tv_>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        _first.template scatter_axpy_dual<Algo_>(vector.first(), buffer, alpha, buffer_offset);
      }
#endif // FEAST_COMPILER_MICROSOFT
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUPLE_MIRROR_HPP
