#pragma once
#ifndef KERNEL_LAFEM_POWER_VECTOR_HPP
#define KERNEL_LAFEM_POWER_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// includes, system
#include <iostream>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Power-Vector meta class template
     *
     * This class template implements a composition of \e n sub-vectors of the same class.
     *
     * \note For a composed meta-vector of different sub-vector classes, see the TupleVector class template.
     *
     * \tparam SubType_
     * The type of the sub-vector.
     *
     * \tparam count_
     * The number of sub-vector blocks.
     *
     * \author Peter Zajac
     */
    template<
      typename SubType_,
      Index count_>
    class PowerVector :
      protected PowerVector<SubType_, count_ - 1>
    {
      // declare this class template as a friend for recursive inheritance
      template<typename, Index>
      friend class PowerVector;

      /// base-class typedef
      typedef PowerVector<SubType_, count_ - 1> BaseClass;

    public:
      /// sub-vector type
      typedef SubType_ SubVectorType;
      /// sub-vector memory type
      typedef typename SubVectorType::MemType MemType;
      /// sub-vector data type
      typedef typename SubVectorType::DataType DataType;

      /// dummy enum
      enum
      {
        /// number of vector blocks
        num_blocks = count_
      };

    protected:
      /// the last sub-vector
      SubVectorType _sub_vector;

    public:
      /// default CTOR
      PowerVector()
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] sub_size
       * The size of any sub-vector of this power-vector.
       */
      explicit PowerVector(Index sub_size) :
        BaseClass(sub_size),
        _sub_vector(sub_size)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] sub_size
       * The size of any sub-vector of this power-vector.
       *
       * \param[in] value
       * The value that the sub-vectors are to be initialised with.
       */
      explicit PowerVector(Index sub_size, DataType value) :
        BaseClass(sub_size, value),
        _sub_vector(sub_size, value)
      {
      }

      /// base-class constructor; for internal use only
      explicit PowerVector(BaseClass&& other_base, SubVectorType&& last_sub) :
        BaseClass(std::move(other_base)),
        _sub_vector(std::move(last_sub))
      {
      }

      /// move CTOR
      PowerVector(PowerVector&& other) :
        BaseClass(static_cast<BaseClass&&>(other)),
        _sub_vector(std::move(other._sub_vector))
      {
      }

      /// move-assign operator
      PowerVector& operator=(PowerVector&& other)
      {
        if(this != &other)
        {
          base().operator=(static_cast<BaseClass&&>(other));
          _sub_vector = std::move(other._sub_vector);
        }
        return *this;
      }

      /// deleted copy-ctor
      PowerVector(const PowerVector&) = delete;
      /// deleted copy-assign operator
      PowerVector& operator=(const PowerVector&) = delete;

      /// virtual destructor
      virtual ~PowerVector()
      {
      }

      /**
       * \brief Creates and returns a deep copy of this vector.
       */
      PowerVector clone() const
      {
        return PowerVector(base().clone(), _sub_vector.clone());
      }

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
      SubVectorType& at()
      {
        static_assert(i_ < count_, "invalid sub-vector index");
        return static_cast<PowerVector<SubType_, i_+1>&>(*this)._sub_vector;
      }

      /** \copydoc at() */
      template<Index i_>
      const SubVectorType& at() const
      {
        static_assert(i_ < count_, "invalid sub-vector index");
        return static_cast<const PowerVector<SubType_, i_+1>&>(*this)._sub_vector;
      }

      /// \cond internal
      SubVectorType& last()
      {
        return _sub_vector;
      }

      const SubVectorType& last() const
      {
        return _sub_vector;
      }

      PowerVector<SubType_, count_ - 1>& base()
      {
        return static_cast<BaseClass&>(*this);
      }

      const PowerVector<SubType_, count_ - 1>& base() const
      {
        return static_cast<const BaseClass&>(*this);
      }
      /// \endcond

      /// Returns the number of blocks in this power-vector.
      Index blocks() const
      {
        return Index(num_blocks);
      }

      /// Returns the total number of scalar entries of this power-vector.
      Index size() const
      {
        return base().size() + last().size();
      }

      /**
       * \brief Clears this vector.
       *
       * \param[in] value
       * The value to which the vector's entries are to be set to.
       */
      void clear(DataType value = DataType(0))
      {
        base().clear(value);
        last().clear(value);
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("PowerVector<") + SubVectorType::name() + "," + stringify(num_blocks) + ">";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The vector to be copied.
       */
      template<typename SubType2_>
      void copy(const PowerVector<SubType2_, count_>& x)
      {
        base().copy(x.base());
        last().copy(x.last());
      }

      /**
       * \brief Performs \f$ this \leftarrow \alpha\cdot x + y \f$
       *
       * \param[in] x
       * The first summand vector.
       *
       * \param[in] y
       * The second summand vector.
       *
       * \param[in] alpha
       * The scaling factor for \p x.
       */
      template<typename Algo_>
      void axpy(const PowerVector& x, const PowerVector& y, DataType alpha = DataType(1))
      {
        base().template axpy<Algo_>(x.base(), y.base(), alpha);
        last().template axpy<Algo_>(x.last(), y.last(), alpha);
      }

      /**
       * \brief Performs \f$this_i \leftarrow x_i \cdot y_i\f$
       *
       * \param[in] x The first factor.
       * \param[in] y The second factor.
       */
      template <typename Algo_>
      void component_product(const PowerVector & x, const PowerVector & y)
      {
        base().template component_product<Algo_>(x.base(), y.base());
        last().template component_product<Algo_>(x.last(), y.last());
      }

      /**
       * \brief Performs \f$this_i \leftarrow x_i \cdot y_i + z_i\f$
       *
       * \param[in] x The first factor.
       * \param[in] y The second factor.
       * \param[in] z The second summand.
       */
      template <typename Algo_>
      void component_product(const PowerVector & x, const PowerVector & y, const PowerVector& z)
      {
        base().template component_product<Algo_>(x.base(), y.base(), z.base());
        last().template component_product<Algo_>(x.last(), y.last(), z.last());
      }

      /**
       * \brief Performs \f$ this \leftarrow \alpha\cdot x \f$
       *
       * \param[in] x
       * The vector to be scaled.
       *
       * \param[in] alpha
       * The scaling factor for \p x.
       */
      template<typename Algo_>
      void scale(const PowerVector& x, DataType alpha)
      {
        base().template scale<Algo_>(x.base(), alpha);
        last().template scale<Algo_>(x.last(), alpha);
      }

      /**
       * \brief Returns the dot-product of this vector and another vector.
       *
       * \param[in] x
       * The second vector for the dot-product.
       */
      template<typename Algo_>
      DataType dot(const PowerVector& x) const
      {
        return base().template dot<Algo_>(x.base()) + last().template dot<Algo_>(x.last());
      }

      /**
       * \brief Returns the squared euclid norm of this vector.
       */
      template<typename Algo_>
      DataType norm2sqr() const
      {
        return base().template norm2sqr<Algo_>() + last().template norm2sqr<Algo_>();
      }

      /**
       * \brief Returns the euclid norm of this vector.
       */
      template<typename Algo_>
      DataType norm2() const
      {
        return Math::sqrt(norm2sqr<Algo_>());
      }
    }; // class PowerVector<...>

    /// \cond internal
    /**
     * \brief Specialisation of PowerVector for only one block
     *
     * \author Peter Zajac
     */
    template<typename SubType_>
    class PowerVector<SubType_, 1>
    {
      template<typename, Index>
      friend class PowerVector;

    public:
      typedef SubType_ SubVectorType;

      typedef typename SubVectorType::MemType MemType;
      typedef typename SubVectorType::DataType DataType;

      enum
      {
        num_blocks = 1
      };

    protected:
      SubVectorType _sub_vector;

    public:
      PowerVector()
      {
      }

      explicit PowerVector(Index sub_size)
        : _sub_vector(sub_size)
      {
      }

      explicit PowerVector(Index sub_size, DataType value)
        : _sub_vector(sub_size, value)
      {
      }

      explicit PowerVector(SubVectorType&& sub_vector) :
        _sub_vector(sub_vector)
      {
      }

      PowerVector(PowerVector&& other) :
        _sub_vector(std::move(other._sub_vector))
      {
      }

      PowerVector& operator=(PowerVector&& other)
      {
        if(this != &other)
        {
          _sub_vector = std::move(other._sub_vector);
        }
        return *this;
      }

      /// deleted copy-ctor
      PowerVector(const PowerVector&) = delete;
      /// deleted copy-assign operator
      PowerVector& operator=(const PowerVector&) = delete;

      virtual ~PowerVector()
      {
      }

      PowerVector clone() const
      {
        return PowerVector(_sub_vector.clone());
      }

      template<Index i>
      SubVectorType& at()
      {
        static_assert(i == 0, "invalid sub-vector index");
        return _sub_vector;
      }

      template<Index i>
      const SubVectorType& at() const
      {
        static_assert(i == 0, "invalid sub-vector index");
        return _sub_vector;
      }

      SubVectorType& last()
      {
        return _sub_vector;
      }

      const SubVectorType& last() const
      {
        return _sub_vector;
      }

      Index blocks() const
      {
        return Index(1);
      }

      Index size() const
      {
        return _sub_vector.size();
      }

      void clear(DataType value = DataType(0))
      {
        _sub_vector.clear(value);
      }

      static String name()
      {
        return String("PowerVector<") + SubVectorType::name() + ",1>";
      }

      template<typename SubType2_>
      void copy(const PowerVector<SubType2_, 1>& x)
      {
        last().copy(x.last());
      }

      template<typename Algo_>
      void axpy(const PowerVector& x, const PowerVector& y, DataType alpha = DataType(1))
      {
        last().template axpy<Algo_>(x.last(), y.last(), alpha);
      }

      template <typename Algo_>
      void component_product(const PowerVector & x, const PowerVector & y)
      {
        last().template component_product<Algo_>(x.last(), y.last());
      }

      template <typename Algo_>
      void component_product(const PowerVector & x, const PowerVector & y, const PowerVector& z)
      {
        last().template component_product<Algo_>(x.last(), y.last(), z.last());
      }

      template<typename Algo_>
      void scale(const PowerVector& x, DataType alpha)
      {
        last().template scale<Algo_>(x.last(), alpha);
      }

      template<typename Algo_>
      DataType dot(const PowerVector& x) const
      {
        return last().template dot<Algo_>(x.last());
      }

      template<typename Algo_>
      DataType norm2sqr() const
      {
        return last().template norm2sqr<Algo_>();
      }

      template<typename Algo_>
      DataType norm2() const
      {
        return Math::sqrt(norm2sqr<Algo_>());
      }
    }; // class PowerVector<...,1>
    /// \cond internal

    /// \cond internal
    template <typename SubType_, Index count_>
    inline void dump_power_vector(std::ostream & os, const PowerVector<SubType_, count_>& x)
    {
      dump_power_vector(os, x.base());
      os << "," << x.last();
    }

    template <typename SubType_>
    inline void dump_power_vector(std::ostream & os, const PowerVector<SubType_, 1>& x)
    {
      os << x.last();
    }

    template <typename SubType_, Index count_>
    inline std::ostream& operator<< (std::ostream & os, const PowerVector<SubType_, count_>& x)
    {
      os << "[";
      dump_power_vector(os, x);
      os << "]";
      return os;
    }
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_VECTOR_HPP
