#pragma once
#ifndef KERNEL_LAFEM_POWER_VECTOR_HPP
#define KERNEL_LAFEM_POWER_VECTOR_HPP 1

// includes, FEAST
#include <kernel/lafem/meta_element.hpp>

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
    class PowerVector
    {
      // declare this class template as a friend for recursive inheritance
      template<typename, Index>
      friend class PowerVector;

      /// base-class typedef
      typedef PowerVector<SubType_, count_ - 1> RestClass;

    public:
      /// sub-vector type
      typedef SubType_ SubVectorType;
      /// sub-vector memory type
      typedef typename SubVectorType::MemType MemType;
      /// sub-vector data type
      typedef typename SubVectorType::DataType DataType;
      /// sub-vector index type
      typedef typename SubVectorType::IndexType IndexType;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_, typename IT2_ = IndexType>
      using ContainerType = class PowerVector<typename SubType_::template ContainerType<Mem2_, DT2_, IT2_>, count_>;

      /// number of vector blocks
      enum
      {
        num_blocks = count_
      };

    protected:
      /// the first sub-vector
      SubVectorType _first;
      /// the remaining part
      RestClass _rest;

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
        _first(sub_size),
        _rest(sub_size)
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
        _first(sub_size, value),
        _rest(sub_size, value)
      {
      }

      /// base-class constructor; for internal use only
      explicit PowerVector(SubVectorType&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

      /// move CTOR
      PowerVector(PowerVector&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      PowerVector& operator=(PowerVector&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
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
        return PowerVector(_first.clone(), _rest.clone());
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
        return PowerElement<i_, SubVectorType>::get(*this);
      }

      /** \copydoc at() */
      template<Index i_>
      const SubVectorType& at() const
      {
        static_assert(i_ < count_, "invalid sub-vector index");
        return PowerElement<i_, SubVectorType>::get(*this);
      }

      /// \cond internal
      SubVectorType& first()
      {
        return _first;
      }

      const SubVectorType& first() const
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

      /// Returns the number of blocks in this power-vector.
      Index blocks() const
      {
        return Index(num_blocks);
      }

      /// Returns the total number of scalar entries of this power-vector.
      Index size() const
      {
        return first().size() + rest().size();
      }

      /**
       * \brief Clears this vector.
       *
       * \param[in] value
       * The value to which the vector's entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
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
        first().copy(x.first());
        rest().copy(x.rest());
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
        first().template axpy<Algo_>(x.first(), y.first(), alpha);
        rest().template axpy<Algo_>(x.rest(), y.rest(), alpha);
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
        first().template component_product<Algo_>(x.first(), y.first());
        rest().template component_product<Algo_>(x.rest(), y.rest());
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
        first().template component_product<Algo_>(x.first(), y.first(), z.first());
        rest().template component_product<Algo_>(x.rest(), y.rest(), z.rest());
      }

      /**
       * \brief Performs \f$ this_i \leftarrow \alpha / x_i \f$
       *
       * \param[in] x
       * The vector whose components serve as denominators.
       *
       * \param[in] alpha
       * The nominator.
       */
      template<typename Algo_>
      void component_invert(const PowerVector& x, DataType alpha = DataType(1))
      {
        first().template component_invert<Algo_>(x.first(), alpha);
        rest().template component_invert<Algo_>(x.rest(), alpha);
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
        first().template scale<Algo_>(x.first(), alpha);
        rest().template scale<Algo_>(x.rest(), alpha);
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
        return first().template dot<Algo_>(x.first()) + rest().template dot<Algo_>(x.rest());
      }

      /**
       * \copydoc LAFEM::DenseVector::triple_dot()
       **/
      template<typename Algo_>
      DataType triple_dot(const PowerVector& x, const PowerVector& y) const
      {
        return first().template triple_dot<Algo_>(x.first(), y.first())
          + rest().template triple_dot<Algo_>(x.rest(), y.rest());
      }

      /**
       * \brief Returns the squared euclid norm of this vector.
       */
      template<typename Algo_>
      DataType norm2sqr() const
      {
        return first().template norm2sqr<Algo_>() + rest().template norm2sqr<Algo_>();
      }

      /**
       * \brief Returns the euclid norm of this vector.
       */
      template<typename Algo_>
      DataType norm2() const
      {
        return Math::sqrt(norm2sqr<Algo_>());
      }

      /**
       * \brief Retrieve specific PowerVector element.
       *
       * \param[in] index The index of the vector element.
       *
       * \returns Specific vector element.
       */
      const DataType operator()(Index index) const
      {
        CONTEXT("When retrieving PowerVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds power vector size " + stringify(size()) + " !");

        if (index < first().size())
        {
          return first()(index);
        }
        else
        {
          return rest()(index - first().size());
        }
      }

      /**
       * \brief Set specific PowerVector element.
       *
       * \param[in] index The index of the vector element.
       * \param[in] value The value to be set.
       */
      void operator()(Index index, DataType value)
      {
        CONTEXT("When retrieving PowerVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds power vector size " + stringify(size()) + " !");

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
#ifdef FEAST_COMPILER_MICROSOFT
      template<typename SubType2_>
      void convert(const PowerVector<SubType2_, count_>& other)
#else
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const ContainerType<Mem2_, DT2_, IT2_> & other)
#endif
      {
        CONTEXT("When converting PowerVector");

        this->first().convert(other.first());
        this->rest().convert(other.rest());
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
      typedef typename SubVectorType::IndexType IndexType;
      template <typename Mem2_, typename DT2_, typename IT2_ = IndexType>
      using ContainerType = class PowerVector<typename SubType_::template ContainerType<Mem2_, DT2_, IT2_>, Index(1)>;

      enum
      {
        num_blocks = 1
      };

    protected:
      SubVectorType _first;

    public:
      PowerVector()
      {
      }

      explicit PowerVector(Index sub_size)
        : _first(sub_size)
      {
      }

      explicit PowerVector(Index sub_size, DataType value)
        : _first(sub_size, value)
      {
      }

      explicit PowerVector(SubVectorType&& the_first) :
        _first(std::move(the_first))
      {
      }

      PowerVector(PowerVector&& other) :
        _first(std::move(other._first))
      {
      }

      PowerVector& operator=(PowerVector&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
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
        return PowerVector(_first.clone());
      }

      template<Index i>
      SubVectorType& at()
      {
        static_assert(i == 0, "invalid sub-vector index");
        return _first;
      }

      template<Index i>
      const SubVectorType& at() const
      {
        static_assert(i == 0, "invalid sub-vector index");
        return _first;
      }

      SubVectorType& first()
      {
        return _first;
      }

      const SubVectorType& first() const
      {
        return _first;
      }

      Index blocks() const
      {
        return Index(1);
      }

      Index size() const
      {
        return _first.size();
      }

      void format(DataType value = DataType(0))
      {
        _first.format(value);
      }

      static String name()
      {
        return String("PowerVector<") + SubVectorType::name() + ",1>";
      }

      template<typename SubType2_>
      void copy(const PowerVector<SubType2_, 1>& x)
      {
        first().copy(x.first());
      }

      template<typename Algo_>
      void axpy(const PowerVector& x, const PowerVector& y, DataType alpha = DataType(1))
      {
        first().template axpy<Algo_>(x.first(), y.first(), alpha);
      }

      template <typename Algo_>
      void component_product(const PowerVector & x, const PowerVector & y)
      {
        first().template component_product<Algo_>(x.first(), y.first());
      }

      template <typename Algo_>
      void component_product(const PowerVector & x, const PowerVector & y, const PowerVector& z)
      {
        first().template component_product<Algo_>(x.first(), y.first(), z.first());
      }

      template<typename Algo_>
      void component_invert(const PowerVector& x, DataType alpha = DataType(1))
      {
        first().template component_invert<Algo_>(x.first(), alpha);
      }

      template<typename Algo_>
      void scale(const PowerVector& x, DataType alpha)
      {
        first().template scale<Algo_>(x.first(), alpha);
      }

      template<typename Algo_>
      DataType dot(const PowerVector& x) const
      {
        return first().template dot<Algo_>(x.first());
      }

      template<typename Algo_>
      DataType triple_dot(const PowerVector& x, const PowerVector& y) const
      {
        return first().template triple_dot<Algo_>(x.first(), y.first());
      }

      template<typename Algo_>
      DataType triple_dot_i(const PowerVector& x, const PowerVector& y) const
      {
        return first().template triple_dot_i<Algo_>(x.first(), y.first());
      }

      template<typename Algo_>
      DataType norm2sqr() const
      {
        return first().template norm2sqr<Algo_>();
      }

      template<typename Algo_>
      DataType norm2() const
      {
        return Math::sqrt(norm2sqr<Algo_>());
      }

      const DataType operator()(Index index) const
      {
        CONTEXT("When retrieving PowerVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds power vector size " + stringify(size()) + " !");

        return first()(index);
      }

      void operator()(Index index, DataType value)
      {
        CONTEXT("When retrieving PowerVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds power vector size " + stringify(size()) + " !");

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
      template<typename SubType2_>
      void convert(const PowerVector<SubType2_, Index(1)>& other)
#else
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const ContainerType<Mem2_, DT2_, IT2_> & other)
#endif
      {
        CONTEXT("When converting PowerVector");

        this->first().convert(other.first());
      }
    }; // class PowerVector<...,1>
    /// \cond internal

    /// \cond internal
    template <typename SubType_, Index count_>
    inline void dump_power_vector(std::ostream & os, const PowerVector<SubType_, count_>& x)
    {
      os << "," << x.first();
      dump_power_vector(os, x.rest());
    }

    template <typename SubType_>
    inline void dump_power_vector(std::ostream & os, const PowerVector<SubType_, 1>& x)
    {
      os << x.first();
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
