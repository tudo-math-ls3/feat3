#pragma once
#ifndef KERNEL_GLOBAL_VECTOR_HPP
#define KERNEL_GLOBAL_VECTOR_HPP 1

#include <kernel/global/gate.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/container.hpp> // required for LAFEM::CloneMode

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Global vector wrapper class template
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename Mirror_>
    class Vector
    {
    public:
      typedef Gate<LocalVector_, Mirror_> GateType;

      typedef typename LocalVector_::MemType MemType;
      typedef typename LocalVector_::DataType DataType;
      typedef typename LocalVector_::IndexType IndexType;
      typedef LocalVector_ LocalVectorType;

      /// Our 'base' class type
      template <typename LocalVector2_, typename Mirror2_ = Mirror_>
      using ContainerType = class Vector<LocalVector2_, Mirror2_>;

      /// this typedef lets you create a vector container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = class Vector<typename LocalVector_::template ContainerType<Mem2_, DataType2_, IndexType2_>, typename Mirror_::template MirrorType<Mem2_, DataType2_, IndexType2_> >;

    protected:
      const GateType* _gate;
      LocalVector_ _vector;

    public:
      Vector() :
        _gate(nullptr),
        _vector()
      {
      }

      template<typename... Args_>
      explicit Vector(const GateType* gate, Args_&&... args) :
        _gate(gate),
        _vector(std::forward<Args_>(args)...)
      {
      }

      LocalVector_& operator*()
      {
        return _vector;
      }

      const LocalVector_& operator*() const
      {
        return _vector;
      }

      LocalVector_& local()
      {
        return _vector;
      }

      const LocalVector_& local() const
      {
        return _vector;
      }

      template<typename OtherGlobalVector_>
      void convert(const GateType* gate, const OtherGlobalVector_ & other)
      {
        this->_gate = gate;
        this->_vector.convert(*other);
      }

      const GateType* get_gate() const
      {
        return _gate;
      }

      const Dist::Comm* get_comm() const
      {
        return (_gate != nullptr ? _gate->get_comm() : nullptr);
      }

      void from_1_to_0()
      {
        if(_gate != nullptr)
        {
          _gate->from_1_to_0(this->_vector);
        }
      }

      void sync_0()
      {
        if(_gate != nullptr)
          _gate->sync_0(_vector);
      }

      auto sync_0_async() -> decltype(_gate->sync_0_async(_vector))
      //decltype(_gate->sync_0(_vector)) sync_0_async()
      {
        return _gate->sync_0_async(_vector);
      }

      void sync_1()
      {
        if(_gate != nullptr)
          _gate->sync_1(_vector);
      }

      auto sync_1_async() -> decltype(_gate->sync_1_async(_vector))
      //decltype(_gate->sync_1(_vector)) sync_1_async()
      {
        return _gate->sync_1_async(_vector);
      }

      /**
       * \brief Gets the total number of elements in this vector
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of elements
       */
      Index rows() const
      {
        // Compute total number of rows
        auto vec_l = clone();
        vec_l.format(DataType(1));

        return Index(vec_l.norm2sqr());
      }

      Vector clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return Vector(_gate, _vector.clone(mode));
      }

      void clone(const Vector& other, LAFEM::CloneMode mode = LAFEM::CloneMode::Weak)
      {
        if(&(*other) == &(*(*this)))
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Trying to self-clone a Global::Vector!");
        }

        *(*this) = (*other).clone(mode);
      }

      void clear()
      {
        _vector.clear();
      }

      /**
       * \brief Reset all elements of the container to a given value or zero if missing.
       *
       * \param[in] alpha The value to be set (defaults to 0)
       */
      void format(DataType alpha = DataType(0))
      {
        _vector.format(alpha);
      }

      /**
       * \brief Reset all elements of the container to random values.
       *
       * \param[in] rng The random number generator.
       * \param[in] min Lower rng bound.
       * \param[in] max Upper rng bound.
       */
      void format(Random & rng, DataType min, DataType max)
      {
        _vector.format(rng, min, max);
        sync_1();
      }

      void copy(const Vector& x)
      {
        // avoid self-copy
        if(this != &x)
        {
          _vector.copy(*x);
        }
      }

      void axpy(const Vector& x, const Vector& y, const DataType alpha = DataType(1))
      {
        _vector.axpy(*x, *y, alpha);
      }

      void scale(const Vector& x, const DataType alpha)
      {
        _vector.scale(*x, alpha);
      }

      DataType dot(const Vector& x) const
      {
        if(_gate != nullptr)
          return _gate->dot(_vector, *x);
        return _vector.dot(*x);
      }

      std::shared_ptr<SynchScalarTicket<DataType>> dot_async(const Vector& x) const
      {
        return _gate->dot_async(_vector, *x);
      }

      DataType norm2sqr() const
      {
        return dot(*this);
      }

      std::shared_ptr<SynchScalarTicket<DataType>> norm2sqr_async() const
      {
        return dot_async(*this);
      }

      DataType norm2() const
      {
        return Math::sqrt(norm2sqr());
      }

      std::shared_ptr<SynchScalarTicket<DataType>> norm2_async() const
      {
        return _gate->dot_async(_vector, _vector, true);
      }

      void component_invert(const Vector& x, const DataType alpha = DataType(1))
      {
        _vector.component_invert(*x, alpha);
      }

      void component_product(const Vector& x, const Vector& y)
      {
        _vector.component_product(*x, *y);
      }

      DataType max_element() const
      {
        if (_gate != nullptr)
          return _gate->max_element(_vector);
        return _vector.max_element();
      }

      std::shared_ptr<SynchScalarTicket<DataType>> max_element_async() const
      {
        return _gate->max_element_async(_vector);
      }
    };
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_VECTOR_HPP
