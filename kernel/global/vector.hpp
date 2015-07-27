#pragma once
#ifndef KERNEL_GLOBAL_VECTOR_HPP
#define KERNEL_GLOBAL_VECTOR_HPP 1

#include <kernel/global/gate.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp> // required for LAFEM::CloneMode

namespace FEAST
{
  namespace Global
  {
    /**
     * \brief Global vector wrapper class template
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_>
    class Vector
    {
    public:
      typedef Gate<LocalVector_> GateType;

      typedef typename LocalVector_::MemType MemType;
      typedef typename LocalVector_::DataType DataType;
      typedef typename LocalVector_::IndexType IndexType;

    protected:
      GateType* _gate;
      LocalVector_ _vector;

    public:
      Vector() :
        _gate(nullptr),
        _vector()
      {
      }

      template<typename... Args_>
      explicit Vector(GateType* gate, Args_&&... args) :
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

      template<typename OtherLocalVector_>
      void convert(GateType* gate, const Global::Vector<OtherLocalVector_>& other)
      {
        this->_gate = gate;
        this->_vector.convert(other._vector);
      }

      void sync_0()
      {
        if(_gate != nullptr)
          _gate->sync_0(_vector);
      }

      void sync_1()
      {
        if(_gate != nullptr)
          _gate->sync_1(_vector);
      }

      Vector clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return Vector(_gate, _vector.clone(mode));
      }

      void clear()
      {
        _vector.clear();
      }

      void format(DataType alpha = DataType(0))
      {
        _vector.format(alpha);
      }

      void copy(const Vector& x)
      {
        _vector.copy(*x);
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

      DataType norm2sqr() const
      {
        return dot(*this);
      }

      DataType norm2() const
      {
        return Math::sqrt(norm2sqr());
      }
    };
  } // namespace Global
} // namespace FEAST

#endif // KERNEL_GLOBAL_VECTOR_HPP
