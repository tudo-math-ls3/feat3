#pragma once
#ifndef KERNEL_UTIL_CPP11_SMART_POINTER_HPP
#define KERNEL_UTIL_CPP11_SMART_POINTER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

// Note: The macro 'HAVE_CPP11_SMART_POINTER' is defined by the compiler detection headers, which are
//       included by the base-header included above.
#ifdef HAVE_CPP11_SMART_POINTER
#  include <memory>
#else

// The compiler's STL implementation does not support smart pointers, so let's use our own
// fallback implementation defined below.

namespace std
{
  /**
   * \brief smart pointer proxy pattern implementation
   *
   * SmartPointer is a counted pointer and designed for use with STL or foundation containers.
   *
   * \tparam T_
   * actual object type
   *
   * \author Markus Geveler
   */
  template <class T_>
  class shared_ptr
  {
  private:

    struct Counter
    {
      Counter(T_* p = nullptr, unsigned c = 1) :
        _ptr(p),
        _count(c)
      {
      }

      T_* _ptr;
      unsigned _count;
    }* _count;

    void _acquire(Counter* c) throw()
    {
      _count = c;
      if(c != nullptr)
        ++c->_count;
    }

    void _release()
    {
      if (_count != nullptr)
      {
        if (--_count->_count == 0)
        {
          delete _count->_ptr;
          delete _count;
        }
        _count = nullptr;
      }
    }

  public:
    typedef T_ element_type;

    explicit shared_ptr(T_* p = nullptr) :
      _count(0)
    {
      if (p)
        _count = new Counter(p);
    }

    ~shared_ptr()
    {
      _release();
    }

    shared_ptr(const shared_ptr& r) throw()
    {
      _acquire(r._count);
    }

    shared_ptr& operator=(const shared_ptr& r)
    {
      if (this != &r)
      {
        _release();
        _acquire(r._count);
      }
      return *this;
    }

#if 0
    template <class Y_>
    friend class shared_ptr<Y_>;

    template <class Y_>
    shared_ptr(const shared_ptr<Y_>& r) throw()
    {
      _acquire(r._count);
    }

    template <class Y_>
    shared_ptr& operator=(const shared_ptr<Y_>& r)
    {
      if (this != &r)
      {
        _release();
        _acquire(r._count);
      }
      return *this;
    }
#endif // NO_MEMBER_TEMPLATES

    T_& operator*() const throw()
    {
      return *_count->_tr;
    }

    T_* operator->() const throw()
    {
      return _count->_ptr;
    }

    T_* get() const throw()
    {
      return (_count != nullptr ? _count->_ptr : nullptr);
    }

    bool unique() const throw()
    {
      return (_count != nullptr ? _count->_count == 1 : true);
    }
  }; // class shared_ptr<...>
} // namespace std

#endif // HAVE_CPP11_SMART_POINTER

#endif // KERNEL_UTIL_CPP11_SMART_POINTER_HPP
