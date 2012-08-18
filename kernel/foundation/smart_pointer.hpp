#pragma once
#ifndef KERNEL_FOUNDATION_SMARTPOINTER_HH
#define KERNEL_FOUNDATION_SMARTPOINTER_HH 1

#include<kernel/base_header.hpp>

namespace FEAST
{
  namespace Foundation
  {

#define NO_MEMBER_TEMPLATES 1

    template <class T_> class SmartPointer
    {
      public:
        typedef T_ element_type;

        explicit SmartPointer(T_* p = nullptr) :
          _count(0)
        {
          if (p)
            _count = new _Counter(p);
        }

        ~SmartPointer()
        {
          _release();
        }

        SmartPointer(const SmartPointer& r) throw()
        {
          _acquire(r._count);
        }

        SmartPointer& operator=(const SmartPointer& r)
        {
          if (this != &r)
          {
            _release();
            _acquire(r._count);
          }
          return *this;
        }

#ifndef NO_MEMBER_TEMPLATES
        template <class Y> friend class SmartPointer<Y>;
        template <class Y> SmartPointer(const SmartPointer<Y>& r) throw()
        {
          _acquire(r._count);
        }

        template <class Y> SmartPointer& operator=(const SmartPointer<Y>& r)
        {
          if (this != &r)
          {
            _release();
            _acquire(r._count);
          }
          return *this;
        }
#endif // NO_MEMBER_TEMPLATES

        T_& operator*()  const throw()
        {
          return *_count->ptr;
        }

        T_* operator->() const throw()
        {
          return _count->ptr;
        }

        T_* get() const throw()
        {
          return _count ? _count->ptr : 0;
        }

        bool unique()   const throw()
        {
          return (_count ? _count->count == 1 : true);
        }

      private:

        struct _Counter
        {
          _Counter(T_* p = 0, unsigned c = 1) :
            __ptr(p),
            __count(c)
          {
          }

          T_*          __ptr;
          unsigned    __count;
        }* _count;

        void _acquire(_Counter* c) throw()
        { // increment the count
          _count = c;
          if (c)
            ++c->__count;
        }

        void _release()
        {
          // decrement the count, delete if it is 0
          if (_count)
          {
            if (--_count->__count == 0) {
              delete _count->__ptr;
              delete _count;
            }
            _count = nullptr;
          }
        }
    };
  }

}
#endif
