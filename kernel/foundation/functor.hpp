#pragma once
#ifndef KERNEL_FOUNDATION_FUNCTOR_HH
#define KERNEL_FOUNDATION_FUNCTOR_HH 1

#include <vector>
#include<kernel/base_header.hpp>
#include<kernel/foundation/functor_error.hpp>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief FunctorBase wraps Foundation operation functors
     *
     * \author Markus Geveler
     */
    class FunctorBase
    {
      public:
        virtual void execute() = 0;
        virtual void undo() = 0;

      protected:
        bool _executed;
        bool _undone;
    };


    /**
     * \brief STL conformal push_back(i) functor
     *
     * \author Markus Geveler
     */
    template<
      typename ContainerType_,
      typename IndexType_,
      typename ValueType_>
    class PushBackFunctor : public FunctorBase
    {
      public:
        PushBackFunctor(ContainerType_& target, IndexType_ position, ValueType_ value) :
          _target(target),
          _position(position),
          _value(value)
        {
          this->_executed = false;
          this->_undone = false;
        }

        virtual void execute()
        {
          if(this->_executed)
            throw FunctorError("Already executed!");

          _target.push_back(_value);
          this->_executed = true;
        }

        virtual void undo()
        {
          if(this->_undone)
            throw FunctorError("Already undone!");

          _target.erase(_target.begin() + _position);
          this->_undone = true;
        }

        ContainerType_& get_target()
        {
          return _target;
        }

        IndexType_ get_position()
        {
          return _position;
        }

        ValueType_ get_value()
        {
          return _value;
        }

      private:
        ContainerType_& _target;
        IndexType_ _position;
        ValueType_ _value;
    };

    /**
     * \brief push_back() functor
     *
     * \author Markus Geveler
     */
    template<typename ContainerType_>
    class EmptyPushBackFunctor : public FunctorBase
    {
      public:
        EmptyPushBackFunctor(ContainerType_ target, typename ContainerType_::index_type_ position) :
          _target(target),
          _position(position)
        {
          this->_executed = false;
          this->_undone = false;
        }

        virtual void execute()
        {
          if(this->_executed)
            throw FunctorError("Already executed!");

          _target.push_back();
          this->_executed = true;
        }

        virtual void undo()
        {
          if(this->_undone)
            throw FunctorError("Already undone!");

          _target.erase(_target.begin() + _position);
          this->_undone = true;
        }

      private:
        ContainerType_& _target;
        typename ContainerType_::index_type_ _position;
    };
  }
}
#endif
