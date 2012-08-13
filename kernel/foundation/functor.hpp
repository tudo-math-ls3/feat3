#pragma once
#ifndef KERNEL_FOUNDATION_FUNCTOR_HH
#define KERNEL_FOUNDATION_FUNCTOR_HH 1

#include <vector>
#include<kernel/base_header.hpp>
#include<kernel/foundation/topology_operations.hpp>
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
     * \brief STL conformal wrapper for push_back(i) functor
     *
     * \author Markus Geveler
     */
    template<
      typename ContainerType_,
      typename IndexType_,
      typename ValueType_>
    class ContainerPushBackFunctor : public FunctorBase
    {
      public:
        ContainerPushBackFunctor(ContainerType_& target, IndexType_ position, ValueType_ value) :
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
     * \brief wrapper for Topology<...>::push_back() functor
     *
     * \author Markus Geveler
     */
    template<typename TopologyType_>
    class TopologyEmptyPushBackFunctor : public FunctorBase
    {
      public:
        TopologyEmptyPushBackFunctor(TopologyType_ target, typename TopologyType_::index_type_ position) :
          _target(target),
          _position(position)
        {
        }

        virtual void execute()
        {
          //TODO
          _target.push_back();
        }

        virtual void undo()
        {
          //TODO
          TopologyElementErasure::execute(_target, _position);
        }

      private:
        TopologyType_& _target;
        typename TopologyType_::index_type_ _position;
    };
  }
}
#endif
