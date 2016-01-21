#pragma once
#ifndef KERNEL_FOUNDATION_FUNCTOR_HH
#define KERNEL_FOUNDATION_FUNCTOR_HH 1

#include <kernel/base_header.hpp>
#include <kernel/foundation/functor_error.hpp>

#include <vector>
#include <memory>

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
        virtual const std::string type_name() = 0;

        virtual ~FunctorBase()
        {
        }

        bool is_executed()
        {
          return _executed;
        }

      protected:
        bool _executed;
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
        }

        virtual const std::string type_name() override
        {
          return "PushBackFunctor";
        }

        virtual void execute() override
        {
          if(this->_executed)
            throw FunctorError("Already executed!");

          _target.push_back(_value);
          this->_executed = true;
        }

        virtual void undo() override
        {
          if(!(this->_executed))
            throw FunctorError("Already undone!");

          _target.pop_back();
          this->_executed = false;
        }

        PushBackFunctor& operator=(const PushBackFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_target = rhs._target;
          this->_position = rhs._position;
          this->_value = rhs._value;
          this->_executed = rhs._executed;
          return *this;
        }

        PushBackFunctor(const PushBackFunctor& other) :
          _target(other._target),
          _position(other._position),
          _value(other._value)
        {
          this->_executed = other._executed;
        }

      private:
        ContainerType_& _target ;
        IndexType_ _position;
        ValueType_ _value;
    };

    /**
     * \brief STL conformal erase(i) functor
     *
     * \author Markus Geveler
     */
    template<
      typename ContainerType_,
      typename IndexType_,
      typename ValueType_>
    class EraseFunctor : public FunctorBase
    {
      public:
        EraseFunctor(ContainerType_& target, IndexType_ position, ValueType_ value) :
          _target(target),
          _position(position),
          _value(value)
        {
          this->_executed = false;
        }

        virtual const std::string type_name() override
        {
          return "EraseFunctor";
        }

        virtual void execute() override
        {
          if(this->_executed)
            throw FunctorError("Already executed!");

          _target.erase(_target.begin() + (typename ContainerType_::difference_type)_position);
          this->_executed = true;
        }

        virtual void undo() override
        {
          if(!(this->_executed))
            throw FunctorError("Already undone!");

          _target.insert(_target.begin() + (typename ContainerType_::difference_type)_position, _value);
          this->_executed = false;
        }

        EraseFunctor& operator=(const EraseFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_target = rhs._target;
          this->_position = rhs._position;
          this->_value = rhs._value;
          this->_executed = rhs._executed;
          return *this;
        }

        EraseFunctor(const EraseFunctor& other) :
          _target(other._target),
          _position(other._position),
          _value(other._value)
        {
          this->_executed = other._executed;
        }

      private:
        ContainerType_& _target;
        IndexType_ _position;
        ValueType_ _value;
    };

    /**
     * \brief compound functor
     *
     * \author Markus Geveler
     */
    template<template<typename, typename> class StorageType_ = std::vector>
    class CompoundFunctor : public FunctorBase
    {
      public:
        typedef StorageType_<std::shared_ptr<FunctorBase>, std::allocator<std::shared_ptr<FunctorBase> > > storage_type_;

        CompoundFunctor(bool executed = false) :
          _functors()
        {
          this->_executed = executed;
        }

        virtual const std::string type_name() override
        {
          return "CompoundFunctor";
        }

        void add_functor(FunctorBase* functor)
        {
          _functors.push_back(std::shared_ptr<FunctorBase>(functor));
        }

        void add_functor(const std::shared_ptr<FunctorBase>& functor)
        {
          _functors.push_back(functor);
        }

        virtual void execute() override
        {
          if(this->_executed)
            throw FunctorError("Already executed!");

          for(Index i(0) ; i < _functors.size() ; ++i)
          {
            _functors.at(i)->execute();
          }
          this->_executed = true;
        }

        virtual void undo() override
        {
          if(!(this->_executed))
            throw FunctorError("Already undone!");

          if(_functors.size() == 0)
            throw FunctorError("Already cleared!");

          if(_functors.size() != 0)
          {
            Index i(Index(_functors.size()) - 1);
            while(true)
            {
              _functors.at(i)->undo();

              if(i != 0)
                --i;
              else
                break;
            }
            this->_executed = false;
          }
        }

        storage_type_& get_functors()
        {
          return _functors;
        }

        Index size()
        {
          return Index(_functors.size());
        }

        const storage_type_& get_functors() const
        {
          return _functors;
        }

        Index size() const
        {
          return Index(_functors.size());
        }

        CompoundFunctor& operator=(const CompoundFunctor& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_functors = rhs._functors;
          this->_executed = rhs._executed;
          return *this;
        }

        CompoundFunctor(const CompoundFunctor& other) :
          _functors(other._functors)
        {
          this->_executed = other._executed;
        }

      protected:
        storage_type_ _functors;
    };
  }
}
#endif
