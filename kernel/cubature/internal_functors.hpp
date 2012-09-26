#pragma once
#ifndef KERNEL_CUBATURE_INTERNAL_FUNCTORS_HPP
#define KERNEL_CUBATURE_INTERNAL_FUNCTORS_HPP 1

// includes, FEAST
#include <kernel/util/string.hpp>

// includes, STL
#include <set>

namespace FEAST
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Factory_,
        bool variadic_ = (Factory_::variadic != 0)>
      class DynamicFactoryAvailFunctorHelper;

      template<typename Factory_>
      class DynamicFactoryAvailFunctorHelper<Factory_, false>
      {
      private:
        std::set<String>& _names;
        bool _map_aliases;

      public:
        explicit DynamicFactoryAvailFunctorHelper(std::set<String>& names, bool map_aliases) :
          _names(names),
          _map_aliases(map_aliases)
        {
          // Insert factory name
          _names.insert(Factory_::name());
        }

        void alias(const String& name)
        {
          if(_map_aliases)
            _names.insert(String(name + " [" + Factory_::name() + "]"));
          else
            _names.insert(name);
        }
      };

      template<typename Factory_>
      class DynamicFactoryAvailFunctorHelper<Factory_, true>
      {
      private:
        std::set<String>& _names;
        bool _map_aliases;

      public:
        explicit DynamicFactoryAvailFunctorHelper(std::set<String>& names, bool map_aliases) :
          _names(names),
          _map_aliases(map_aliases)
        {
          // Insert factory name
          _names.insert(Factory_::name() + ":<"
            + stringify(int(Factory_::min_points)) + "-"
            + stringify(int(Factory_::max_points)) + ">");
        }

        void alias(const String& name, Index num_points)
        {
          if(_map_aliases)
            _names.insert(String(name + " [" + Factory_::name() + ":" + stringify(num_points) + "]"));
          else
            _names.insert(name);
        }
      };

      template<typename RuleType_>
      class CreateFunctor
      {
      private:
        const String& _name;
        RuleType_& _rule;
        bool _okay;

      public:
        CreateFunctor(const String& name, RuleType_& rule) :
          _name(name),
          _rule(rule),
          _okay(false)
        {
        }

        template<typename Factory_>
        void factory()
        {
          if(!_okay)
          {
            _okay = Factory_::create(_name, _rule);
          }
        }

        bool okay() const
        {
          return _okay;
        }
      };

      class AvailFunctor
      {
      private:
        std::set<String>& _names;
        bool _list_aliases;
        bool _map_aliases;

      public:
        explicit AvailFunctor(std::set<String>& names, bool list_aliases, bool map_aliases) :
          _names(names),
          _list_aliases(list_aliases),
          _map_aliases(map_aliases)
        {
        }

        template<typename Factory_>
        void factory()
        {
          DynamicFactoryAvailFunctorHelper<Factory_> functor(_names, _map_aliases);
          if(_list_aliases)
            Factory_::alias(functor);
        }
      };

    } // namespace Intern
    /// \endcond
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
