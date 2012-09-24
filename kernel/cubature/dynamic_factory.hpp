#pragma once
#ifndef KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/factory_wrapper.hpp>

namespace FEAST
{
  namespace Cubature
  {
    template<typename Policy_>
    class DynamicFactory :
      public Rule<Policy_>::Factory
    {
    public:
      typedef Policy_ Policy;
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      typedef typename RuleType::PointType PointType;
      typedef typename RuleType::Factory BaseClass;

    private:
      String _name;

      class CreateFunctor
      {
      private:
        const String& _name;
        RuleType& _rule;
        bool _okay;

      public:
        CreateFunctor(const String& name, RuleType& rule) :
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

      public:
        explicit AvailFunctor(std::set<String>& names) :
          _names(names)
        {
        }

        template<typename Factory_>
        void factory()
        {
          _names.insert(Factory_::avail_name());
        }
      };

    public:
      explicit DynamicFactory(const String& name) :
        BaseClass(),
        _name(name)
      {
      }

      virtual RuleType produce() const
      {
        return create(_name);
      }

      static bool create(const String& name, RuleType& rule)
      {
        CreateFunctor functor(name, rule);
        FactoryWrapper<Policy_>::factory(functor);
        return functor.okay();
      }

      static RuleType create(const String& name)
      {
        RuleType rule;
        if(create(name, rule))
          return rule;

        // 'name' does not match any known cubature rule
        throw InternalError("Unrecognised cubature rule name: '" + name + "'");
      }

      static void avail(std::set<String>& names)
      {
        AvailFunctor functor(names);
        FactoryWrapper<Policy_>::factory(functor);
      }
    }; // class DynamicFactory<...>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
