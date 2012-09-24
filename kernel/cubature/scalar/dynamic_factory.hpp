#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/factory_wrapper.hpp>

// includes, STL
#include <set>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      template<
        typename Weight_ = Real,
        typename Coord_ = Real>
      class DynamicFactory :
        public Rule<Weight_, Coord_>::Factory
      {
      public:
        typedef Weight_ WeightType;
        typedef Coord_ CoordType;
        typedef Rule<Weight_, Coord_> RuleType;
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
          FactoryWrapper<WeightType, CoordType>::factory(functor);
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
          FactoryWrapper<WeightType, CoordType>::factory(functor);
        }
      }; // class DynamicFactory<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP
