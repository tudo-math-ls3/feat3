#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/factory_wrapper.hpp>
#include <kernel/cubature/avail_functor.hpp>

// includes, STL
#include <set>
#include <map>

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

        static bool create(RuleType& rule, const String& name)
        {
          CreateFunctor functor(rule, name);
          FactoryWrapper<WeightType, CoordType>::factory(functor);
          return functor.okay();
        }

        static RuleType create(const String& name)
        {
          RuleType rule;
          if(create(rule, name))
            return rule;

          // 'name' does not match any known cubature rule
          throw InternalError("Unrecognised cubature rule name: '" + name + "'");
        }

        static void avail(std::set<String>& names, bool aliases = true)
        {
          Cubature::Intern::AvailSetFunctor functor(names, aliases);
          FactoryWrapper<WeightType, CoordType>::factory(functor);
        }

        static void avail(std::map<String,String>& names)
        {
          // list all factories except for the refine factory
          Cubature::Intern::AvailMapFunctor functor(names);
          FactoryWrapper<WeightType, CoordType>::factory(functor);
        }

        /// \cond internal
      private:
        class CreateFunctor
        {
        private:
          RuleType& _rule;
          const String& _name;
          bool _okay;

        public:
          CreateFunctor(RuleType& rule, const String& name) :
            _rule(rule),
            _name(name),
            _okay(false)
          {
          }

          template<typename Factory_>
          void factory()
          {
            if(!_okay)
            {
              _okay = Factory_::create(_rule, _name);
            }
          }

          bool okay() const
          {
            return _okay;
          }
        };
        /// \endcond
      }; // class DynamicFactory<...>

      template<typename Rule_>
      class DynamicFactorySelect DOXY({});

      template<
        typename Weight_,
        typename Coord_>
      class DynamicFactorySelect< Rule<Weight_, Coord_> >
      {
      public:
        typedef DynamicFactory<Weight_, Coord_> Type;
      };
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP
