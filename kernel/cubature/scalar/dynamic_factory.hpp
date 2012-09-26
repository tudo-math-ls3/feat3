#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/factory_wrapper.hpp>
#include <kernel/cubature/internal_functors.hpp>

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
          Cubature::Intern::CreateFunctor<RuleType> functor(rule, name);
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

        static void avail(
          std::set<String>& names,
          bool list_aliases = true,
          bool map_aliases = true)
        {
          Cubature::Intern::AvailFunctor functor(names, list_aliases, map_aliases);
          FactoryWrapper<WeightType, CoordType>::factory(functor);
        }

        static void print_avail(
          bool list_aliases = true,
          bool map_aliases = true,
          std::ostream& stream = std::cout)
        {
          std::set<String> names;
          avail(names, list_aliases, map_aliases);
          std::set<String>::iterator it(names.begin()), jt(names.end());
          for(; it != jt; ++it)
          {
            stream << *it << std::endl;
          }
        }
      }; // class DynamicFactory<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP
