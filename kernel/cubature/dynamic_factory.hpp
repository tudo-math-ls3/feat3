#pragma once
#ifndef KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/factory_wrapper.hpp>
#include <kernel/cubature/internal_functors.hpp>

// includes, STL
#include <set>

namespace FEAST
{
  namespace Cubature
  {
    template<
      typename Shape_,
      typename Weight_ = Real,
      typename Coord_ = Real,
      typename Point_ = Coord_[Shape_::dimension]>
    class DynamicFactory :
      public Rule<Shape_, Weight_, Coord_, Point_>::Factory
    {
    public:
      typedef Rule<Shape_, Weight_, Coord_, Point_> RuleType;
      typedef typename RuleType::Factory BaseClass;
      typedef Shape_ ShapeType;
      typedef Weight_ WeightType;
      typedef Coord_ CoordType;
      typedef Point_ PointType;

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
        Intern::CreateFunctor<RuleType> functor(rule, name);
        FactoryWrapper<ShapeType, WeightType, CoordType, PointType>::factory(functor);
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
        Intern::AvailFunctor functor(names, list_aliases, map_aliases);
        FactoryWrapper<ShapeType, WeightType, CoordType, PointType>::factory(functor);
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
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
