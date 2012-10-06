#pragma once
#ifndef KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/factory_wrapper.hpp>
#include <kernel/cubature/auto_alias.hpp>
#include <kernel/cubature/avail_functor.hpp>

// includes, STL
#include <set>
#include <map>

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
        // map auto-aliases
        String mapped_name(AutoAlias<ShapeType>::map(name));
        CreateFunctor functor(rule, mapped_name);
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

      static void avail(std::set<String>& names, bool aliases = true)
      {
        // list all factories except for the refine factory
        Intern::AvailSetFunctor functor(names, aliases);
        FactoryWrapper<ShapeType, WeightType, CoordType, PointType>::factory_no_refine(functor);
      }

      static void avail(std::map<String,String>& names)
      {
        // list all factories except for the refine factory
        Intern::AvailMapFunctor functor(names);
        FactoryWrapper<ShapeType, WeightType, CoordType, PointType>::factory_no_refine(functor);
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
    class DynamicFactorySelect;

    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class DynamicFactorySelect< Rule<Shape_, Weight_, Coord_, Point_> >
    {
    public:
      typedef DynamicFactory<Shape_, Weight_, Coord_, Point_> Type;
    };
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
