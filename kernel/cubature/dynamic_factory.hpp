#pragma once
#ifndef KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/factory_wrapper.hpp>

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
        FactoryWrapper<Shape_, Weight_, Coord_, Point_>::factory(functor);
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
        FactoryWrapper<Shape_, Weight_, Coord_, Point_>::factory(functor);
      }
    }; // class DynamicFactory<...>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
