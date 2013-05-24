#pragma once
#ifndef KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/factory_wrapper.hpp>
#include <kernel/cubature/auto_alias.hpp>

namespace FEAST
{
  namespace Cubature
  {
    class DynamicFactory
    {
    private:
      String _name;

    public:
      explicit DynamicFactory(const String& name) :
        _name(name)
      {
      }

      template<typename Shape_, typename Weight_, typename Coord_, typename Point_>
      bool create(Rule<Shape_, Weight_, Coord_, Point_>& rule) const
      {
        return create(rule, _name);
      }

      template<typename Shape_, typename Weight_, typename Coord_, typename Point_>
      static bool create(Rule<Shape_, Weight_, Coord_, Point_>& rule, const String& name)
      {
        // map auto-aliases
        String mapped_name(AutoAlias<Shape_>::map(name));
        CreateFunctor<Rule<Shape_, Weight_, Coord_, Point_> > functor(rule, mapped_name);
        FactoryWrapper<Shape_>::factory(functor);
        return functor.okay();
      }

      /// \cond internal
    private:
      template<typename Rule_>
      class CreateFunctor
      {
      private:
        Rule_& _rule;
        const String& _name;
        bool _okay;

      public:
        CreateFunctor(Rule_& rule, const String& name) :
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
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
