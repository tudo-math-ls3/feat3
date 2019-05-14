// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP 1

// includes, FEAT
#include <kernel/cubature/factory_wrapper.hpp>
#include <kernel/cubature/auto_alias.hpp>
#include <kernel/util/exception.hpp>

namespace FEAT
{
  namespace Cubature
  {
    class UnknownRule :
      public Exception
    {
    public:
      explicit UnknownRule(const String & rule, const String& shape) :
        Exception(String("Unknown cubature rule '") + rule + "' for shape " + shape)
      {
      }
      virtual ~UnknownRule() throw()
      {
      }
    };

    class DynamicFactory
    {
    private:
      String _name;

    public:
      DynamicFactory()
      {
      }

      explicit DynamicFactory(String name_) :
        _name(name_)
      {
      }

      DynamicFactory(const DynamicFactory& other) :
        _name(other._name)
      {
      }

      DynamicFactory& operator=(const DynamicFactory& other)
      {
        _name = other._name;
        return *this;
      }

      const String& name() const
      {
        return _name;
      }

      template<typename Shape_, typename Weight_, typename Coord_, typename Point_>
      void create_throw(Rule<Shape_, Weight_, Coord_, Point_>& rule) const
      {
        if(!create(rule, _name))
          throw UnknownRule(_name, Shape_::name());
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

#ifdef FEAT_EICKT
    extern template bool DynamicFactory::create<Shape::Simplex<2>, Real, Real, Tiny::Vector<Real, 2, 2>>
      (Rule<Shape::Simplex<2>, Real, Real, Tiny::Vector<Real, 2, 2>>&) const;
    extern template bool DynamicFactory::create<Shape::Simplex<3>, Real, Real, Tiny::Vector<Real, 3, 3>>
      (Rule<Shape::Simplex<3>, Real, Real, Tiny::Vector<Real, 3, 3>>&) const;
    extern template bool DynamicFactory::create<Shape::Hypercube<2>, Real, Real, Tiny::Vector<Real, 2, 2>>
      (Rule<Shape::Hypercube<2>, Real, Real, Tiny::Vector<Real, 2, 2>>&) const;
    extern template bool DynamicFactory::create<Shape::Hypercube<3>, Real, Real, Tiny::Vector<Real, 3, 3>>
      (Rule<Shape::Hypercube<3>, Real, Real, Tiny::Vector<Real, 3, 3>>&) const;
#endif // FEAT_EICKT
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_DYNAMIC_FACTORY_HPP
