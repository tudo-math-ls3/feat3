// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_AVAIL_FUNCTOR_HPP
#define KERNEL_CUBATURE_AVAIL_FUNCTOR_HPP 1

// includes, FEAT
#include <kernel/util/string.hpp>

// includes, STL
#include <set>
#include <map>

namespace FEAT
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Factory_,
        typename Functor_,
        bool variadic_ = Factory_::variadic>
      class AvailFunctorHelper;

      template<
        typename Factory_,
        typename Functor_>
      class AvailFunctorHelper<Factory_, Functor_, false>
      {
      private:
        Functor_& _functor;

      public:
        explicit AvailFunctorHelper(Functor_& functor) :
          _functor(functor)
        {
          _functor.add_name(Factory_::name());
        }

        void alias(const String& name)
        {
          _functor.add_alias(name, Factory_::name());
        }
      };

      template<
        typename Factory_,
        typename Functor_>
      class AvailFunctorHelper<Factory_, Functor_, true>
      {
      private:
        Functor_& _functor;

      public:
        explicit AvailFunctorHelper(Functor_& functor) :
          _functor(functor)
        {
          _functor.add_name(Factory_::name() + ":<"
            + stringify(int(Factory_::min_points)) + "-"
            + stringify(int(Factory_::max_points)) + ">");
        }

        void alias(const String& name, int num_points)
        {
          _functor.add_alias(name, Factory_::name() + ":" + stringify(num_points));
        }
      };

      class AvailSetFunctor
      {
      private:
        std::set<String>& _names;
        bool _aliases;

      public:
        explicit AvailSetFunctor(std::set<String>& names, bool aliases) :
          _names(names),
          _aliases(aliases)
        {
        }

        template<typename Factory_>
        void factory()
        {
          AvailFunctorHelper<Factory_, AvailSetFunctor> functor(*this);
          if(_aliases)
            Factory_::alias(functor);
        }

        void add_name(const String& name)
        {
          _names.insert(name);
        }

        void add_alias(const String& alias, const String& /*name*/)
        {
          if(_aliases)
            _names.insert(alias);
        }
      };

      class AvailMapFunctor
      {
      private:
        std::map<String,String>& _names;

      public:
        explicit AvailMapFunctor(std::map<String,String>& names) :
          _names(names)
        {
        }

        template<typename Factory_>
        void factory()
        {
          AvailFunctorHelper<Factory_, AvailMapFunctor> functor(*this);
          Factory_::alias(functor);
        }

        void add_name(const String& name)
        {
          _names.insert(std::make_pair(name, String()));
        }

        void add_alias(const String& alias, const String& name)
        {
          _names.insert(std::make_pair(alias, name));
        }
      };
    } // namespace internal
    /// \endcond
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_AVAIL_FUNCTOR_HPP
