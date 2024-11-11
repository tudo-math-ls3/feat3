// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/cubature/scalar/factory_wrapper.hpp>

namespace FEAT
{
  namespace Cubature
  {
    namespace Scalar
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

        template<typename Weight_, typename Coord_>
        bool create(Rule<Weight_, Coord_>& rule) const
        {
          return create(rule, _name);
        }

        template<typename Weight_, typename Coord_>
        static bool create(Rule<Weight_, Coord_>& rule, const String& name)
        {
          CreateFunctor<Rule<Weight_, Coord_> > functor(rule, name);
          FactoryWrapper::factory(functor);
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
      }; // class DynamicFactory
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT
