#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP
#define KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/factory_wrapper.hpp>

// includes, STL
#include <set>
#include <map>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      /// \cond internal
      namespace Intern
      {
        template<
          typename Factory_,
          bool variadic_ = (Factory_::variadic != 0)>
        class DynamicFactoryAvailFunctor;
        template<
          typename Factory_,
          bool variadic_ = (Factory_::variadic != 0)>
        class DynamicFactoryAliasFunctor;
      } // namespace internal
      /// \endcond

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
          AvailFunctor functor(names, aliases);
          FactoryWrapper<WeightType, CoordType>::factory(functor);
        }

        static void alias(std::map<String,String>& names)
        {
          // list all factories except for the refine factory
          AliasFunctor functor(names);
          FactoryWrapper<WeightType, CoordType>::factory(functor);
        }

        static void print_avail(
          bool list_aliases = true,
          bool map_aliases = true,
          std::ostream& stream = std::cout)
        {
          if(!list_aliases || (list_aliases && !map_aliases))
          {
            std::set<String> names;
            avail(names, list_aliases);
            std::set<String>::iterator it(names.begin()), jt(names.end());
            for(; it != jt; ++it)
            {
              stream << *it << std::endl;
            }
          }
          else
          {
            std::map<String,String> names;
            alias(names);
            {
              std::set<String> names2;
              avail(names2, false);
              std::set<String>::iterator it(names2.begin()), jt(names2.end());
              for(; it != jt; ++it)
              {
                names.insert(std::make_pair(*it, String()));
              }
            }
            std::map<String,String>::iterator it(names.begin()), jt(names.end());
            for(; it != jt; ++it)
            {
              stream << it->first;
              if(!it->second.empty())
                stream << " [" << it->second << "]";
              stream << std::endl;
            }
          }
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

        class AvailFunctor
        {
        private:
          std::set<String>& _names;
          bool _aliases;

        public:
          explicit AvailFunctor(std::set<String>& names, bool aliases) :
            _names(names),
            _aliases(aliases)
          {
          }

          template<typename Factory_>
          void factory()
          {
            Intern::DynamicFactoryAvailFunctor<Factory_> functor(_names);
            if(_aliases)
              Factory_::alias(functor);
          }
        };

        class AliasFunctor
        {
        private:
          std::map<String,String>& _names;

        public:
          explicit AliasFunctor(std::map<String,String>& names) :
            _names(names)
          {
          }

          template<typename Factory_>
          void factory()
          {
            Intern::DynamicFactoryAliasFunctor<Factory_> functor(_names);
            Factory_::alias(functor);
          }
        };
        /// \endcond
      }; // class DynamicFactory<...>

      template<typename Rule_>
      class DynamicFactorySelect;

      template<
        typename Weight_,
        typename Coord_>
      class DynamicFactorySelect< Rule<Weight_, Coord_> >
      {
      public:
        typedef DynamicFactory<Weight_, Coord_> Type;
      };


      /// \cond internal
      namespace Intern
      {
        template<typename Factory_>
        class DynamicFactoryAvailFunctor<Factory_, false>
        {
        private:
          std::set<String>& _names;

        public:
          explicit DynamicFactoryAvailFunctor(std::set<String>& names) :
            _names(names)
          {
            _names.insert(Factory_::name());
          }

          void alias(const String& name)
          {
            _names.insert(name);
          }
        };

        template<typename Factory_>
        class DynamicFactoryAvailFunctor<Factory_, true>
        {
        private:
          std::set<String>& _names;

        public:
          explicit DynamicFactoryAvailFunctor(std::set<String>& names) :
            _names(names)
          {
            _names.insert(Factory_::name() + ":<"
              + stringify(int(Factory_::min_points)) + "-"
              + stringify(int(Factory_::max_points)) + ">");
          }

          void alias(const String& name, Index /*num_points*/)
          {
            _names.insert(name);
          }
        };

        template<typename Factory_>
        class DynamicFactoryAliasFunctor<Factory_, false>
        {
        private:
          std::map<String,String>& _names;

        public:
          explicit DynamicFactoryAliasFunctor(std::map<String,String>& names) :
            _names(names)
          {
          }

          void alias(const String& name)
          {
            _names.insert(std::make_pair(name, Factory_::name()));
          }
        };

        template<typename Factory_>
        class DynamicFactoryAliasFunctor<Factory_, true>
        {
        private:
          std::map<String,String>& _names;

        public:
          explicit DynamicFactoryAliasFunctor(std::map<String,String>& names) :
            _names(names)
          {
            _names.insert(std::make_pair(String(Factory_::name() + ":<"
              + stringify(int(Factory_::min_points)) + "-"
              + stringify(int(Factory_::max_points)) + ">"), String()));
          }

          void alias(const String& name, Index num_points)
          {
            _names.insert(std::make_pair(name, String(Factory_::name() + ":" + stringify(num_points))));
          }
        };
      } // namespace internal
      /// \endcond
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DYNAMIC_FACTORY_HPP
