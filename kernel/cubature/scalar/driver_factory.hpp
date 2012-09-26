#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_DRIVER_FACTORY_HPP
#define KERNEL_CUBATURE_SCALAR_DRIVER_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/rule.hpp>

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
          typename Driver_,
          bool variadic_ = (Driver_::variadic != 0)>
        class DriverFactoryAliasFunctor;
      } // namespace Intern
      /// \endcond

      template<
        template<typename, typename> class Driver_,
        typename Weight_ = Real,
        typename Coord_ = Real,
        bool variadic_ = (Driver_<Weight_, Coord_>::variadic != 0)>
      class DriverFactory;

      template<
        template<typename, typename> class Driver_,
        typename Weight_,
        typename Coord_>
      class DriverFactory<Driver_, Weight_, Coord_, false> :
        public Rule<Weight_, Coord_>::Factory
      {
      public:
        typedef Weight_ WeightType;
        typedef Coord_ CoordType;
        typedef Rule<Weight_, Coord_> RuleType;
        typedef Driver_<Weight_, Coord_> DriverType;
        enum
        {
          variadic = 0,
          tensorise = DriverType::tensorise,
          num_points = DriverType::num_points
        };

      public:
        DriverFactory()
        {
        }

        virtual RuleType produce() const
        {
          return create();
        }

        static RuleType create()
        {
          RuleType rule(DriverType::num_points, DriverType::name());
          DriverType::fill(rule);
          return rule;
        }

        static String name()
        {
          return DriverType::name();
        }

        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          DriverType::alias(functor);
        }

        static bool create(RuleType& rule, const String& name)
        {
          // map alias names
          Intern::DriverFactoryAliasFunctor<DriverType> functor(name);
          DriverType::alias(functor);
          String mapped_name(functor.name());

          // check mapped name
          if(mapped_name.trim().compare_no_case(DriverType::name()) != 0)
            return false;

          // create the rule
          rule = create();
          return true;
        }
      }; // class DriverFactory<...>

      template<
        template<typename,typename> class Driver_,
        typename Weight_,
        typename Coord_>
      class DriverFactory<Driver_, Weight_, Coord_, true> :
        public Rule<Weight_, Coord_>::Factory
      {
      public:
        typedef Weight_ WeightType;
        typedef Coord_ CoordType;
        typedef Rule<Weight_, Coord_> RuleType;
        typedef Driver_<Weight_, Coord_> DriverType;
        enum
        {
          variadic = 1,
          tensorise = DriverType::tensorise,
          min_points = DriverType::min_points,
          max_points = DriverType::max_points
        };

      protected:
        Index _num_points;

      public:
        explicit DriverFactory(Index num_points) :
          _num_points(num_points)
        {
          ASSERT_(num_points >= DriverType::min_points);
          ASSERT_(num_points <= DriverType::max_points);
        }

        virtual RuleType produce() const
        {
          return create(_num_points);
        }

        static RuleType create(Index num_points)
        {
          ASSERT_(num_points >= DriverType::min_points);
          ASSERT_(num_points <= DriverType::max_points);

          RuleType rule(num_points, (DriverType::name() + ":" + stringify(num_points)));
          DriverType::fill(rule, num_points);
          return rule;
        }

        static String name()
        {
          return DriverType::name();
        }

        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          DriverType::alias(functor);
        }

        static bool create(RuleType& rule, const String& name)
        {
          // map alias names
          Intern::DriverFactoryAliasFunctor<DriverType> functor(name);
          DriverType::alias(functor);
          String mapped_name(functor.name());

          // try to find a colon within the string
          String::size_type k = mapped_name.find_first_of(':');
          if(k == mapped_name.npos)
            return false;

          // extract substrings until the colon
          String head(mapped_name.substr(0, k));
          String tail(mapped_name.substr(k + 1));

          // check head - this is the name of the formula
          if(head.trim().compare_no_case(DriverType::name()) != 0)
            return false;

          // check substring
          Index num_points = 0;
          if(!parse_tail(tail.trim(), num_points))
            return false;

          // try to create the rule
          rule = create(num_points);
          return true;
        }

        static bool parse_tail(const String& tail, Index& num_points)
        {
          // must be a numerical string
          if(tail.empty() || (tail.find_first_not_of("0123456789") != tail.npos))
            return false;

          // fetch point count
          num_points = Index(atoi(tail.c_str()));
          return true;
        }
      }; // class DriverFactory<...>

      /// \cond internal
      namespace Intern
      {
        template<typename Driver_>
        class DriverFactoryAliasFunctor<Driver_, false>
        {
        private:
          String _name;
          bool _mapped;

        public:
          explicit DriverFactoryAliasFunctor(const String& name) :
            _name(name),
            _mapped(false)
          {
          }

          void alias(const String& alias_name)
          {
            if(!_mapped)
            {
              if(_name.compare_no_case(alias_name) == 0)
              {
                _name = Driver_::name();
                _mapped = true;
              }
            }
          }

          String name()
          {
            return _name;
          }
        };

        template<typename Driver_>
        class DriverFactoryAliasFunctor<Driver_, true>
        {
        private:
          String _name;
          bool _mapped;

        public:
          explicit DriverFactoryAliasFunctor(const String& name) :
            _name(name),
            _mapped(false)
          {
          }

          void alias(const String& alias_name, Index num_points)
          {
            if(!_mapped)
            {
              if(_name.compare_no_case(alias_name) == 0)
              {
                _name = Driver_::name() + ":" + stringify(num_points);
                _mapped = true;
              }
            }
          }

          String name()
          {
            return _name;
          }
        };
      } // namespace Intern
      /// \endcond
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DRIVER_FACTORY_HPP
