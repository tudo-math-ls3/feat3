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
          tensorise = DriverType::tensorise
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
          DriverType::create(rule);
          return rule;
        }

        static bool create(const String& name, RuleType& rule)
        {
          if(name.trim().compare_no_case(DriverType::name()) != 0)
            return false;
          rule = create();
          return true;
        }

        static String avail_name()
        {
          return DriverType::name();
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
          tensorise = DriverType::tensorise
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
          DriverType::create(num_points, rule);
          return rule;
        }

        static bool create(const String& name, RuleType& rule)
        {
          // try to find a colon within the string
          String::size_type k = name.find_first_of(':');
          if(k == name.npos)
            return false;

          // extract substrings until the colon
          String head(name.substr(0, k));
          String tail(name.substr(k + 1));

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

        static String avail_name()
        {
          return DriverType::name() + ":<"
            + stringify(int(DriverType::min_points)) + "-"
            + stringify(int(DriverType::max_points)) + ">";
        }
      }; // class DriverFactory<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DRIVER_FACTORY_HPP
