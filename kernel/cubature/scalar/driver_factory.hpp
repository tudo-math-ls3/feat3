#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_DRIVER_FACTORY_HPP
#define KERNEL_CUBATURE_SCALAR_DRIVER_FACTORY_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/rule.hpp>

namespace FEAT
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
          bool variadic_ = Driver_::variadic>
        class DriverFactoryAliasMapper;
      } // namespace Intern
      /// \endcond

      /**
       * \brief Scalar Cubature Driver-Factory class template
       *
       * This template acts as a wrapper around a scalar cubature driver to implement
       * the scalar cubature factory interface.
       *
       * \tparam Driver_
       * The scalar cubature driver class to be used.
       *
       * \tparam variadic_
       * Specifies whether the driver is variadic. This template parameter is chosen automatically
       * based on the driver's \p variadic member and shall therefore not be specified explicitly.
       *
       * \author Peter Zajac
       */
      template<
        typename Driver_,
        bool variadic_ = Driver_::variadic>
      class DriverFactory DOXY({});

      /**
       * \brief DriverFactory specialisation for non-variadic drivers.
       *
       * \author Peter Zajac
       */
      template<typename Driver_>
      class DriverFactory<Driver_, false>
      {
      public:
        typedef Driver_ DriverType;

        static constexpr bool variadic = false;
        static constexpr bool tensorise = DriverType::tensorise;
        static constexpr int num_points = DriverType::num_points;

      public:
        DriverFactory()
        {
        }

        /**
         * \brief Creates the cubature rule.
         *
         * \param[out] rule
         * The rule to be created.
         *
         * \returns
         * \c true
         */
        template<typename Weight_, typename Coord_>
        static bool create(Rule<Weight_, Coord_>& rule)
        {
          // create the rule
          rule = Rule<Weight_, Coord_>(DriverType::num_points, DriverType::name());
          DriverType::fill(rule);
          return true;
        }

        /**
         * \brief Creates the cubature rule from a string.
         *
         * This function creates the cubature rule if the string in \p name represents a valid
         * representation of the cubature rule's name.
         *
         * \param[out] rule
         * The rule to be created.
         *
         * \param[in] name
         * The name of the cubature rule to create.
         *
         * \returns
         * \c true, if the cubature rule was created successfully, or \c false, if \p name is not
         * a valid name of the cubature rule implemented by this factory.
         */
        template<typename Weight_, typename Coord_>
        static bool create(Rule<Weight_, Coord_>& rule, const String& name)
        {
          // map alias names
          Intern::DriverFactoryAliasMapper<DriverType> mapper(name);
          DriverType::alias(mapper);
          String mapped_name(mapper.name());

          // check mapped name
          if(mapped_name.trim().compare_no_case(DriverType::name()) != 0)
            return false;

          // create the rule
          return create(rule);
        }

        /**
         * \brief Returns the name of the cubature rule.
         */
        static String name()
        {
          return DriverType::name();
        }

        /**
         * \brief Calls the driver's alias function.
         */
        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          DriverType::alias(functor);
        }
      }; // class DriverFactory<...>

      /**
       * \brief DriverFactory specialisation for variadic drivers.
       *
       * \author Peter Zajac
       */
      template<typename Driver_>
      class DriverFactory<Driver_, true>
      {
      public:
        typedef Driver_ DriverType;

        static constexpr bool variadic = true;
        static constexpr bool tensorise = DriverType::tensorise;
        static constexpr int min_points = DriverType::min_points;
        static constexpr int max_points = DriverType::max_points;

      protected:
        int _num_points;

      public:
        explicit DriverFactory(int num_points) :
          _num_points(num_points)
        {
          XASSERT(num_points >= DriverType::min_points);
          XASSERT(num_points <= DriverType::max_points);
        }

        /**
         * \brief Creates the cubature rule.
         *
         * \param[out] rule
         * The cubature rule to be created.
         */
        template<typename Weight_, typename Coord_>
        bool create(Rule<Weight_, Coord_>& rule)
        {
          return create(rule, _num_points);
        }

        /**
         * \brief Creates the cubature rule.
         *
         * \param[out] rule
         * The cubature rule to be created.
         *
         * \param[in] num_points
         * The number of points for the cubature rule.
         */
        template<typename Weight_, typename Coord_>
        static bool create(Rule<Weight_, Coord_>& rule, int num_points)
        {
          if((num_points < DriverType::min_points) || (num_points > DriverType::max_points))
            return false;

          rule = Rule<Weight_, Coord_>(num_points, (DriverType::name() + ":" + stringify(num_points)));
          DriverType::fill(rule, num_points);
          return true;
        }

        /**
         * \brief Creates the cubature rule from a string.
         *
         * This function creates the cubature rule if the string in \p name represents a valid
         * representation of the cubature rule's name.
         *
         * \param[out] rule
         * The rule to be created.
         *
         * \param[in] name
         * The name of the cubature rule to create.
         *
         * \returns
         * \c true, if the cubature rule was created successfully, or \c false, if \p name is not
         * a valid name of the cubature rule implemented by this factory.
         */
        template<typename Weight_, typename Coord_>
        static bool create(Rule<Weight_, Coord_>& rule, const String& name_in)
        {
          // map alias names
          Intern::DriverFactoryAliasMapper<DriverType> mapper(name_in);
          DriverType::alias(mapper);
          String mapped_name(mapper.name());

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
          int num_points = 0;
          if(!tail.trim().parse(num_points))
            return false;

          // try to create the rule
          return create(rule, num_points);
        }

        /**
         * \brief Returns the name of the cubature rule.
         */
        static String name()
        {
          return DriverType::name();
        }

        /**
         * \brief Calls the driver's alias function.
         */
        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          DriverType::alias(functor);
        }
      }; // class DriverFactory<...>

      /// \cond internal
      namespace Intern
      {
        template<typename Driver_>
        class DriverFactoryAliasMapper<Driver_, false>
        {
        private:
          String _name;
          bool _mapped;

        public:
          explicit DriverFactoryAliasMapper(const String& name_in) :
            _name(name_in),
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
        class DriverFactoryAliasMapper<Driver_, true>
        {
        private:
          String _name;
          bool _mapped;

        public:
          explicit DriverFactoryAliasMapper(const String& name_in) :
            _name(name_in),
            _mapped(false)
          {
          }

          void alias(const String& alias_name, int num_points)
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
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_DRIVER_FACTORY_HPP
