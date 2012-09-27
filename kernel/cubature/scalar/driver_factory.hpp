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

      /**
       * \brief Scalar Cubature Driver-Factory class template
       *
       * This template acts as a wrapper around a scalar cubature driver to implement
       * the scalar cubature factory interface.
       *
       * \tparam Driver_
       * The scalar cubature driver class template to be used.
       *
       * \tparam Weight_
       * The data type for the cubature weights.
       *
       * \tparam Coord_
       * The data type for the cubature point coordinates.
       *
       * \tparam variadic_
       * Specifies whether the driver is variadic. This template parameter is chosen automatically
       * based on the driver's \p variadic member and shall therefore not be specified explicitly.
       *
       * \author Peter Zajac
       */
      template<
        template<typename, typename> class Driver_,
        typename Weight_ = Real,
        typename Coord_ = Real,
        bool variadic_ = (Driver_<Weight_, Coord_>::variadic != 0)>
      class DriverFactory;

      /**
       * \brief DriverFactory specialisation for non-variadic drivers.
       *
       * \author Peter Zajac
       */
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

        /**
         * \brief Produces the cubature rule.
         *
         * This is an implementation of the scalar cubature factory interface.
         *
         * \returns
         * The cubature rule.
         */
        virtual RuleType produce() const
        {
          return create();
        }

        /**
         * \brief Creates the cubature rule.
         *
         * \returns
         * The cubature rule.
         */
        static RuleType create()
        {
          RuleType rule(DriverType::num_points, DriverType::name());
          DriverType::fill(rule);
          return rule;
        }

        /**
         * \brief Creates the cubature rule from a string.
         *
         * This function creates the cubature rule if the string in \p name represents a valid
         * representation of the cubature rule's name.
         *
         * \param[in] name
         * The name of the cubature rule to create.
         *
         * \param[out] rule
         * The rule to be created.
         *
         * \returns
         * \c true, if the cubature rule was created successfully, or \c false, if \p name is not
         * a valid name of the cubature rule implemented by this factory.
         */
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

        /**
         * \brief Produces the cubature rule.
         *
         * This is an implementation of the scalar cubature factory interface.
         *
         * \returns
         * The cubature rule.
         */
        virtual RuleType produce() const
        {
          return create(_num_points);
        }

        /**
         * \brief Creates the cubature rule.
         *
         * \param[in] num_points
         * The number of points for the cubature rule.
         *
         * \returns
         * The cubature rule.
         */
        static RuleType create(Index num_points)
        {
          ASSERT_(num_points >= DriverType::min_points);
          ASSERT_(num_points <= DriverType::max_points);

          RuleType rule(num_points, (DriverType::name() + ":" + stringify(num_points)));
          DriverType::fill(rule, num_points);
          return rule;
        }

        /**
         * \brief Creates the cubature rule from a string.
         *
         * This function creates the cubature rule if the string in \p name represents a valid
         * representation of the cubature rule's name.
         *
         * \param[in] name
         * The name of the cubature rule to create.
         *
         * \param[out] rule
         * The rule to be created.
         *
         * \returns
         * \c true, if the cubature rule was created successfully, or \c false, if \p name is not
         * a valid name of the cubature rule implemented by this factory.
         */
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
          if(!tail.trim().parse(num_points))
            return false;

          // try to create the rule
          rule = create(num_points);
          return true;
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
