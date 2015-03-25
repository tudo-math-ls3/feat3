#pragma once
#ifndef KERNEL_CUBATURE_DRIVER_FACTORY_HPP
#define KERNEL_CUBATURE_DRIVER_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/rule.hpp>

namespace FEAST
{
  namespace Cubature
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
     * \tparam Shape_
     * The shape for the cubature driver.
     *
     * \tparam variadic_
     * Specifies whether the driver is variadic. This template parameter is chosen automatically
     * based on the driver's \p variadic member and shall therefore not be specified explicitly.
     *
     * \author Peter Zajac
     */
    template<
      template<typename> class Driver_,
      typename Shape_,
      bool variadic_ = Driver_<Shape_>::variadic>
    class DriverFactory DOXY({});

    /**
     * \brief DriverFactory specialisation for non-variadic drivers.
     *
     * \author Peter Zajac
     */
    template<
      template<typename> class Driver_,
      typename Shape_>
    class DriverFactory<Driver_, Shape_, false>
    {
    public:
      typedef Driver_<Shape_> DriverType;
      typedef Shape_ ShapeType;
      static constexpr bool variadic = false;
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
      template<typename Weight_, typename Coord_, typename Point_>
      static bool create(Rule<Shape_, Weight_, Coord_, Point_>& rule)
      {
        rule = Rule<Shape_, Weight_, Coord_, Point_>(DriverType::num_points, DriverType::name());
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
      template<typename Weight_, typename Coord_, typename Point_>
      static bool create(Rule<Shape_, Weight_, Coord_, Point_>& rule, const String& name)
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
    }; // class DriverFactory<...,false>

    /**
     * \brief DriverFactory specialisation for variadic drivers.
     *
     * \author Peter Zajac
     */
    template<
      template<typename> class Driver_,
      typename Shape_>
    class DriverFactory<Driver_, Shape_, true>
    {
    public:
      typedef Driver_<Shape_> DriverType;
      typedef Shape_ ShapeType;
      static constexpr bool variadic = true;
      static constexpr int min_points = DriverType::min_points;
      static constexpr int max_points = DriverType::max_points;

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
       * \brief Creates the cubature rule.
       *
       * \param[out] rule
       * The cubature rule to be created.
       */
      template<typename Weight_, typename Coord_, typename Point_>
      bool create(Rule<Shape_, Weight_, Coord_, Point_>& rule)
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
      template<typename Weight_, typename Coord_, typename Point_>
      static bool create(Rule<Shape_, Weight_, Coord_, Point_>& rule, Index num_points)
      {
        if((num_points < DriverType::min_points) || (num_points > DriverType::max_points))
          return false;

        rule = Rule<Shape_, Weight_, Coord_, Point_>(DriverType::count(num_points),
          (DriverType::name() + ":" + stringify(num_points)));
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
      template<typename Weight_, typename Coord_, typename Point_>
      static bool create(Rule<Shape_, Weight_, Coord_, Point_>& rule, const String& name)
      {
        // map alias names
        Intern::DriverFactoryAliasMapper<DriverType> mapper(name);
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
        Index num_points = 0;
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
    }; // class DriverFactory<...,true>

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
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DRIVER_FACTORY_HPP
