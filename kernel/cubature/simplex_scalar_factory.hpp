#pragma once
#ifndef KERNEL_CUBATURE_SIMPLEX_SCALAR_FACTORY_HPP
#define KERNEL_CUBATURE_SIMPLEX_SCALAR_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_factory.hpp>
#include <kernel/cubature/rule.hpp>

namespace FEAST
{
  namespace Cubature
  {
    template<typename ScalarDriver_>
    class SimplexScalarFactoryBase
    {
    public:
      typedef Scalar::DriverFactory<ScalarDriver_> ScalarFactoryType;

    public:
      template<typename Weight_, typename Coord_, typename Point_>
      static void create(
        Rule<Shape::Simplex<1>, Weight_, Coord_, Point_>& rule,
        const Scalar::Rule<Weight_, Coord_>& scalar_rule)
      {
        int num_points = scalar_rule.get_num_points();
#ifdef FEAST_CUBATURE_SCALAR_PREFIX
        rule = Rule<Shape::Simplex<1>, Weight_, Coord_, Point_>(num_points, "scalar:" + scalar_rule.get_name());
#else
        rule = Rule<Shape::Simplex<1>, Weight_, Coord_, Point_>(num_points, scalar_rule.get_name());
#endif // FEAST_CUBATURE_SCALAR_PREFIX
        for(int i(0); i < num_points; ++i)
        {
          rule.get_weight(i) = scalar_rule.get_weight(i) * Weight_(0.5);
          rule.get_coord(i, 0) = (scalar_rule.get_coord(i) + Coord_(1)) * Coord_(0.5);
        }
      }

      template<typename Weight_, typename Coord_, typename Point_>
      static bool create(Rule<Shape::Simplex<1>, Weight_, Coord_, Point_>& rule, const String& name)
      {
        typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;

#ifdef FEAST_CUBATURE_SCALAR_PREFIX
        // try to find a colon within the string
        String::size_type k = name.find_first_of(':');
        if(k == name.npos)
          return false;

        // extract substrings until the colon
        String head(name.substr(0, k));
        String tail(name.substr(k + 1));

        // check head - this is the name of the formula
        if(head.trim().compare_no_case("scalar") != 0)
          return false;

        // call scalar factory to create the scalar rule
        ScalarRuleType scalar_rule;
        if(!ScalarFactoryType::create(scalar_rule, tail.trim()))
          return false;
#else
        // call scalar factory to create the scalar rule
        ScalarRuleType scalar_rule;
        if(!ScalarFactoryType::create(scalar_rule, name))
          return false;
#endif // FEAST_CUBATURE_SCALAR_PREFIX

        // convert scalar rule
        create(rule, scalar_rule);
        return true;
      }

      static String name()
      {
#ifdef FEAST_CUBATURE_SCALAR_PREFIX
        return "scalar:" + ScalarFactoryType::name();
#else
        return ScalarFactoryType::name();
#endif // FEAST_CUBATURE_SCALAR_PREFIX
      }

      template<typename Functor_>
      static void alias(Functor_& functor)
      {
#ifdef FEAST_CUBATURE_SCALAR_PREFIX
        AliasScalarPrefixFunctor<Functor_> prefix_functor(functor);
        ScalarFactoryType::alias(prefix_functor);
#else
        ScalarFactoryType::alias(functor);
#endif // FEAST_CUBATURE_SCALAR_PREFIX
      }

      /// \cond internal
    private:
#ifdef FEAST_CUBATURE_SCALAR_PREFIX
      template<typename Functor_>
      class AliasScalarPrefixFunctor
      {
      private:
        Functor_& _functor;

      public:
        explicit AliasScalarPrefixFunctor(Functor_& functor) :
          _functor(functor)
        {
        }

        void alias(const String& name)
        {
          _functor.alias("scalar:" + name);
        }

        void alias(const String& name, int num_points)
        {
          _functor.alias("scalar:" + name, num_points);
        }
      };
#endif // FEAST_CUBATURE_SCALAR_PREFIX
      /// \endcond
    };

    template<
      typename ScalarDriver_,
      bool variadic_ = ScalarDriver_::variadic>
    class SimplexScalarFactory DOXY({});

    template<typename ScalarDriver_>
    class SimplexScalarFactory<ScalarDriver_, false> :
      public SimplexScalarFactoryBase<ScalarDriver_>
    {
    public:
      typedef Shape::Simplex<1> ShapeType;
      typedef SimplexScalarFactoryBase<ScalarDriver_> BaseClass;
      typedef Scalar::DriverFactory<ScalarDriver_> ScalarFactoryType;
      static constexpr bool variadic = false;
      static constexpr int num_points = ScalarFactoryType::num_points;

    public:
      SimplexScalarFactory()
      {
      }

      using BaseClass::create;

      template<typename Weight_, typename Coord_, typename Point_>
      static void create(Rule<Shape::Simplex<1>, Weight_, Coord_, Point_>& rule)
      {
        // call scalar factory to create the scalar rule
        Scalar::Rule<Weight_, Coord_> scalar_rule;
        ScalarFactoryType::create(scalar_rule);

        // convert scalar rule
        create(rule, scalar_rule);
      }
    };

    template<typename ScalarDriver_>
    class SimplexScalarFactory<ScalarDriver_, true> :
      public SimplexScalarFactoryBase<ScalarDriver_>
    {
    public:
      typedef Shape::Simplex<1> ShapeType;
      typedef SimplexScalarFactoryBase<ScalarDriver_> BaseClass;
      typedef Scalar::DriverFactory<ScalarDriver_> ScalarFactoryType;
      static constexpr bool variadic = true;
      static constexpr int min_points = ScalarFactoryType::min_points;
      static constexpr int max_points = ScalarFactoryType::max_points;

    protected:
      int _num_points;

    public:
      explicit SimplexScalarFactory(int num_points) :
        _num_points(num_points)
      {
      }

      using BaseClass::create;

      template<typename Weight_, typename Coord_, typename Point_>
      void create(Rule<Shape::Simplex<1>, Weight_, Coord_, Point_>& rule)
      {
        create(rule, _num_points);
      }

      template<typename Weight_, typename Coord_, typename Point_>
      static void create(Rule<Shape::Simplex<1>, Weight_, Coord_, Point_>& rule, int num_points)
      {
        // call scalar factory to create the scalar rule
        Scalar::Rule<Weight_, Coord_> scalar_rule;
        ScalarFactoryType::create(scalar_rule, num_points);

        // convert scalar rule
        create(rule, scalar_rule);
      }
    };
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SIMPLEX_SCALAR_FACTORY_HPP
