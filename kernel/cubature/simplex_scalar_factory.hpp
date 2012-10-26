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
    template<
      template<typename,typename> class ScalarDriver_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class SimplexScalarFactoryBase
    {
    public:
      typedef Rule<Shape::Simplex<1>, Weight_, Coord_, Point_> RuleType;
      typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;

      typedef Scalar::DriverFactory<ScalarDriver_, Weight_, Coord_> ScalarFactoryType;

    public:
      static bool create(RuleType& rule, const String& name)
      {
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
        rule = create(scalar_rule);
        return true;
      }

      static RuleType create(const ScalarRuleType& scalar_rule)
      {
        Index num_points = scalar_rule.get_num_points();
#ifdef FEAST_CUBATURE_SCALAR_PREFIX
        RuleType rule(num_points, "scalar:" + scalar_rule.get_name());
#else
        RuleType rule(num_points, scalar_rule.get_name());
#endif // FEAST_CUBATURE_SCALAR_PREFIX
        for(Index i(0); i < num_points; ++i)
        {
          rule.get_weight(i) = scalar_rule.get_weight(i) * Weight_(0.5);
          rule.get_coord(i, 0) = (scalar_rule.get_coord(i) + Coord_(1)) * Coord_(0.5);
        }
        return rule;
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

        void alias(const String& name, Index num_points)
        {
          _functor.alias("scalar:" + name, num_points);
        }
      };
#endif // FEAST_CUBATURE_SCALAR_PREFIX
      /// \endcond
    };

    template<
      template<typename,typename> class ScalarDriver_,
      typename Weight_,
      typename Coord_,
      typename Point_,
      bool variadic_ = (ScalarDriver_<Weight_, Coord_>::variadic != 0)>
    class SimplexScalarFactory DOXY({});

    template<
      template<typename,typename> class ScalarDriver_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class SimplexScalarFactory<ScalarDriver_, Weight_, Coord_, Point_, false> :
      public SimplexScalarFactoryBase<ScalarDriver_, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Simplex<1> ShapeType;
      typedef Rule<ShapeType, Weight_, Coord_, Point_> RuleType;
      typedef SimplexScalarFactoryBase<ScalarDriver_, Weight_, Coord_, Point_> BaseClass;
      typedef Weight_ WeightType;
      typedef Coord_ CoordType;
      typedef Point_ PointType;
      typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;
      typedef Scalar::DriverFactory<ScalarDriver_, Weight_, Coord_> ScalarFactoryType;
      enum
      {
        variadic = 0,
        num_points = ScalarFactoryType::num_points
      };

    public:
      SimplexScalarFactory()
      {
      }

      virtual RuleType produce() const
      {
        return create();
      }

      using BaseClass::create;

      static RuleType create()
      {
        return create(ScalarFactoryType::create());
      }
    };

    template<
      template<typename,typename> class ScalarDriver_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class SimplexScalarFactory<ScalarDriver_, Weight_, Coord_, Point_, true> :
      public SimplexScalarFactoryBase<ScalarDriver_, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Simplex<1> ShapeType;
      typedef Rule<ShapeType, Weight_, Coord_, Point_> RuleType;
      typedef SimplexScalarFactoryBase<ScalarDriver_, Weight_, Coord_, Point_> BaseClass;
      typedef Weight_ WeightType;
      typedef Coord_ CoordType;
      typedef Point_ PointType;
      typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;
      typedef Scalar::DriverFactory<ScalarDriver_, Weight_, Coord_> ScalarFactoryType;
      enum
      {
        variadic = 1,
        min_points = ScalarFactoryType::min_points,
        max_points = ScalarFactoryType::max_points
      };

    protected:
      Index _num_points;

    public:
      explicit SimplexScalarFactory(Index num_points) :
        _num_points(num_points)
      {
      }

      virtual RuleType produce() const
      {
        return create(_num_points);
      }

      using BaseClass::create;

      static RuleType create(Index num_points)
      {
        return create(ScalarFactoryType::create(num_points));
      }

    };
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SIMPLEX_SCALAR_FACTORY_HPP
