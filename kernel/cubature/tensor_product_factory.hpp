#pragma once
#ifndef KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP
#define KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/rule.hpp>
#include <kernel/cubature/rule.hpp>

namespace FEAST
{
  namespace Cubature
  {
    template<
      typename Policy_,
      typename ScalarFactory_,
      bool variadic = (ScalarFactory_::variadic != 0)>
    class TensorProductFactory;

    template<typename Policy_>
    class TensorProductDriverBase
    {
    public:
      template<
        typename ScalarFactory_,
        typename Functor_>
      static void factory(Functor_& functor)
      {
        functor.template factory<TensorProductFactory<Policy_, ScalarFactory_> >();
      }
    };

    template<
      typename Policy_,
      typename Shape_ = typename Policy_::ShapeType>
    class TensorProductDriver
    {
    public:
      template<
        typename ScalarFactory_,
        typename Functor_>
      static void factory(Functor_&)
      {
        // do nothing
      }
    };

    template<typename Policy_>
    class TensorProductDriver<Policy_, Shape::Simplex<1> > :
      public TensorProductDriverBase<Policy_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      typedef Scalar::Rule<WeightType, CoordType> ScalarRuleType;

      static Index count(Index num_points)
      {
        return num_points;
      }

      static void create(const ScalarRuleType& scalar_rule, RuleType& rule)
      {
        for(Index i(0); i < scalar_rule.get_num_points(); ++i)
        {
          rule.get_weight(i) = scalar_rule.get_weight(i) * WeightType(0.5);
          rule.get_coord(i, 0) = (scalar_rule.get_coord(i) + CoordType(1)) * CoordType(0.5);
        }
      }
    };

    template<typename Policy_>
    class TensorProductDriver<Policy_, Shape::Hypercube<1> > :
      public TensorProductDriverBase<Policy_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      typedef Scalar::Rule<WeightType, CoordType> ScalarRuleType;

      static Index count(Index num_points)
      {
        return num_points;
      }

      static void create(const ScalarRuleType& scalar_rule, RuleType& rule)
      {
        for(Index i(0); i < scalar_rule.get_num_points(); ++i)
        {
          rule.get_weight(i) = scalar_rule.get_weight(i);
          rule.get_coord(i, 0) = scalar_rule.get_coord(i);
        }
      }
    };

    template<typename Policy_>
    class TensorProductDriver<Policy_, Shape::Hypercube<2> > :
      public TensorProductDriverBase<Policy_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      typedef Scalar::Rule<WeightType, CoordType> ScalarRuleType;

      static Index count(Index num_points)
      {
        return num_points * num_points;
      }

      static void create(const ScalarRuleType& scalar_rule, RuleType& rule)
      {
        Index num_points = scalar_rule.get_num_points();
        for(Index i(0); i < num_points; ++i)
        {
          for(Index j(0); j < num_points; ++j)
          {
            Index l = i*num_points + j;
            rule.get_weight(l) = scalar_rule.get_weight(i) * scalar_rule.get_weight(j);
            rule.get_coord(l, 0) = scalar_rule.get_coord(i);
            rule.get_coord(l, 1) = scalar_rule.get_coord(j);
          }
        }
      }
    };

    template<typename Policy_>
    class TensorProductDriver<Policy_, Shape::Hypercube<3> > :
      public TensorProductDriverBase<Policy_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      typedef Scalar::Rule<WeightType, CoordType> ScalarRuleType;

      static Index count(Index num_points)
      {
        return num_points * num_points * num_points;
      }

      static void create(const ScalarRuleType& scalar_rule, RuleType& rule)
      {
        Index num_points = scalar_rule.get_num_points();
        for(Index i(0); i < num_points; ++i)
        {
          for(Index j(0); j < num_points; ++j)
          {
            for(Index k(0); k < num_points; ++k)
            {
              Index l = (i*num_points + j)*num_points + k;
              rule.get_weight(l) = scalar_rule.get_weight(i) * scalar_rule.get_weight(j) * scalar_rule.get_weight(k);
              rule.get_coord(l, 0) = scalar_rule.get_coord(i);
              rule.get_coord(l, 1) = scalar_rule.get_coord(j);
              rule.get_coord(l, 2) = scalar_rule.get_coord(k);
            }
          }
        }
      }
    };

    template<
      typename Policy_,
      typename ScalarFactory_>
    class TensorProductFactoryBase
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      typedef typename RuleType::PointType PointType;
      typedef Scalar::Rule<WeightType, CoordType> ScalarRuleType;
      typedef ScalarFactory_ ScalarFactoryType;

    public:
      static bool create(const String& name, RuleType& rule)
      {
#ifdef TENSOR_PREFIX
        // try to find a colon within the string
        String::size_type k = name.find_first_of(':');
        if(k == name.npos)
          return false;

        // extract substrings until the colon
        String head(name.substr(0, k));
        String tail(name.substr(k + 1));

        // check head - this is the name of the formula
        if(head.trim().compare_no_case("tensor") != 0)
          return false;

        // call scalar factory to create the scalar rule
        ScalarRuleType scalar_rule;
        if(!ScalarFactoryType::create(tail.trim(), scalar_rule))
          return false;
#else
        // call scalar factory to create the scalar rule
        ScalarRuleType scalar_rule;
        if(!ScalarFactoryType::create(name, scalar_rule))
          return false;
#endif

        // convert scalar rule
        rule = create(scalar_rule);
        return true;
      }

      static RuleType create(const ScalarRuleType& scalar_rule)
      {
        Index num_points = TensorProductDriver<Policy_>::count(scalar_rule.get_num_points());
#ifdef TENSOR_PREFIX
        RuleType rule(num_points, "tensor:" + scalar_rule.get_name());
#else
        RuleType rule(num_points, scalar_rule.get_name());
#endif
        TensorProductDriver<Policy_>::create(scalar_rule, rule);
        return rule;
      }

      static String avail_name()
      {
#ifdef TENSOR_PREFIX
        return "tensor:" + ScalarFactoryType::avail_name();
#else
        return ScalarFactoryType::avail_name();
#endif
      }
    };

    template<
      typename Policy_,
      typename ScalarFactory_>
    class TensorProductFactory<Policy_, ScalarFactory_, false> :
      public TensorProductFactoryBase<Policy_, ScalarFactory_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef TensorProductFactoryBase<Policy_, ScalarFactory_> BaseClass;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      typedef typename RuleType::PointType PointType;
      typedef Scalar::Rule<WeightType, CoordType> ScalarRuleType;
      typedef ScalarFactory_ ScalarFactoryType;
      enum
      {
        variadic = 0
      };

    public:
      TensorProductFactory()
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
      typename Policy_,
      typename ScalarFactory_>
    class TensorProductFactory<Policy_, ScalarFactory_, true> :
      public TensorProductFactoryBase<Policy_, ScalarFactory_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef TensorProductFactoryBase<Policy_, ScalarFactory_> BaseClass;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      typedef typename RuleType::PointType PointType;
      typedef Scalar::Rule<WeightType, CoordType> ScalarRuleType;
      typedef ScalarFactory_ ScalarFactoryType;
      enum
      {
        variadic = 1
      };

    protected:
      Index _num_points;

    public:
      explicit TensorProductFactory(Index num_points) :
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

#endif // KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP
