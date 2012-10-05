#pragma once
#ifndef KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP
#define KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/rule.hpp>
#include <kernel/cubature/scalar/driver_factory.hpp>
#include <kernel/cubature/rule.hpp>

namespace FEAST
{
  namespace Cubature
  {
    template<
      template<typename,typename> class ScalarDriver_,
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_,
      bool variadic_ = (ScalarDriver_<Weight_, Coord_>::variadic != 0)>
    class TensorProductFactory;

    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductDriverBase
    {
    public:
      template<
        template<typename,typename> class ScalarDriver_,
        typename Functor_>
      static void scalar_driver(Functor_& functor)
      {
        functor.template factory<TensorProductFactory<ScalarDriver_, Shape_, Weight_, Coord_, Point_> >();
      }
    };

    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductDriver
    {
    public:
      template<
        template<typename,typename> class ScalarDriver_,
        typename Functor_>
      static void scalar_driver(Functor_&)
      {
        // do nothing
      }
    };

    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductDriver<Shape::Simplex<1>, Weight_, Coord_, Point_> :
      public TensorProductDriverBase<Shape::Simplex<1>, Weight_, Coord_, Point_>
    {
    public:
      typedef Rule<Shape::Simplex<1>, Weight_, Coord_, Point_> RuleType;
      typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;

      static Index count(Index num_points)
      {
        return num_points;
      }

      static void fill(RuleType& rule, const ScalarRuleType& scalar_rule)
      {
        for(Index i(0); i < scalar_rule.get_num_points(); ++i)
        {
          rule.get_weight(i) = scalar_rule.get_weight(i) * Weight_(0.5);
          rule.get_coord(i, 0) = (scalar_rule.get_coord(i) + Coord_(1)) * Coord_(0.5);
        }
      }
    };

    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductDriver<Shape::Hypercube<1>, Weight_, Coord_, Point_> :
      public TensorProductDriverBase<Shape::Hypercube<1>, Weight_, Coord_, Point_>
    {
    public:
      typedef Rule<Shape::Hypercube<1>, Weight_, Coord_, Point_> RuleType;
      typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;

      static Index count(Index num_points)
      {
        return num_points;
      }

      static void fill(RuleType& rule, const ScalarRuleType& scalar_rule)
      {
        for(Index i(0); i < scalar_rule.get_num_points(); ++i)
        {
          rule.get_weight(i) = scalar_rule.get_weight(i);
          rule.get_coord(i, 0) = scalar_rule.get_coord(i);
        }
      }
    };

    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductDriver<Shape::Hypercube<2>, Weight_, Coord_, Point_> :
      public TensorProductDriverBase<Shape::Hypercube<2>, Weight_, Coord_, Point_>
    {
    public:
      typedef Rule<Shape::Hypercube<2>, Weight_, Coord_, Point_> RuleType;
      typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;

      static Index count(Index num_points)
      {
        return num_points * num_points;
      }

      static void fill(RuleType& rule, const ScalarRuleType& scalar_rule)
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

    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductDriver<Shape::Hypercube<3>, Weight_, Coord_, Point_> :
      public TensorProductDriverBase<Shape::Hypercube<3>, Weight_, Coord_, Point_>
    {
    public:
      typedef Rule<Shape::Hypercube<3>, Weight_, Coord_, Point_> RuleType;
      typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;

      static Index count(Index num_points)
      {
        return num_points * num_points * num_points;
      }

      static void fill(RuleType& rule, const ScalarRuleType& scalar_rule)
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
      template<typename,typename> class ScalarDriver_,
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductFactoryBase
    {
    public:
      typedef Rule<Shape_, Weight_, Coord_, Point_> RuleType;
      typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;

      typedef Scalar::DriverFactory<ScalarDriver_, Weight_, Coord_> ScalarFactoryType;
      typedef TensorProductDriver<Shape_, Weight_, Coord_, Point_> TensorProductDriverType;

    public:
      static bool create(RuleType& rule, const String& name)
      {
#ifdef FEAST_CUBATURE_TENSOR_PREFIX
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
        if(!ScalarFactoryType::create(scalar_rule, tail.trim()))
          return false;
#else
        // call scalar factory to create the scalar rule
        ScalarRuleType scalar_rule;
        if(!ScalarFactoryType::create(scalar_rule, name))
          return false;
#endif // FEAST_CUBATURE_TENSOR_PREFIX

        // convert scalar rule
        rule = create(scalar_rule);
        return true;
      }

      static RuleType create(const ScalarRuleType& scalar_rule)
      {
        Index num_points = TensorProductDriverType::count(scalar_rule.get_num_points());
#ifdef FEAST_CUBATURE_TENSOR_PREFIX
        RuleType rule(num_points, "tensor:" + scalar_rule.get_name());
#else
        RuleType rule(num_points, scalar_rule.get_name());
#endif // FEAST_CUBATURE_TENSOR_PREFIX
        TensorProductDriverType::fill(rule, scalar_rule);
        return rule;
      }

      static String name()
      {
#ifdef FEAST_CUBATURE_TENSOR_PREFIX
        return "tensor:" + ScalarFactoryType::name();
#else
        return ScalarFactoryType::name();
#endif // FEAST_CUBATURE_TENSOR_PREFIX
      }

      template<typename Functor_>
      static void alias(Functor_& functor)
      {
#ifdef FEAST_CUBATURE_TENSOR_PREFIX
        AliasTensorPrefixFunctor<Functor_> prefix_functor(functor);
        ScalarFactoryType::alias(prefix_functor);
#else
        ScalarFactoryType::alias(functor);
#endif
      }

      /// \cond internal
    private:
      template<typename Functor_>
      class AliasTensorPrefixFunctor
      {
      private:
        Functor_& _functor;

      public:
        explicit AliasTensorPrefixFunctor(Functor_& functor) :
          _functor(functor)
        {
        }

        void alias(const String& name)
        {
          _functor.alias("tensor:" + name);
        }

        void alias(const String& name, Index num_points)
        {
          _functor.alias("tensor:" + name, num_points);
        }
      };
      /// \endcond
    };

    template<
      template<typename,typename> class ScalarDriver_,
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductFactory<ScalarDriver_, Shape_, Weight_, Coord_, Point_, false> :
      public TensorProductFactoryBase<ScalarDriver_, Shape_, Weight_, Coord_, Point_>
    {
    public:
      typedef Rule<Shape_, Weight_, Coord_, Point_> RuleType;
      typedef TensorProductFactoryBase<ScalarDriver_, Shape_, Weight_, Coord_, Point_> BaseClass;
      typedef TensorProductDriver<Shape_, Weight_, Coord_, Point_> TensorProductDriverType;
      typedef Shape_ ShapeType;
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
      template<typename,typename> class ScalarDriver_,
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TensorProductFactory<ScalarDriver_, Shape_, Weight_, Coord_, Point_, true> :
      public TensorProductFactoryBase<ScalarDriver_, Shape_, Weight_, Coord_, Point_>
    {
    public:
      typedef Rule<Shape_, Weight_, Coord_, Point_> RuleType;
      typedef TensorProductFactoryBase<ScalarDriver_, Shape_, Weight_, Coord_, Point_> BaseClass;
      typedef Shape_ ShapeType;
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
