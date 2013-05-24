#pragma once
#ifndef KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP
#define KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_factory.hpp>
#include <kernel/cubature/rule.hpp>

namespace FEAST
{
  namespace Cubature
  {
    template<typename Shape_>
    class TensorProductDriver DOXY({});

    template<>
    class TensorProductDriver<Shape::Hypercube<1> >
    {
    public:
      static Index count(Index num_points)
      {
        return num_points;
      }

      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static void fill(
        Rule<Shape::Hypercube<1>, Weight_, Coord_, Point_>& rule,
        const Scalar::Rule<Weight_, Coord_>& scalar_rule)
      {
        for(Index i(0); i < scalar_rule.get_num_points(); ++i)
        {
          rule.get_weight(i) = scalar_rule.get_weight(i);
          rule.get_coord(i, 0) = scalar_rule.get_coord(i);
        }
      }
    };

    template<>
    class TensorProductDriver<Shape::Hypercube<2> >
    {
    public:
      static Index count(Index num_points)
      {
        return num_points * num_points;
      }

      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static void fill(
        Rule<Shape::Hypercube<2>, Weight_, Coord_, Point_>& rule,
        const Scalar::Rule<Weight_, Coord_>& scalar_rule)
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

    template<>
    class TensorProductDriver<Shape::Hypercube<3> >
    {
    public:
      static Index count(Index num_points)
      {
        return num_points * num_points * num_points;
      }

      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static void fill(
        Rule<Shape::Hypercube<3>, Weight_, Coord_, Point_>& rule,
        const Scalar::Rule<Weight_, Coord_>& scalar_rule)
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
      typename ScalarDriver_,
      typename Shape_>
    class TensorProductFactoryBase
    {
    public:
      typedef Scalar::DriverFactory<ScalarDriver_> ScalarFactoryType;
      typedef TensorProductDriver<Shape_> TensorProductDriverType;

    public:

      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static void create(
        Rule<Shape_, Weight_, Coord_, Point_>& rule,
        const Scalar::Rule<Weight_, Coord_>& scalar_rule)
      {
        Index num_points = TensorProductDriverType::count(scalar_rule.get_num_points());
#ifdef FEAST_CUBATURE_TENSOR_PREFIX
        rule.create(num_points, "tensor:" + scalar_rule.get_name());
#else
        rule.create(num_points, scalar_rule.get_name());
#endif // FEAST_CUBATURE_TENSOR_PREFIX
        TensorProductDriverType::fill(rule, scalar_rule);
      }

      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static bool create(
        Rule<Shape_, Weight_, Coord_, Point_>& rule,
        const String& name)
      {
        typedef Scalar::Rule<Weight_, Coord_> ScalarRuleType;

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
        create(rule, scalar_rule);
        return true;
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
#endif // FEAST_CUBATURE_TENSOR_PREFIX
      }

      /// \cond internal
    private:
#ifdef FEAST_CUBATURE_TENSOR_PREFIX
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
#endif // FEAST_CUBATURE_TENSOR_PREFIX
      /// \endcond
    };

    template<
      typename ScalarDriver_,
      typename Shape_,
      bool variadic_ = (ScalarDriver_::variadic != 0)>
    class TensorProductFactory DOXY({});

    template<
      typename ScalarDriver_,
      typename Shape_>
    class TensorProductFactory<ScalarDriver_, Shape_, false> :
      public TensorProductFactoryBase<ScalarDriver_, Shape_>
    {
    public:
      typedef TensorProductFactoryBase<ScalarDriver_, Shape_> BaseClass;
      typedef Shape_ ShapeType;
      typedef Scalar::DriverFactory<ScalarDriver_> ScalarFactoryType;
      enum
      {
        variadic = 0,
        num_points = ScalarFactoryType::num_points
      };

    public:
      TensorProductFactory()
      {
      }

      using BaseClass::create;

      template<typename Weight_, typename Coord_, typename Point_>
      static void create(Rule<Shape_, Weight_, Coord_, Point_>& rule)
      {
        // call scalar factory to create the scalar rule
        Scalar::Rule<Weight_, Coord_> scalar_rule;
        ScalarFactoryType::create(scalar_rule);

        // convert scalar rule
        create(rule, scalar_rule);
      }
    };

    template<
      typename ScalarDriver_,
      typename Shape_>
    class TensorProductFactory<ScalarDriver_, Shape_, true> :
      public TensorProductFactoryBase<ScalarDriver_, Shape_>
    {
    public:
      typedef TensorProductFactoryBase<ScalarDriver_, Shape_> BaseClass;
      typedef Shape_ ShapeType;
      typedef Scalar::DriverFactory<ScalarDriver_> ScalarFactoryType;
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

      using BaseClass::create;

      template<typename Weight_, typename Coord_, typename Point_>
      void create(Rule<Shape_, Weight_, Coord_, Point_>& rule)
      {
        create(rule, _num_points);
      }

      template<typename Weight_, typename Coord_, typename Point_>
      static void create(Rule<Shape_, Weight_, Coord_, Point_>& rule, Index num_points)
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

#endif // KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP
