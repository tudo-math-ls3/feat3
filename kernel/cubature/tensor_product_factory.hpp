// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP
#define KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/driver_factory.hpp>
#include <kernel/cubature/rule.hpp>

namespace FEAT
{
  namespace Cubature
  {
    template<typename Shape_>
    class TensorProductDriver DOXY({});

    template<>
    class TensorProductDriver<Shape::Hypercube<1> >
    {
    public:
      static int count(int num_points)
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
        for(int i(0); i < scalar_rule.get_num_points(); ++i)
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
      static int count(int num_points)
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
        int num_points = scalar_rule.get_num_points();
        for(int i(0); i < num_points; ++i)
        {
          for(int j(0); j < num_points; ++j)
          {
            int l = i*num_points + j;
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
      static int count(int num_points)
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
        int num_points = scalar_rule.get_num_points();
        for(int i(0); i < num_points; ++i)
        {
          for(int j(0); j < num_points; ++j)
          {
            for(int k(0); k < num_points; ++k)
            {
              int l = (i*num_points + j)*num_points + k;
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
        int num_points = TensorProductDriverType::count(scalar_rule.get_num_points());
#ifdef FEAT_CUBATURE_TENSOR_PREFIX
        rule = Rule<Shape_, Weight_, Coord_, Point_>(num_points, "tensor:" + scalar_rule.get_name());
#else
        rule = Rule<Shape_, Weight_, Coord_, Point_>(num_points, scalar_rule.get_name());
#endif // FEAT_CUBATURE_TENSOR_PREFIX
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

#ifdef FEAT_CUBATURE_TENSOR_PREFIX
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
#endif // FEAT_CUBATURE_TENSOR_PREFIX

        // convert scalar rule
        create(rule, scalar_rule);
        return true;
      }

      static String name()
      {
#ifdef FEAT_CUBATURE_TENSOR_PREFIX
        return "tensor:" + ScalarFactoryType::name();
#else
        return ScalarFactoryType::name();
#endif // FEAT_CUBATURE_TENSOR_PREFIX
      }

      template<typename Functor_>
      static void alias(Functor_& functor)
      {
#ifdef FEAT_CUBATURE_TENSOR_PREFIX
        AliasTensorPrefixFunctor<Functor_> prefix_functor(functor);
        ScalarFactoryType::alias(prefix_functor);
#else
        ScalarFactoryType::alias(functor);
#endif // FEAT_CUBATURE_TENSOR_PREFIX
      }

      /// \cond internal
    private:
#ifdef FEAT_CUBATURE_TENSOR_PREFIX
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

        void alias(const String& name, int num_points)
        {
          _functor.alias("tensor:" + name, num_points);
        }
      };
#endif // FEAT_CUBATURE_TENSOR_PREFIX
      /// \endcond
    };

    template<
      typename ScalarDriver_,
      typename Shape_,
      bool variadic_ = ScalarDriver_::variadic>
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
      static constexpr bool variadic = false;
      static constexpr int num_points = ScalarFactoryType::num_points;

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
      static constexpr bool variadic = true;
      static constexpr int min_points = ScalarFactoryType::min_points;
      static constexpr int max_points = ScalarFactoryType::max_points;

    protected:
      int _num_points;

    public:
      explicit TensorProductFactory(int num_points) :
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
      static void create(Rule<Shape_, Weight_, Coord_, Point_>& rule, int num_points)
      {
        // call scalar factory to create the scalar rule
        Scalar::Rule<Weight_, Coord_> scalar_rule;
        ScalarFactoryType::create(scalar_rule, num_points);

        // convert scalar rule
        create(rule, scalar_rule);
      }
    };
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_TENSOR_PRODUCT_FACTORY_HPP
