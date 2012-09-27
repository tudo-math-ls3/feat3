#pragma once
#ifndef KERNEL_CUBATURE_REFINE_FACTORY_HPP
#define KERNEL_CUBATURE_REFINE_FACTORY_HPP 1

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
        typename Rule_,
        typename Shape_ = typename Rule_::ShapeType>
      class RuleRefinery;
    } // namespace Intern
    /// \endcond

    template<typename Factory_>
    class RefineFactoryBase :
      public Rule<
        typename Factory_::ShapeType,
        typename Factory_::WeightType,
        typename Factory_::CoordType,
        typename Factory_::PointType>::Factory
    {
    public:
      typedef Factory_ FactoryType;
      typedef typename Factory_::ShapeType ShapeType;
      typedef typename Factory_::WeightType WeightType;
      typedef typename Factory_::CoordType CoordType;
      typedef typename Factory_::PointType PointType;
      typedef Rule<ShapeType, WeightType, CoordType, PointType> RuleType;
      typedef Intern::RuleRefinery<RuleType> RefineryType;

      static bool create(RuleType& rule, const String& name)
      {
        // try to find a colon within the name string
        String::size_type k = name.find_first_of(':');
        if(k == name.npos)
          return false;

        // extract substrings until the colon
        String head(name.substr(0, k));
        String tail(name.substr(k + 1));

        // let's assume we refine 1 time
        Index num_refines = 1;

        // try to find a star within the head string
        k = head.find_first_of('*');
        if(k != head.npos)
        {
          // try to parse refine count
          if(!String(head.substr(k+1)).parse(num_refines))
            return false;

          // trim head of refine count
          head = head.substr(0, k);
        }

        // check head - this is the name of the formula
        if(head.trim().compare_no_case("refine") != 0)
          return false;

        // call factory to create the input rule
        RuleType rule_in;
        if(!FactoryType::create(rule_in, tail.trim()))
          return false;

        // convert input rule
        rule = create(rule_in, num_refines);
        return true;
      }

      static RuleType create(const RuleType& rule_in, Index num_refines = 1)
      {
        if(num_refines == 0)
        {
          return rule_in;
        }
        if(num_refines == 1)
        {
          RuleType rule(rule_in.get_num_points() * RefineryType::count, "refine:" + rule_in.get_name());
          RefineryType::refine(rule, rule_in);
          return rule;
        }

        // copy input rule
        RuleType rule(rule_in);
        for(Index i(0); i < num_refines; ++i)
        {
          RuleType rule_tmp(rule.get_num_points() * RefineryType::count,
            "refine*" + stringify(i+1) + ":" + rule_in.get_name());
          RefineryType::refine(rule_tmp, rule);
          rule = rule_tmp;
        }
        return rule;
      }
    };

    template<
      typename Factory_,
      bool variadic_ = (Factory_::variadic != 0)>
    class RefineFactory;

    template<typename Factory_>
    class RefineFactory<Factory_, false> :
      public RefineFactoryBase<Factory_>
    {
    public:
      typedef Factory_ FactoryType;
      typedef typename Factory_::ShapeType ShapeType;
      typedef typename Factory_::WeightType WeightType;
      typedef typename Factory_::CoordType CoordType;
      typedef typename Factory_::PointType PointType;
      typedef Rule<ShapeType, WeightType, CoordType, PointType> RuleType;
      typedef Intern::RuleRefinery<RuleType> RefineryType;
      typedef RefineFactoryBase<FactoryType> BaseClass;
      enum
      {
        variadic = 0,
        num_points = FactoryType::num_points
      };

    protected:
      Index _num_refines;

    public:
      explicit RefineFactory(Index num_refines = 1) :
        _num_refines(num_refines)
      {
      }

      virtual RuleType produce() const
      {
        return create(_num_refines);
      }

      using BaseClass::create;

      static RuleType create(Index num_refines = 1)
      {
        return BaseClass::create(FactoryType::create(), num_refines);
      }
    };

    template<typename Factory_>
    class RefineFactory<Factory_, true> :
      public RefineFactoryBase<Factory_>
    {
    public:
      typedef Factory_ FactoryType;
      typedef typename Factory_::ShapeType ShapeType;
      typedef typename Factory_::WeightType WeightType;
      typedef typename Factory_::CoordType CoordType;
      typedef typename Factory_::PointType PointType;
      typedef Rule<ShapeType, WeightType, CoordType, PointType> RuleType;
      typedef Intern::RuleRefinery<RuleType> RefineryType;
      typedef RefineFactoryBase<FactoryType> BaseClass;
      enum
      {
        variadic = 1,
        min_points = FactoryType::min_points,
        max_points = FactoryType::max_points
      };

    protected:
      Index _num_points;
      Index _num_refines;

    public:
      explicit RefineFactory(Index num_points, Index num_refines = 1) :
        _num_points(num_points),
        _num_refines(num_refines)
      {
      }

      virtual RuleType produce() const
      {
        return create(_num_points, _num_refines);
      }

      using BaseClass::create;

      static RuleType create(Index num_points, Index num_refines = 1)
      {
        return BaseClass::create(FactoryType::create(num_points), num_refines);
      }
    };

    /// \cond internal
    namespace Intern
    {
      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Simplex<1> >
      {
      public:
        typedef typename Rule_::WeightType WeightType;
        typedef typename Rule_::CoordType CoordType;
        enum
        {
          count = 2
        };

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          Index n = rule_in.get_num_points();
          for(Index i(0); i < n; ++i)
          {
            rule.get_coord(  i, 0) = CoordType(0.5) *  rule_in.get_coord(i,0);
            rule.get_coord(n+i, 0) = CoordType(0.5) * (rule_in.get_coord(i,0) + CoordType(1));
            rule.get_weight(  i) = WeightType(0.5) * rule_in.get_weight(i);
            rule.get_weight(n+i) = WeightType(0.5) * rule_in.get_weight(i);
          }
        };
      };

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Simplex<2> >
      {
      public:
        enum
        {
          count = 4
        };

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          // todo
        }
      };

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Simplex<3> >
      {
      public:
        enum
        {
          count = 12
        };

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          // todo
        }
      };

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Hypercube<1> >
      {
      public:
        enum
        {
          count = 2
        };

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          // todo
        }
      };

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Hypercube<2> >
      {
      public:
        enum
        {
          count = 4
        };

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          // todo
        }
      };

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Hypercube<3> >
      {
      public:
        enum
        {
          count = 8
        };

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          // todo
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_REFINE_FACTORY_HPP
