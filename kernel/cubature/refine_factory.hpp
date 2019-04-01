// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_REFINE_FACTORY_HPP
#define KERNEL_CUBATURE_REFINE_FACTORY_HPP 1

// includes, FEAT
#include <kernel/cubature/rule.hpp>

namespace FEAT
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

    class RefineFactoryCore
    {
    public:
      template<typename Shape_, typename Weight_, typename Coord_, typename Point_>
      static void create(
        Rule<Shape_, Weight_, Coord_, Point_>& rule,
        const Rule<Shape_, Weight_, Coord_, Point_>& rule_in,
        Index num_refines = 1)
      {
        typedef  Rule<Shape_, Weight_, Coord_, Point_> RuleType;
        typedef Intern::RuleRefinery<RuleType> RefineryType;

        if(num_refines == 0)
        {
          rule = std::move(rule_in.clone());
        }
        else if(num_refines == 1)
        {
          rule = Rule<Shape_, Weight_, Coord_, Point_>(rule_in.get_num_points() * RefineryType::count,
            "refine:" + rule_in.get_name());
          RefineryType::refine(rule, rule_in);
        }
        else
        {
          // copy input rule
          rule = std::move(rule_in.clone());
          for(Index i(0); i < num_refines; ++i)
          {
            RuleType rule_tmp(rule.get_num_points() * RefineryType::count,
              "refine*" + stringify(i+1) + ":" + rule_in.get_name());
            RefineryType::refine(rule_tmp, rule);
            rule = std::move(rule_tmp);
          }
        }
      }
    };

    template<typename Factory_>
    class RefineFactoryBase :
      public RefineFactoryCore
    {
    public:
      typedef Factory_ FactoryType;
      typedef typename Factory_::ShapeType ShapeType;
      typedef RefineFactoryCore BaseClass;

/*
      template<typename Weight_, typename Coord_, typename Point_>
      static void create(
        Rule<ShapeType, Weight_, Coord_, Point_>& rule,
        const Rule<ShapeType, Weight_, Coord_, Point_>& rule_in,
        Index num_refines = 1)
      {
        typedef  Rule<ShapeType, Weight_, Coord_, Point_> RuleType;
        typedef Intern::RuleRefinery<RuleType> RefineryType;

        if(num_refines == 0)
        {
          rule = rule_in;
        }
        if(num_refines == 1)
        {
          RuleType rule(rule_in.get_num_points() * RefineryType::count, "refine:" + rule_in.get_name());
          RefineryType::refine(rule, rule_in);
        }

        // copy input rule
        rule = rule_in;
        for(Index i(0); i < num_refines; ++i)
        {
          RuleType rule_tmp(rule.get_num_points() * RefineryType::count,
            "refine*" + stringify(i+1) + ":" + rule_in.get_name());
          RefineryType::refine(rule_tmp, rule);
          rule = rule_tmp;
        }
      }*/

      using BaseClass::create;

      template<typename Weight_, typename Coord_, typename Point_>
      static bool create(Rule<ShapeType, Weight_, Coord_, Point_>& rule, const String& name)
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
        Rule<ShapeType, Weight_, Coord_, Point_> rule_in;
        if(!FactoryType::create(rule_in, tail.trim()))
          return false;

        // convert input rule
        create(rule, rule_in, num_refines);
        return true;
      }
    };

    template<
      typename Factory_,
      bool variadic_ = Factory_::variadic>
    class RefineFactory DOXY({});

    template<typename Factory_>
    class RefineFactory<Factory_, false> :
      public RefineFactoryBase<Factory_>
    {
    public:
      typedef Factory_ FactoryType;
      typedef typename Factory_::ShapeType ShapeType;
      typedef RefineFactoryBase<FactoryType> BaseClass;
      static constexpr bool variadic = false;
      static constexpr int num_points = FactoryType::num_points;

    protected:
      Index _num_refines;

    public:
      explicit RefineFactory(Index num_refines = 1) :
        _num_refines(num_refines)
      {
      }

      using BaseClass::create;

      template<typename Weight_, typename Coord_, typename Point_>
      void create(Rule<ShapeType, Weight_, Coord_, Point_>& rule)
      {
        create(rule, _num_refines);
      }

      template<typename Weight_, typename Coord_, typename Point_>
      static void create(Rule<ShapeType, Weight_, Coord_, Point_>& rule, Index num_refines)
      {
        Rule<ShapeType, Weight_, Coord_, Point_> rule_in;
        FactoryType::create(rule_in);
        create(rule, rule_in, num_refines);
      }
    };

    template<typename Factory_>
    class RefineFactory<Factory_, true> :
      public RefineFactoryBase<Factory_>
    {
    public:
      typedef Factory_ FactoryType;
      typedef typename Factory_::ShapeType ShapeType;
      typedef RefineFactoryBase<FactoryType> BaseClass;
      static constexpr bool variadic = true;
      static constexpr int min_points = FactoryType::min_points;
      static constexpr int max_points = FactoryType::max_points;

    protected:
      int _num_points;
      Index _num_refines;

    public:
      explicit RefineFactory(int num_points, Index num_refines = 1) :
        _num_points(num_points),
        _num_refines(num_refines)
      {
      }

      using BaseClass::create;

      template<typename Weight_, typename Coord_, typename Point_>
      void create(Rule<ShapeType, Weight_, Coord_, Point_>& rule)
      {
        create(rule, _num_points, _num_refines);
      }

      template<typename Weight_, typename Coord_, typename Point_>
      static void create(Rule<ShapeType, Weight_, Coord_, Point_>& rule, int num_points, Index num_refines)
      {
        Rule<ShapeType, Weight_, Coord_, Point_> rule_in;
        FactoryType::create(rule_in, num_points);
        create(rule, rule_in, num_refines);
      }
    };

    /// \cond internal
    namespace Intern
    {
      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Simplex<1> >
      {
      public:
        typedef typename Rule_::WeightType Weight_;
        typedef typename Rule_::CoordType Coord_;
        static constexpr int count = 2;

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          int n = rule_in.get_num_points();
          for(int i(0); i < n; ++i)
          {
            rule.get_coord(  i, 0) = Coord_(0.5) *  rule_in.get_coord(i,0);
            rule.get_coord(n+i, 0) = Coord_(0.5) * (rule_in.get_coord(i,0) + Coord_(1));
            rule.get_weight(  i) = Weight_(0.5) * rule_in.get_weight(i);
            rule.get_weight(n+i) = Weight_(0.5) * rule_in.get_weight(i);
          }
        };
      }; //RuleRefinery<Rule_, Shape::Simplex<1> >

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Simplex<2> >
      {
      public:
        typedef typename Rule_::WeightType Weight_;
        typedef typename Rule_::CoordType Coord_;
        static constexpr int count = 4;

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          int n = rule_in.get_num_points();

          // points of the child triangles
          Coord_ v[4][3][2];
          v[0][0][0] = Coord_(0);
          v[0][0][1] = Coord_(0);
          v[0][1][0] = Coord_(1)/Coord_(2);
          v[0][1][1] = Coord_(0);
          v[0][2][0] = Coord_(0);
          v[0][2][1] = Coord_(1)/Coord_(2);

          v[1][0][0] = Coord_(1)/Coord_(2);
          v[1][0][1] = Coord_(0);
          v[1][1][0] = Coord_(1);
          v[1][1][1] = Coord_(0);
          v[1][2][0] = Coord_(1)/Coord_(2);
          v[1][2][1] = Coord_(1)/Coord_(2);

          v[2][0][0] = Coord_(0);
          v[2][0][1] = Coord_(1)/Coord_(2);
          v[2][1][0] = Coord_(1)/Coord_(2);
          v[2][1][1] = Coord_(1)/Coord_(2);
          v[2][2][0] = Coord_(0);
          v[2][2][1] = Coord_(1);

          v[3][0][0] = Coord_(1)/Coord_(2);
          v[3][0][1] = Coord_(1)/Coord_(2);
          v[3][1][0] = Coord_(0);
          v[3][1][1] = Coord_(1)/Coord_(2);
          v[3][2][0] = Coord_(1)/Coord_(2);
          v[3][2][1] = Coord_(0);

          // refining
          for(int i(0); i <= 3; ++i)
          {
            for(int j(0); j < n; ++j)
            {
              Coord_ x = rule_in.get_coord(j,0);
              Coord_ y = rule_in.get_coord(j,1);
              Coord_ w = 1 - x - y;
              rule.get_coord(i*n + j, 0) = w*v[i][0][0] + x*v[i][1][0] + y*v[i][2][0];
              rule.get_coord(i*n + j, 1) = w*v[i][0][1] + x*v[i][1][1] + y*v[i][2][1];
              rule.get_weight(i*n + j) = Weight_(0.25) * rule_in.get_weight(j);
            }
          }
        }
      }; //RuleRefinery<Rule_, Shape::Simplex<2> >

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Simplex<3> >
      {
      public:
        typedef typename Rule_::WeightType Weight_;
        typedef typename Rule_::CoordType Coord_;
        static constexpr int count = 12;

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          // points of the child triangles
          Coord_ v[12][4][3];

          // v_0-child nr.0

          //v_4
          v[0][0][0] = Coord_(1)/Coord_(2);
          v[0][0][1] = Coord_(0);
          v[0][0][2] = Coord_(0);

          //v_6
          v[0][1][0] = Coord_(0);
          v[0][1][1] = Coord_(0);
          v[0][1][2] = Coord_(1)/Coord_(2);

          //v_5
          v[0][2][0] = Coord_(0);
          v[0][2][1] = Coord_(1)/Coord_(2);
          v[0][2][2] = Coord_(0);

          //v_0
          v[0][3][0] = Coord_(0);
          v[0][3][1] = Coord_(0);
          v[0][3][2] = Coord_(0);


          // v_0-child nr.1

          //v_4
          v[1][0][0] = Coord_(1)/Coord_(2);
          v[1][0][1] = Coord_(0);
          v[1][0][2] = Coord_(0);

          //v_5
          v[1][1][0] = Coord_(0);
          v[1][1][1] = Coord_(1)/Coord_(2);
          v[1][1][2] = Coord_(0);

          //v_6
          v[1][2][0] = Coord_(0);
          v[1][2][1] = Coord_(0);
          v[1][2][2] = Coord_(1)/Coord_(2);

          //v_10
          v[1][3][0] = Coord_(1)/Coord_(4);
          v[1][3][1] = Coord_(1)/Coord_(4);
          v[1][3][2] = Coord_(1)/Coord_(4);


          // v_1-child nr.0

          //v_4
          v[2][0][0] = Coord_(1)/Coord_(2);
          v[2][0][1] = Coord_(0);
          v[2][0][2] = Coord_(0);

          //v_7
          v[2][1][0] = Coord_(1)/Coord_(2);
          v[2][1][1] = Coord_(1)/Coord_(2);
          v[2][1][2] = Coord_(0);

          //v_8
          v[2][2][0] = Coord_(1)/Coord_(2);
          v[2][2][1] = Coord_(0);
          v[2][2][2] = Coord_(1)/Coord_(2);

          //v_1
          v[2][3][0] = Coord_(1);
          v[2][3][1] = Coord_(0);
          v[2][3][2] = Coord_(0);


          // v_1-child nr.1

          //v_4
          v[3][0][0] = Coord_(1)/Coord_(2);
          v[3][0][1] = Coord_(0);
          v[3][0][2] = Coord_(0);

          //v_8
          v[3][1][0] = Coord_(1)/Coord_(2);
          v[3][1][1] = Coord_(0);
          v[3][1][2] = Coord_(1)/Coord_(2);

          //v_7
          v[3][2][0] = Coord_(1)/Coord_(2);
          v[3][2][1] = Coord_(1)/Coord_(2);
          v[3][2][2] = Coord_(0);

          //v_10
          v[3][3][0] = Coord_(1)/Coord_(4);
          v[3][3][1] = Coord_(1)/Coord_(4);
          v[3][3][2] = Coord_(1)/Coord_(4);


          // v_2-child nr.0

          //v_5
          v[4][0][0] = Coord_(0);
          v[4][0][1] = Coord_(1)/Coord_(2);
          v[4][0][2] = Coord_(0);

          //v_9
          v[4][1][0] = Coord_(0);
          v[4][1][1] = Coord_(1)/Coord_(2);
          v[4][1][2] = Coord_(1)/Coord_(2);

          //v_7
          v[4][2][0] = Coord_(1)/Coord_(2);
          v[4][2][1] = Coord_(1)/Coord_(2);
          v[4][2][2] = Coord_(0);

          //v_2
          v[4][3][0] = Coord_(0);
          v[4][3][1] = Coord_(1);
          v[4][3][2] = Coord_(0);


          // v_2-child nr.1

          //v_5
          v[5][0][0] = Coord_(0);
          v[5][0][1] = Coord_(1)/Coord_(2);
          v[5][0][2] = Coord_(0);

          //v_7
          v[5][1][0] = Coord_(1)/Coord_(2);
          v[5][1][1] = Coord_(1)/Coord_(2);
          v[5][1][2] = Coord_(0);

          //v_9
          v[5][2][0] = Coord_(0);
          v[5][2][1] = Coord_(1)/Coord_(2);
          v[5][2][2] = Coord_(1)/Coord_(2);

          //v_10
          v[5][3][0] = Coord_(1)/Coord_(4);
          v[5][3][1] = Coord_(1)/Coord_(4);
          v[5][3][2] = Coord_(1)/Coord_(4);


          // v_3-child nr.0

          //v_6
          v[6][0][0] = Coord_(0);
          v[6][0][1] = Coord_(0);
          v[6][0][2] = Coord_(1)/Coord_(2);

          //v_8
          v[6][1][0] = Coord_(1)/Coord_(2);
          v[6][1][1] = Coord_(0);
          v[6][1][2] = Coord_(1)/Coord_(2);

          //v_9
          v[6][2][0] = Coord_(0);
          v[6][2][1] = Coord_(1)/Coord_(2);
          v[6][2][2] = Coord_(1)/Coord_(2);

          //v_3
          v[6][3][0] = Coord_(0);
          v[6][3][1] = Coord_(0);
          v[6][3][2] = Coord_(1);


          // v_3-child nr.1

          //v_6
          v[7][0][0] = Coord_(0);
          v[7][0][1] = Coord_(0);
          v[7][0][2] = Coord_(1)/Coord_(2);

          //v_9
          v[7][1][0] = Coord_(0);
          v[7][1][1] = Coord_(1)/Coord_(2);
          v[7][1][2] = Coord_(1)/Coord_(2);

          //v_8
          v[7][2][0] = Coord_(1)/Coord_(2);
          v[7][2][1] = Coord_(0);
          v[7][2][2] = Coord_(1)/Coord_(2);

          //v_10
          v[7][3][0] = Coord_(1)/Coord_(4);
          v[7][3][1] = Coord_(1)/Coord_(4);
          v[7][3][2] = Coord_(1)/Coord_(4);


          //t_0-child

          //v_7
          v[8][0][0] = Coord_(1)/Coord_(2);
          v[8][0][1] = Coord_(1)/Coord_(2);
          v[8][0][2] = Coord_(0);

          //v_8
          v[8][1][0] = Coord_(1)/Coord_(2);
          v[8][1][1] = Coord_(0);
          v[8][1][2] = Coord_(1)/Coord_(2);

          //v_9
          v[8][2][0] = Coord_(0);
          v[8][2][1] = Coord_(1)/Coord_(2);
          v[8][2][2] = Coord_(1)/Coord_(2);

          //v_10
          v[8][3][0] = Coord_(1)/Coord_(4);
          v[8][3][1] = Coord_(1)/Coord_(4);
          v[8][3][2] = Coord_(1)/Coord_(4);


          //t_1-child

          //v_5
          v[9][0][0] = Coord_(0);
          v[9][0][1] = Coord_(1)/Coord_(2);
          v[9][0][2] = Coord_(0);

          //v_9
          v[9][1][0] = Coord_(0);
          v[9][1][1] = Coord_(1)/Coord_(2);
          v[9][1][2] = Coord_(1)/Coord_(2);

          //v_6
          v[9][2][0] = Coord_(0);
          v[9][2][1] = Coord_(0);
          v[9][2][2] = Coord_(1)/Coord_(2);

          //v_10
          v[9][3][0] = Coord_(1)/Coord_(4);
          v[9][3][1] = Coord_(1)/Coord_(4);
          v[9][3][2] = Coord_(1)/Coord_(4);


          // t_2-child

          //v_4
          v[10][0][0] = Coord_(1)/Coord_(2);
          v[10][0][1] = Coord_(0);
          v[10][0][2] = Coord_(0);

          //v_6
          v[10][1][0] = Coord_(0);
          v[10][1][1] = Coord_(0);
          v[10][1][2] = Coord_(1)/Coord_(2);

          //v_8
          v[10][2][0] = Coord_(1)/Coord_(2);
          v[10][2][1] = Coord_(0);
          v[10][2][2] = Coord_(1)/Coord_(2);

          //v_10
          v[10][3][0] = Coord_(1)/Coord_(4);
          v[10][3][1] = Coord_(1)/Coord_(4);
          v[10][3][2] = Coord_(1)/Coord_(4);


          //t_3-child

          //v_4
          v[11][0][0] = Coord_(1)/Coord_(2);
          v[11][0][1] = Coord_(0);
          v[11][0][2] = Coord_(0);

          //v_7
          v[11][1][0] = Coord_(1)/Coord_(2);
          v[11][1][1] = Coord_(1)/Coord_(2);
          v[11][1][2] = Coord_(0);

          //v_5
          v[11][2][0] = Coord_(0);
          v[11][2][1] = Coord_(1)/Coord_(2);
          v[11][2][2] = Coord_(0);

          //v_10
          v[11][3][0] = Coord_(1)/Coord_(4);
          v[11][3][1] = Coord_(1)/Coord_(4);
          v[11][3][2] = Coord_(1)/Coord_(4);


          int n = rule_in.get_num_points();
          Weight_ vol;

          // refining
          for(int i(0); i <= 11; ++i)
          {
            for(int j(0); j < n; ++j)
            {
              Coord_ x = rule_in.get_coord(j,0);
              Coord_ y = rule_in.get_coord(j,1);
              Coord_ z = rule_in.get_coord(j,2);
              Coord_ w = 1 - x - y - z;
              rule.get_coord(i*n + j, 0) = w*v[i][0][0] + x*v[i][1][0] + y*v[i][2][0] + z*v[i][3][0];
              rule.get_coord(i*n + j, 1) = w*v[i][0][1] + x*v[i][1][1] + y*v[i][2][1] + z*v[i][3][1];
              rule.get_coord(i*n + j, 2) = w*v[i][0][2] + x*v[i][1][2] + y*v[i][2][2] + z*v[i][3][2];

              // adjust the weights
              /*
              vol = (v[i][0][0] - v[i][3][0])*(v[i][1][1] - v[i][3][1])*(v[i][2][2] - v[i][3][2])
                  + (v[i][1][0] - v[i][3][0])*(v[i][2][1] - v[i][3][1])*(v[i][0][2] - v[i][3][2])
                  + (v[i][2][0] - v[i][3][0])*(v[i][0][1] - v[i][3][1])*(v[i][1][2] - v[i][3][2])
                  - (v[i][2][0] - v[i][3][0])*(v[i][1][1] - v[i][3][1])*(v[i][0][2] - v[i][3][2])
                  - (v[i][1][0] - v[i][3][0])*(v[i][0][1] - v[i][3][1])*(v[i][2][2] - v[i][3][2])
                  - (v[i][0][0] - v[i][3][0])*(v[i][2][1] - v[i][3][1])*(v[i][1][2] - v[i][3][2]);
              */
              if(i == 0 || i == 2 || i == 4|| i == 6)
              {
                vol = Weight_(1)/Weight_(8);
              }
              else
              {
                vol = Weight_(1)/Weight_(16);
              }

              rule.get_weight(i*n + j) = rule_in.get_weight(j)*Math::abs(vol);
            }
          }
        }
      }; //RuleRefinery<Rule_, Shape::Simplex<3> >

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Hypercube<1> >
      {
      public:
        typedef typename Rule_::WeightType Weight_;
        typedef typename Rule_::CoordType Coord_;
        static constexpr int count = 2;

        static void refine(Rule_& rule, const Rule_& rule_in)
        {
          int n = rule_in.get_num_points();
          for(int i(0); i < n; ++i)
          {
            rule.get_coord(  i, 0) = Coord_(0.5) * rule_in.get_coord(i,0) - Coord_(0.5);
            rule.get_coord(n+i, 0) = Coord_(0.5) * rule_in.get_coord(i,0) + Coord_(0.5);
            rule.get_weight(  i) = Weight_(0.5) * rule_in.get_weight(i);
            rule.get_weight(n+i) = Weight_(0.5) * rule_in.get_weight(i);
          }
        };

      };//RuleRefinery<Rule_, Shape::Hypercube<1> >

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Hypercube<2> >
      {
      public:
        typedef typename Rule_::WeightType Weight_;
        typedef typename Rule_::CoordType Coord_;
        static constexpr int count = 4;

        static void refine(Rule_& rule, const Rule_& rule_in)
        {

          int n = rule_in.get_num_points();

          // refining
          for(int i(0); i <= 3; ++i)
          {
            // shift
            Coord_ dx = Coord_(0);
            Coord_ dy = Coord_(0);

            if(i == 0)
            {
              dx = -Coord_(1)/Coord_(2);
              dy = -Coord_(1)/Coord_(2);
            }
            else if(i == 1)
            {
              dx = Coord_(1)/Coord_(2);
              dy = -Coord_(1)/Coord_(2);
            }
            else if(i == 2)
            {
              dx = -Coord_(1)/Coord_(2);
              dy = Coord_(1)/Coord_(2);
            }
            else if(i == 3)
            {
              dx = Coord_(1)/Coord_(2);
              dy = Coord_(1)/Coord_(2);
            }

            for(int j(0); j < n; ++j)
            {
              rule.get_coord(i*n + j, 0) = Coord_(0.5) * rule_in.get_coord(j,0) + dx;
              rule.get_coord(i*n + j, 1) = Coord_(0.5) * rule_in.get_coord(j,1) + dy;
              rule.get_weight(i*n + j) = rule_in.get_weight(j)/Weight_(4);
            }
          }
        }
      }; //RuleRefinery<Rule_, Shape::Hypercube<2> >

      template<typename Rule_>
      class RuleRefinery<Rule_, Shape::Hypercube<3> >
      {
      public:
        typedef typename Rule_::WeightType Weight_;
        typedef typename Rule_::CoordType Coord_;
        static constexpr int count = 8;

        static void refine(Rule_& rule, const Rule_& rule_in)
        {

          int n = rule_in.get_num_points();

          // refining
          for(int i(0); i <= 7; ++i)
          {
            // shift
            Coord_ dx = Coord_(0);
            Coord_ dy = Coord_(0);
            Coord_ dz = Coord_(0);

            if(i == 0)
            {
              dx = -Coord_(1)/Coord_(2);
              dy = -Coord_(1)/Coord_(2);
              dz = -Coord_(1)/Coord_(2);
            }
            else if(i == 1)
            {
              dx =  Coord_(1)/Coord_(2);
              dy = -Coord_(1)/Coord_(2);
              dz = -Coord_(1)/Coord_(2);
            }
            else if(i == 2)
            {
              dx = -Coord_(1)/Coord_(2);
              dy =  Coord_(1)/Coord_(2);
              dz = -Coord_(1)/Coord_(2);
            }
            else if(i == 3)
            {
              dx =  Coord_(1)/Coord_(2);
              dy =  Coord_(1)/Coord_(2);
              dz = -Coord_(1)/Coord_(2);
            }
            else if(i == 4)
            {
              dx = -Coord_(1)/Coord_(2);
              dy = -Coord_(1)/Coord_(2);
              dz =  Coord_(1)/Coord_(2);
            }
            else if(i == 5)
            {
              dx =  Coord_(1)/Coord_(2);
              dy = -Coord_(1)/Coord_(2);
              dz =  Coord_(1)/Coord_(2);
            }
            else if(i == 6)
            {
              dx = -Coord_(1)/Coord_(2);
              dy =  Coord_(1)/Coord_(2);
              dz =  Coord_(1)/Coord_(2);
            }
            else if(i == 7)
            {
              dx =  Coord_(1)/Coord_(2);
              dy =  Coord_(1)/Coord_(2);
              dz =  Coord_(1)/Coord_(2);
            }

            for(int j(0); j < n; ++j)
            {
              rule.get_coord(i*n + j, 0) = Coord_(0.5) * rule_in.get_coord(j,0) + dx;
              rule.get_coord(i*n + j, 1) = Coord_(0.5) * rule_in.get_coord(j,1) + dy;
              rule.get_coord(i*n + j, 2) = Coord_(0.5) * rule_in.get_coord(j,2) + dz;
              rule.get_weight(i*n + j) = rule_in.get_weight(j)/Weight_(8);
            }
          }
        }
      }; //RuleRefinery<Rule_, Shape::Hypercube<3> >

    } // namespace Intern
    /// \endcond
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_REFINE_FACTORY_HPP
