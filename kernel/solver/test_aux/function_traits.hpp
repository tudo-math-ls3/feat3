// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_TEST_AUX_FUNCTION_TRAITS
#define KERNEL_SOLVER_TEST_AUX_FUNCTION_TRAITS 1
#include <kernel/base_header.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <deque>
namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Class holding additional information about certain AnalyticFunction used for optimisation tests
     *
     * \tparam DT_
     * Floating point precision
     *
     * \tparam Function_
     * The class of the AnalyticFunction
     *
     * \note This generic class template is only for documentation purposes. All relevant data is to be contained in
     * specialisations in Function_.
     */
    template<typename DT_, typename Function_>
#ifndef DOXYGEN
    struct OptimisationTestTraits;
#else
    struct OptimisationTestTraits
    {
      /// The AnalyticFunction i.e. from FEAT::Analytic::Common
      typedef Function_ FunctionType;
      /// The evaluation traits for this function
      typedef Analytic::EvalTraits<DT_, FunctionType> FuncEvalTraits;
      /// Type of the evaluator
      typedef typename FunctionType::template Evaluator<FuncEvalTraits> EvalType;
      /// Type the function maps to
      typedef typename EvalType::ValueType ValueType;
      /// Type the function maps from
      typedef typename EvalType::PointType PointType;

      /**
       * \brief Gets an initial guess to start from
       *
       * \param[out] start
       * Domain point that gets set to the initial guess
       *
       */
      static void get_starting_point(PointType& start)
      {
      }

      /**
       * \brief Gets a list of all minimal points
       *
       * \param[out] min_points
       * All the domain points representing local minima of the function
       *
       */
        static void get_minimal_points(std::deque<PointType>& min_points)
        {
        }

        /**
         * \brief Gets the domain bounds
         *
         * \tparam sn_
         * stride for the Tiny::Vector holding the domain bounds
         *
         * \param[out] domain
         * Vector with 2 entries containing the lower and the upper bound of the domain
         *
         * \param[in] dim
         * Coordinate direction for which we want to have the domain bounds
         *
         * For debugging purposes, it is often handy to plot the target function on a domain. For out purposes, this
         * is always a cartesian product of intervals and this routine returns those intervals.
         *
         */
        template<int sn_>
        static void get_domain_bounds(Tiny::Vector<DT_, 2, sn_>& domain, int dim)
        {
        }

        /**
         * \brief Returns the name of the AnalyticFunction
         *
         * \returns A String containing the name of the AnalyticFunction
         */
        static String name()
        {
        }
    };
#endif

    /// \cond internal
    template<typename DT_>
    struct OptimisationTestTraits<DT_, Analytic::Common::BazaraaShettyFunction>
    {
      public:
        typedef Analytic::Common::BazaraaShettyFunction FunctionType;
        typedef Analytic::EvalTraits<DT_, FunctionType> FuncEvalTraits;
        typedef typename FunctionType::template Evaluator<FuncEvalTraits> EvalType;
        typedef typename EvalType::ValueType ValueType;
        typedef typename EvalType::PointType PointType;

        static void get_starting_point(PointType& start)
        {
          start(0) = DT_(4);
          start(1) = DT_(2);
        }

        static void get_minimal_points(std::deque<PointType>& min_points)
        {
          min_points.clear();
          PointType tmp(DT_(0));
          tmp(0) = DT_(2);
          tmp(1) = DT_(1);
          min_points.push_back(tmp);
        }

        template<int sn_>
        static void get_domain_bounds(Tiny::Vector<DT_, 2, sn_>& domain, int dim)
        {
          if( dim < 0 || dim > 2)
            throw InternalError("get_domain_bounds defined up to dim = 2, but got "+stringify(dim));

          domain(0) = -DT_(0.5);
          domain(1) = DT_(4.5);
        }

        static String name()
        {
          return "BazaraaShettyFunction";
        }
    };

    template<typename DT_>
    struct OptimisationTestTraits<DT_, Analytic::Common::GoldsteinPriceFunction>
    {
      public:
        typedef Analytic::Common::GoldsteinPriceFunction FunctionType;
        typedef Analytic::EvalTraits<DT_, FunctionType> FuncEvalTraits;
        typedef typename FunctionType::template Evaluator<FuncEvalTraits> EvalType;
        typedef typename EvalType::ValueType ValueType;
        typedef typename EvalType::PointType PointType;

        static void get_starting_point(PointType& start)
        {
          // This is a nice starting point
          //start(0) = -DT_(0.5);
          //start(1) = DT_(0.25);

          // Is is a hard starting point
          start(0) = DT_(1.5);
          start(1) = DT_(1.17);
        }

        static void get_minimal_points(std::deque<PointType>& min_points)
        {
          min_points.clear();
          PointType tmp(DT_(0));

          tmp(0) = DT_(0);
          tmp(1) = -DT_(1);
          min_points.push_back(tmp);

          tmp(0) = -DT_(0.6);
          tmp(1) = -DT_(0.4);
          min_points.push_back(tmp);

          tmp(0) = DT_(1.2);
          tmp(1) = DT_(0.8);
          min_points.push_back(tmp);

          tmp(0) = DT_(1.8);
          tmp(1) = DT_(0.2);
          min_points.push_back(tmp);
        }

        template<int sn_>
        static void get_domain_bounds(Tiny::Vector<DT_, 2, sn_>& domain, int dim)
        {
          if( dim < 0 || dim > 2)
            throw InternalError("get_domain_bounds defined up to dim = 2, but got "+stringify(dim));

          domain(0) = -DT_(0.2);
          domain(1) = DT_(0.2);
        }

        static String name()
        {
          return "GoldsteinPriceFunction";
        }
    };

    template<typename DT_>
    struct OptimisationTestTraits<DT_, Analytic::Common::HimmelblauFunction>
    {
      public:
        typedef Analytic::Common::HimmelblauFunction FunctionType;
        typedef Analytic::EvalTraits<DT_, FunctionType> FuncEvalTraits;
        typedef typename FunctionType::template Evaluator<FuncEvalTraits> EvalType;
        typedef typename EvalType::ValueType ValueType;
        typedef typename EvalType::PointType PointType;

        static void get_starting_point(PointType& start)
        {
          start(0) = DT_(1.5);
          start(1) = DT_(4);
        }

        static void get_minimal_points(std::deque<PointType>& min_points)
        {
          min_points.clear();

          PointType tmp(DT_(0));

          tmp(0) = -DT_(3.77931025337774689189076584129);
          tmp(1) = -DT_(3.28318599128616941226600051437);
          min_points.push_back(tmp);

          tmp(0) = -DT_(2.80511808695274485305357239809);
          tmp(1) =  DT_(3.13131251825057296580430072341);
          min_points.push_back(tmp);

          tmp(0) =  DT_(3);
          tmp(1) =  DT_(2);
          min_points.push_back(tmp);

          tmp(0) =  DT_(3.58442834033049174494433823938);
          tmp(1) = -DT_(1.84812652696440355353830020904);
          min_points.push_back(tmp);
        }

        template<int sn_>
        static void get_domain_bounds(Tiny::Vector<DT_, 2, sn_>& domain, int dim)
        {
          if( dim < 0 || dim > 2)
            throw InternalError("get_domain_bounds defined up to dim = 2, but got "+stringify(dim));

          domain(0) = DT_(1);
          domain(1) = DT_(4.5);
        }

        static String name()
        {
          return "HimmelblauFunction";
        }
    };

    template<typename DT_>
    struct OptimisationTestTraits<DT_, Analytic::Common::RosenbrockFunction>
    {
      public:
        typedef Analytic::Common::RosenbrockFunction FunctionType;
        typedef Analytic::EvalTraits<DT_, FunctionType> FuncEvalTraits;
        typedef typename FunctionType::template Evaluator<FuncEvalTraits> EvalType;
        typedef typename EvalType::ValueType ValueType;
        typedef typename EvalType::PointType PointType;

        static void get_starting_point(PointType& start)
        {
          start(0) = -DT_(1.9);
          start(1) = DT_(2);
        }

        static void get_minimal_points(std::deque<PointType>& min_points)
        {
          min_points.clear();
          min_points.push_back(PointType(DT_(1)));
        }

        template<int sn_>
        static void get_domain_bounds(Tiny::Vector<DT_, 2, sn_>& domain, int dim)
        {
          if(dim == 0)
          {
            domain(0) = -DT_(3);
            domain(1) = DT_(2);
          }
          else if(dim == 1)
          {
            domain(0) = -DT_(0.5);
            domain(1) = DT_(3.5);
          }
          else
            throw InternalError("get_domain_bounds defined up to dim = 2, but got "+stringify(dim));
        }

        static String name()
        {
          return "RosenbrockFunction";
        }
    };

    /// \endcond
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_TEST_AUX_FUNCTION_TRAITS
