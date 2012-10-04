#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_base.hpp>

// includes, STL
#include <cmath>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Gauss-Lobatto Rule driver class template
       *
       * This driver implements the Gauss-Lobatto rule.
       * \see http://en.wikipedia.org/wiki/Gaussian_quadrature
       *
       * \tparam Weight_
       * The data type for the cubature weights.
       *
       * \tparam Coord_
       * The data type for the cubature point coordinates.
       *
       * \author Constantin Christof
       */
      template<
        typename Weight_,
        typename Coord_>
      class GaussLobattoDriver :
        public DriverBase
      {
      public:
        typedef Weight_ WeightType;
        typedef Coord_ CoordType;
        typedef Rule<WeightType, CoordType> RuleType;

        enum
        {
          /// this rule is variadic
          variadic = 1,
          /// this rule has at least 3 points
          min_points = 3,
          /// this rule has at most 6 points
          max_points = 6
        };

        ///Returns the name of the cubature rule.
        static String name()
        {
          return "gauss-lobatto";
        }

        /**
         * \brief Fills the cubature rule structure.
         *
         * \param[in,out] rule
         * The cubature rule to be filled.
         * \param[in] num_points
         * The number of quadrature points.
         */
        static void fill(RuleType& rule, Index num_points)
        {
          // how many points do we have?
          switch(num_points)
          {
          case 3:
            create_3(rule);
            break;
          case 4:
            create_4(rule);
            break;
          case 5:
            create_5(rule);
            break;
          case 6:
            create_6(rule);
            break;
          }
        }

        /// Create 3-point Gauss-Lobatto rule
        static void create_3(RuleType& rule)
        {
          rule.get_coord(0) = - CoordType(1);
          rule.get_coord(1) = CoordType(0);
          rule.get_coord(2) = CoordType(1);

          rule.get_weight(0) = WeightType(1) / WeightType(3);
          rule.get_weight(1) = WeightType(4) / WeightType(3);
          rule.get_weight(2) = WeightType(1) / WeightType(3);
        }

        /// Create 4-point Gauss-Lobatto rule
        static void create_4(RuleType& rule)
        {
          rule.get_coord(0) = - CoordType(1);
          rule.get_coord(1) = -std::sqrt(CoordType(5)) / CoordType(5);
          rule.get_coord(2) = +std::sqrt(CoordType(5)) / CoordType(5);
          rule.get_coord(3) =  CoordType(1);

          rule.get_weight(0) = WeightType(1) / WeightType(6);
          rule.get_weight(1) = WeightType(5) / WeightType(6);
          rule.get_weight(2) = WeightType(5) / WeightType(6);
          rule.get_weight(3) = WeightType(1) / WeightType(6);
        }

        /// Create 5-point Gauss-Lobatto rule
        static void create_5(RuleType& rule)
        {
          rule.get_coord(0) = - CoordType(1);
          rule.get_coord(1) = -std::sqrt(CoordType(21)) / CoordType(7);
          rule.get_coord(2) = CoordType(0);
          rule.get_coord(3) = +std::sqrt(CoordType(21)) / CoordType(7);
          rule.get_coord(4) = CoordType(1);

          rule.get_weight(0) = WeightType(1) / WeightType(10);
          rule.get_weight(1) = WeightType(49) / WeightType(90);
          rule.get_weight(2) = WeightType(32) / WeightType(45);
          rule.get_weight(3) = WeightType(49) / WeightType(90);
          rule.get_weight(4) = WeightType(1) / WeightType(10);
        }

        /// Create 6-point Gauss-Lobatto rule
        static void create_6(RuleType& rule)
        {
          rule.get_coord(0) = - CoordType(1);
          rule.get_coord(1) = -std::sqrt((CoordType(7) + std::sqrt(CoordType(28)))/ CoordType(21));
          rule.get_coord(2) = -std::sqrt((CoordType(7) - std::sqrt(CoordType(28)))/ CoordType(21));
          rule.get_coord(3) = std::sqrt((CoordType(7) - std::sqrt(CoordType(28)))/ CoordType(21));
          rule.get_coord(4) = std::sqrt((CoordType(7) + std::sqrt(CoordType(28)))/ CoordType(21));
          rule.get_coord(5) = CoordType(1);


          rule.get_weight(0) = WeightType(1) / WeightType(15);
          rule.get_weight(1) = (WeightType(14) - std::sqrt(WeightType(7))) / WeightType(30);
          rule.get_weight(2) = (WeightType(14) + std::sqrt(WeightType(7))) / WeightType(30);
          rule.get_weight(3) = (WeightType(14) + std::sqrt(WeightType(7))) / WeightType(30);
          rule.get_weight(4) = (WeightType(14) - std::sqrt(WeightType(7))) / WeightType(30);
          rule.get_weight(5) = WeightType(1) / WeightType(15);
        }

      }; // class GaussLobattoDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP