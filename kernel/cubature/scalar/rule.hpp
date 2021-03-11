// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_RULE_HPP
#define KERNEL_CUBATURE_SCALAR_RULE_HPP

// includes, FEAT
#include <kernel/util/assertion.hpp>
#include <kernel/util/string.hpp>

namespace FEAT
{
  namespace Cubature
  {
    /**
     * \brief Scalar cubature namespace
     */
    namespace Scalar
    {
      /**
       * \brief Scalar Cubature Rule class template
       *
       * \tparam Weight_
       * The datatype for the cubature weights.
       *
       * \tparam Coord_
       * The datatype for the cubature point coordinates.
       *
       * \author Peter Zajac
       */
      template<
        typename Weight_ = Real,
        typename Coord_ = Real>
      class Rule
      {
      public:
        /// Weight typedef
        typedef Weight_ WeightType;
        /// Coord typedef
        typedef Coord_ CoordType;

      protected:
        /// The name of the cubature rule.
        String _name;
        /// The total number of points in the cubature rule.
        int _num_points;
        /// Cubature weights array.
        WeightType* _weights;
        /// Cubature point coordinates array.
        CoordType* _coords;

      public:
        /// default constructor
        Rule() :
          _name(),
          _num_points(0),
          _weights(nullptr),
          _coords(nullptr)
        {
        }

        /**
         * \brief Constructor
         *
         * This constructor allocates the weights and coords arrays, but does not initialize them.
         *
         * \param[in] num_points
         * The number of points to be allocated.
         *
         * \param[in] name
         * The name of the cubature rule.
         */
        explicit Rule(int num_points, String name) :
          _name(name),
          _num_points(num_points),
          _weights(nullptr),
          _coords(nullptr)
        {
          if(num_points > 0)
          {
            _weights = new WeightType[size_t(num_points)];
            _coords = new CoordType[size_t(num_points)];
          }
        }

        /// move ctor
        Rule(Rule&& other) :
          _name(other._name),
          _num_points(other._num_points),
          _weights(other._weights),
          _coords(other._coords)
        {
          other._name.clear();
          other._num_points = 0;
          other._weights = nullptr;
          other._coords = nullptr;
        }

        /// move-assign operator
        Rule& operator=(Rule&& other)
        {
          // avoid self-move
          if(this == &other)
            return *this;

          if(_weights != nullptr)
            delete [] _weights;
          if(_coords != nullptr)
            delete [] _coords;

          _name = other._name;
          _num_points = other._num_points;
          _weights = other._weights;
          _coords = other._coords;

          other._name.clear();
          other._num_points = 0;
          other._weights = nullptr;
          other._coords = nullptr;

          return *this;
        }

        /// virtual destructor
        virtual ~Rule()
        {
          if(_coords != nullptr)
          {
            delete [] _coords;
          }
          if(_weights != nullptr)
          {
            delete [] _weights;
          }
        }

        Rule clone() const
        {
          Rule rule(_num_points, _name);
          for(int i(0); i < _num_points; ++i)
          {
            rule._weights[i] = _weights[i];
            rule._coords[i] = _coords[i];
          }
          return rule;
        }

        const String& get_name() const
        {
          return _name;
        }

        int get_num_points() const
        {
          return _num_points;
        }

        WeightType& get_weight(int i)
        {
          ASSERT(_weights != nullptr);
          ASSERTM(i < _num_points, "weight index out-of-range");
          return _weights[i];
        }

        const WeightType& get_weight(int i) const
        {
          ASSERT(_weights != nullptr);
          ASSERTM(i < _num_points, "weight index out-of-range");
          return _weights[i];
        }

        CoordType& get_coord(int i)
        {
          ASSERT(_coords != nullptr);
          ASSERTM(i < _num_points, "point index out-of-range");
          return _coords[i];
        }

        const CoordType& get_coord(int i) const
        {
          ASSERT(_coords != nullptr);
          ASSERTM(i < _num_points, "point index out-of-range");
          return _coords[i];
        }
      }; // class Rule<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_RULE_HPP
