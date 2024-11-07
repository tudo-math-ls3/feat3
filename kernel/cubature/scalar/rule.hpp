// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/util/assertion.hpp>
#include <kernel/util/string.hpp>

// includes, system
#include <vector>

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
        std::vector<WeightType> _weights;
        /// Cubature point coordinates array.
        std::vector<CoordType> _coords;

      public:
        /// default constructor
        Rule() :
          _name(),
          _num_points(0)
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
          _num_points(num_points)
        {
          if(num_points > 0)
          {
            _weights.resize(std::size_t(num_points));
            _coords.resize(std::size_t(num_points));
          }
        }

        /// move ctor
        Rule(Rule&& other) :
          _name(other._name),
          _num_points(other._num_points),
          _weights(std::forward<std::vector<WeightType>>(other._weights)),
          _coords(std::forward<std::vector<CoordType>>(other._coords))
        {
          other._name.clear();
          other._num_points = 0;
        }

        /// move-assign operator
        Rule& operator=(Rule&& other)
        {
          // avoid self-move
          if(this == &other)
            return *this;

          _name = other._name;
          _num_points = other._num_points;
          _weights = std::forward<std::vector<WeightType>>(other._weights);
          _coords = std::forward<std::vector<CoordType>>(other._coords);

          other._name.clear();
          other._num_points = 0;

          return *this;
        }

        /// virtual destructor
        virtual ~Rule()
        {
        }

        Rule clone() const
        {
          Rule rule(_num_points, _name);
          rule._weights = this->_weights;
          rule._coords = this->_coords;
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
          ASSERTM(i < _num_points, "weight index out-of-range");
          return _weights[std::size_t(i)];
        }

        const WeightType& get_weight(int i) const
        {
          ASSERTM(i < _num_points, "weight index out-of-range");
          return _weights[std::size_t(i)];
        }

        CoordType& get_coord(int i)
        {
          ASSERTM(i < _num_points, "point index out-of-range");
          return _coords[std::size_t(i)];
        }

        const CoordType& get_coord(int i) const
        {
          ASSERTM(i < _num_points, "point index out-of-range");
          return _coords[std::size_t(i)];
        }
      }; // class Rule<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT
