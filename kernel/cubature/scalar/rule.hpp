#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_RULE_HPP
#define KERNEL_CUBATURE_SCALAR_RULE_HPP

// includes, FEAST
#include <kernel/util/assertion.hpp>

namespace FEAST
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
       * \author Peter Zajac
       */
      template<
        typename WeightType_ = Real,
        typename CoordType_ = Real>
      class Rule
      {
      public:
        typedef WeightType_ WeightType;
        typedef CoordType_ CoordType;

        /**
         * \brief Scalar cubature rule factory interface
         */
        class Factory
        {
        public:
          virtual ~Factory() {}

          virtual Rule produce() const = 0;
        }; // class Rule<...>::Factory

      protected:
        String _name;
        Index _num_points;
        WeightType* _weights;
        CoordType* _coords;

      public:
        Rule() :
          _name(),
          _num_points(0),
          _weights(nullptr),
          _coords(nullptr)
        {
        }

        explicit Rule(Index num_points, const String& name) :
          _name(name),
          _num_points(num_points),
          _weights(nullptr),
          _coords(nullptr)
        {
          if(num_points > 0)
          {
            _weights = new WeightType[num_points];
            _coords = new CoordType[num_points];
          }
        }

        Rule(const Rule& other) :
          _name(),
          _num_points(0),
          _weights(nullptr),
          _coords(nullptr)
        {
          clone(other);
        }

        explicit Rule(const Factory& factory) :
          _name(),
          _num_points(0),
          _weights(nullptr),
          _coords(nullptr)
        {
          clone(factory.produce());
        }

        virtual ~Rule()
        {
          clear();
        }

        Rule& operator=(const Rule& other)
        {
          clone(other);
          return *this;
        }

        void clear()
        {
          if(_coords != nullptr)
          {
            delete [] _coords;
            _coords = nullptr;
          }
          if(_weights != nullptr)
          {
            delete [] _weights;
            _weights = nullptr;
          }
          _num_points = 0;
          _name.clear();
        }

        void create(Index num_points)
        {
          clear();
          if(num_points > 0)
          {
            _num_points = num_points;
            _weights = new WeightType[num_points];
            _coords = new CoordType[num_points];
          }
        }

        void create(Index num_points, const String& name)
        {
          create(num_points);
          _name = name;
        }

        void clone(const Rule& other)
        {
          create(other.get_num_points(), other.get_name());
          for(Index i(0); i < _num_points; ++i)
          {
            _weights[i] = other._weights[i];
            _coords[i] = other._coords[i];
          }
        }

        const String& get_name() const
        {
          return _name;
        }

        Index get_num_points() const
        {
          return _num_points;
        }

        WeightType& get_weight(Index i)
        {
          ASSERT_(_weights != nullptr);
          ASSERT(i < _num_points, "index out-of-range");
          return _weights[i];
        }

        const WeightType& get_weight(Index i) const
        {
          ASSERT_(_weights != nullptr);
          ASSERT(i < _num_points, "index out-of-range");
          return _weights[i];
        }

        CoordType& get_coord(Index i)
        {
          ASSERT_(_coords != nullptr);
          ASSERT(i < _num_points, "index out-of-range");
          return _coords[i];
        }

        const CoordType& get_coord(Index i) const
        {
          ASSERT_(_coords != nullptr);
          ASSERT(i < _num_points, "index out-of-range");
          return _coords[i];
        }
      }; // class Rule<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_RULE_HPP
