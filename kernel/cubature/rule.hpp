#pragma once
#ifndef KERNEL_CUBATURE_RULE_HPP
#define KERNEL_CUBATURE_RULE_HPP

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
{
  /**
   * \brief Cubature namespace
   */
  namespace Cubature
  {
    /**
     * \brief Cubature Rule class template
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename Weight_ = Real,
      typename Coord_ = Real,
      typename Point_ = Coord_[Shape_::dimension]>
    class Rule
    {
    public:
      typedef Shape_ ShapeType;
      typedef Weight_ WeightType;
      typedef Coord_ CoordType;
      typedef Point_ PointType;
      enum
      {
        dimension = ShapeType::dimension
      };

      /**
       * \rief Cubature rule factory interface
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
      PointType* _points;

    public:
      Rule() :
        _name(),
        _num_points(0),
        _weights(nullptr),
        _points(nullptr)
      {
      }

      explicit Rule(Index num_points, const String& name) :
        _name(name),
        _num_points(num_points),
        _weights(nullptr),
        _points(nullptr)
      {
        if(num_points > 0)
        {
          _weights = new WeightType[num_points];
          _points = new PointType[num_points];
        }
      }

      Rule(const Rule& other) :
        _name(),
        _num_points(0),
        _weights(nullptr),
        _points(nullptr)
      {
        clone(other);
      }

      Rule(const Factory& factory) :
        _name(),
        _num_points(0),
        _weights(nullptr),
        _points(nullptr)
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
        if(_points != nullptr)
        {
          delete [] _points;
          _points = nullptr;
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
          _points = new PointType[num_points];
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
          for(int j(0); j < dimension; ++j)
          {
            _points[i][j] = other._points[i][j];
          }
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
        ASSERT(i < _num_points, "index out-of-range");
        return _weights[i];
      }

      const WeightType& get_weight(Index i) const
      {
        ASSERT(i < _num_points, "index out-of-range");
        return _weights[i];
      }

      PointType& get_point(Index i)
      {
        ASSERT(i < _num_points, "index out-of-range");
        return _points[i];
      }

      const PointType& get_point(Index i) const
      {
        ASSERT(i < _num_points, "index out-of-range");
        return _points[i];
      }

      CoordType& get_coord(Index i, int j)
      {
        ASSERT(i < _num_points, "index i out-of-range");
        ASSERT((j >= 0) && (j < dimension), "index j out-of-range");
        return _points[i][j];
      }

      const CoordType& get_coord(Index i, int j) const
      {
        ASSERT(i < _num_points, "index i out-of-range");
        ASSERT((j >= 0) && (j < dimension), "index j out-of-range");
        return _points[i][j];
      }
    }; // class Rule<...>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_RULE_HPP
