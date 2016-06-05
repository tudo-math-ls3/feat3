#pragma once
#ifndef KERNEL_CUBATURE_RULE_HPP
#define KERNEL_CUBATURE_RULE_HPP

// includes, FEAT
#include <kernel/shape.hpp>

// includes, system
#include <utility> // for std::move

namespace FEAT
{
  /**
   * \brief Cubature namespace
   */
  namespace Cubature
  {
    enum CtorFactory
    {
      ctor_factory
    };

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
      static constexpr int dimension = ShapeType::dimension;

    protected:
      String _name;
      int _num_points;
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

      explicit Rule(int num_points, const String& name) :
        _name(name),
        _num_points(num_points),
        _weights(nullptr),
        _points(nullptr)
      {
        if(num_points > 0)
        {
          _weights = new WeightType[size_t(num_points)];
          _points = new PointType[size_t(num_points)];
        }
      }

      template<typename Factory_>
      Rule(CtorFactory, const Factory_& factory) :
        _name(),
        _num_points(0),
        _weights(nullptr),
        _points(nullptr)
      {
        factory.create(*this);
      }

      /// move ctor
      Rule(Rule&& other) :
        _name(other._name),
        _num_points(other._num_points),
        _weights(other._weights),
        _points(other._points)
      {
        other._name.clear();
        other._num_points = 0;
        other._weights = nullptr;
        other._points = nullptr;
      }

      /// move-assign operator
      Rule& operator=(Rule&& other)
      {
        // avoid self-move
        if(this == &other)
          return *this;

        if(_weights != nullptr)
          delete [] _weights;
        if(_points != nullptr)
          delete [] _points;

        _name = other._name;
        _num_points = other._num_points;
        _weights = other._weights;
        _points = other._points;

        other._name.clear();
        other._num_points = 0;
        other._weights = nullptr;
        other._points = nullptr;

        return *this;
      }

      virtual ~Rule()
      {
        if(_points != nullptr)
        {
          delete [] _points;
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
          for(int j(0); j < dimension; ++j)
            rule._points[i][j] = _points[i][j];
        }
        return std::move(rule);
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
        ASSERT(i < _num_points);
        return _weights[i];
      }

      const WeightType& get_weight(int i) const
      {
        ASSERT(i < _num_points);
        return _weights[i];
      }

      PointType& get_point(int i)
      {
        ASSERT(i < _num_points);
        return _points[i];
      }

      const PointType& get_point(int i) const
      {
        ASSERT(i < _num_points);
        return _points[i];
      }

      CoordType& get_coord(int i, int j)
      {
        ASSERTM((i >= 0) && (i < _num_points), "point index i out-of-range");
        ASSERTM((j >= 0) && (j < dimension), "coord index j out-of-range");
        return _points[i][j];
      }

      const CoordType& get_coord(int i, int j) const
      {
        ASSERTM((i >= 0) && (i < _num_points), "point index i out-of-range");
        ASSERTM((j >= 0) && (j < dimension), "coord index j out-of-range");
        return _points[i][j];
      }
    }; // class Rule<...>
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_RULE_HPP
