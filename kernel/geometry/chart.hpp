#pragma once
#ifndef KERNEL_GEOMETRY_CHART_HPP
#define KERNEL_GEOMETRY_CHART_HPP 1

#include<kernel/base_header.hpp>

namespace FEAST
{
  namespace Geometry
  {

    template<Index world_dim_>
    class Chart
    {
      public:
        static constexpr Index world_dim = world_dim_;

        explicit Chart()
        {
        }

        virtual ~Chart()
        {
        }


    }

    template<Index world_dim_>
    class AnalyticChart : public Chart<world_dim_>
    {
      public:
        explicit AnalyticChart()
        {
        }

        virtual ~AnalyticChart()
        {
        }
    }

    template<Index world_dim_, typename Shape_>
    class DiscreteChart : public Chart<world_dim_>
    {
      public :
        typedef Shape_ ShapeType;

        explicit DiscreteChart()
        {
        }

        virtual ~DiscreteChart()
        {
        }

        static constexpr int shape_dim = ShapeType::dimension;
    }

    template<Index world_dim>
    class Polyline : public DiscreteChart<world_dim_, Shape::Edge>
    {

    }
  } // namespace Geometry
} //namespace FEAST

#endif
