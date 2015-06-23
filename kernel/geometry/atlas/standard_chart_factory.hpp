#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_STANDARD_CHART_FACTORY_HPP
#define KERNEL_GEOMETRY_ATLAS_STANDARD_CHART_FACTORY_HPP 1

#include <kernel/geometry/atlas/chart.hpp>
#include <kernel/geometry/atlas/circle.hpp>
#include <kernel/geometry/atlas/polyline.hpp>
#include <kernel/geometry/atlas/sphere.hpp>
#include <kernel/geometry/atlas/tube.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Atlas
    {
      /**
       * \brief Standard Chart factory class template.
       *
       * This class template implements the ChartFactory interface and acts as a wrapper
       * around all charts implemented in the Atlas namespace.
       *
       * \tparam Mesh_
       * The type of the mesh to be parameterised.
       */
      template<typename Mesh_>
      class StandardChartFactory:
        public ChartFactory<Mesh_>
      {
      public:
        /** \copydoc ChartFactory::parse_chart() */
        virtual ChartBase<Mesh_>* parse_chart(const String& type, const std::deque<String>& data, const Index line) override
        {
          // let's check the type
               if(type == "circle")       return Circle  <Mesh_>::parse(data, line);
          else if(type == "polyline")     return Polyline<Mesh_>::parse(data, line);
          else if(type == "sphere")       return Sphere  <Mesh_>::parse(data, line);
          else if(type == "tube")         return Tube    <Mesh_>::parse(data, line);

          // unknown
          return nullptr;
        }
      }; // class StandardChartFactory<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_STANDARD_CHART_FACTORY_HPP
