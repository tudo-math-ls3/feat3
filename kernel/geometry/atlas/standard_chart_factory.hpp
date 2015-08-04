#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_STANDARD_CHART_FACTORY_HPP
#define KERNEL_GEOMETRY_ATLAS_STANDARD_CHART_FACTORY_HPP 1

#include <kernel/util/mesh_streamer.hpp> // for MeshDataContainer
#include <kernel/geometry/atlas/chart.hpp>
#include <kernel/geometry/atlas/circle.hpp>
#include <kernel/geometry/atlas/discrete_chart.hpp>
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
        // Boundary mesh typedefs
        private:
          /// Simplex<1> mesh in 2d
          typedef DiscreteChart<Mesh_, ConformalMesh<Shape::Simplex<1>, 2, 2, Real>> Simplex1_2d;
          /// Simplex<1> mesh in 3d
          typedef DiscreteChart<Mesh_, ConformalMesh<Shape::Simplex<1>, 3, 3, Real>> Simplex1_3d;
          /// Simplex<2> mesh in 3d
          typedef DiscreteChart<Mesh_, ConformalMesh<Shape::Simplex<2>, 3, 3, Real>> Simplex2_3d;
          /// Hypercube<2> mesh in 3d
          typedef DiscreteChart<Mesh_, ConformalMesh<Shape::Hypercube<2>, 3, 3, Real>> Hypercube2_3d;

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

          virtual ChartBase<Mesh_>* parse_discrete_chart(FEAST::MeshStreamer::MeshDataContainer& data) override
          {
            if(data.shape_type == data.st_edge)
            {
              if(data.coord_per_vertex == 2)
                return DiscreteChart<Mesh_, ConformalMesh<Shape::Simplex<1>, 2, 2, Real>>::parse(data);
              if(data.coord_per_vertex == 3)
                return DiscreteChart<Mesh_, ConformalMesh<Shape::Simplex<1>, 3, 3, Real>>::parse(data);
              throw InternalError("Boundary mesh of shape dim 1 needs a word dimension of 2 or 3, but got "+stringify(data.coord_per_vertex));
            }

            if(data.shape_type == data.st_tria)
            {
              if(data.coord_per_vertex == 3)
                return DiscreteChart<Mesh_, ConformalMesh<Shape::Simplex<2>, 3, 3, Real>>::parse(data);
              throw InternalError("Simplex<2> boundary mesh needs a word dimension of 3, but got "+stringify(data.coord_per_vertex));
            }

            if(data.shape_type == data.st_quad)
            {
              if(data.coord_per_vertex == 3)
                return DiscreteChart<Mesh_, ConformalMesh<Shape::Hypercube<2>, 3, 3, Real>>::parse(data);
              throw InternalError("Hypercube<2> boundary mesh needs a word dimension of 3, but got "+stringify(data.coord_per_vertex));
            }

            throw InternalError("Boundary mesh needs a shape type of Simplex<1,2> or Hypercube<1,2>, but got "+stringify(data.shape_type));
          }

      }; // class StandardChartFactory<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_STANDARD_CHART_FACTORY_HPP
