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
      namespace Intern
      {
        template<typename Mesh_, int shape_dim>
        struct DiscreteChartHelper;
      }
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
          //static constexpr int shape_dim = Mesh_::shape_dim;
          //static constexpr int world_dim = Mesh_::world_dim;


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
            return Intern::DiscreteChartHelper<Mesh_, Mesh_::shape_dim>::parse(data);
          }

      }; // class StandardChartFactory<...>

      /// \cond internal
      namespace Intern
      {
        template<typename Mesh_, int shape_dim = Mesh_::shape_dim>
        struct DiscreteChartHelper
        {
          typedef typename Mesh_::VertexSetType VertexSetType;
          static constexpr int stride = VertexSetType::stride;

          typedef typename Shape::FaceTraits<Shape::Simplex<shape_dim>, shape_dim-1>::ShapeType SurfaceShapeType;
          typedef DiscreteChart
          <
            Mesh_,
            ConformalMesh<SurfaceShapeType, Mesh_::world_dim, stride, typename Mesh_::CoordType>
          > ChartType;

          static ChartBase<Mesh_>* parse(FEAST::MeshStreamer::MeshDataContainer& data)
          {
            return ChartType::parse(data);
          }
        };

        template<typename Mesh_>
        struct DiscreteChartHelper<Mesh_, 1>
        {
          static ChartBase<Mesh_>* parse(FEAST::MeshStreamer::MeshDataContainer& DOXY(data))
          {
            return nullptr;
          }
        };
      } // namespace Intern
      /// \endcond

    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_STANDARD_CHART_FACTORY_HPP
