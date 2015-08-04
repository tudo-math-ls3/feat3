#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_DISCRETE_CHART_HPP
#define KERNEL_GEOMETRY_ATLAS_DISCRETE_CHART_HPP 1

#include <kernel/geometry/atlas/chart.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>

#include <kernel/geometry/export_vtk.hpp>
namespace FEAST
{
  namespace Geometry
  {
    namespace Atlas
    {

      /// Discrete chart traits
      template<typename SurfaceMeshType>
      struct DiscreteChartTraits
      {
        /// No explicit map is available in general
        static constexpr bool is_explicit = false;
        /// We support implicit projection
        static constexpr bool is_implicit = true;
        /// This is a world_dim dimensional object
        static constexpr int world_dim = SurfaceMeshType::world_dim;
        /// If there was a parametrisation, it would be the object's shape dim
        static constexpr int param_dim = SurfaceMeshType::shape_dim;
      }; // struct DiscreteChartTraits

      /**
       * \brief Class template for discrete charts
       *
       * A DiscreteChart is just a (possibly very fine) mesh, that may even have a completely different shape type
       * than the mesh referring to it. It will in general never have an explicit form (a 1d DiscreteChart has it,
       * but then one could use a Polyline as well).
       *
       * \tparam Mesh_
       * Type of the (root) mesh this chart refers to.
       *
       * \tparam SurfaceMesh_
       * Type of the boundary mesh defining the DiscreteChart.
       *
       * Mesh_ could be ConformalMesh<Hypercube<3>> and SurfaceMesh_ could be ConformalMesh<Simplex<2>>.
       *
       */
      template<typename Mesh_, typename SurfaceMesh_>
      class DiscreteChart:
        public ChartCRTP<DiscreteChart<Mesh_, SurfaceMesh_>, Mesh_, DiscreteChartTraits<SurfaceMesh_>>
      {
        public:
          typedef SurfaceMesh_ SurfaceMeshType;
          typedef ChartCRTP<DiscreteChart<Mesh_, SurfaceMeshType>, Mesh_, DiscreteChartTraits<SurfaceMesh_>> BaseClass;
          typedef typename BaseClass::CoordType DataType;
          typedef typename BaseClass::WorldPoint WorldPoint;

        protected:
          /// Pointer to the surface mesh object
          SurfaceMeshType* _surface_mesh;

        private:
          /**
           * \brief Constructor
           * See the documentation of this class template for details.
           */
          explicit DiscreteChart(SurfaceMeshType* surface_mesh_) :
            _surface_mesh(surface_mesh_)
          {
          }

        public:
          /**
           * \brief Constructor
           * See the documentation of this class template for details.
           */
          explicit DiscreteChart() :
            _surface_mesh(nullptr)
          {
          }

          /**
           * \brief Destructor
           */
          virtual ~DiscreteChart()
          {
            if(_surface_mesh != nullptr)
              delete _surface_mesh;
          }

          /**
           * \brief Parses a MeshDataContainer
           *
           * \returns
           * A new object of type DiscreteChart containing the mesh data from the container
           */
          static DiscreteChart<Mesh_, SurfaceMeshType>* parse(FEAST::MeshStreamer::MeshDataContainer& data)
          {
            // Create a factory for our mesh
            MeshStreamerFactory<SurfaceMeshType> my_factory(&data);
            // Create our mesh
            SurfaceMeshType* surface_mesh(new SurfaceMeshType(my_factory));

            return new DiscreteChart(surface_mesh);
          }

          /** \copydoc ChartBase::project()
           *
           * \warning Dumb, inefficient and imprecise routine for testing purposes.
           *
           * \todo: Implement something better
           */
          void project(WorldPoint& point) const
          {
            Index min_index(0);
            DataType min_dist(Math::Limits<DataType>::max());

            typename SurfaceMeshType::VertexSetType& vtx(_surface_mesh->get_vertex_set());
            DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.5)));

            for(Index i(0); i < _surface_mesh->get_num_entities(0); ++i)
            {
              typename SurfaceMeshType::VertexSetType::VertexType tmp(point);
              tmp -= vtx[i];

              DataType current_dist(tmp.norm_euclid());

              if(current_dist <= tol)
              {
                min_index = i;
                break;
              }

              if(current_dist < min_dist)
              {
                min_dist = current_dist;
                min_index = i;
              }

            }

            point = vtx[min_index];

          } // void project()

      }; // class DiscreteChart

    } // namespace Atlas

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_ATLAS_DISCRETE_CHART_HPP
