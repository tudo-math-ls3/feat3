#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_SPHERE_HPP
#define KERNEL_GEOMETRY_ATLAS_SPHERE_HPP 1

#include <kernel/geometry/atlas/chart.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Atlas
    {
      /// Sphere chart traits
      struct SphereTraits
      {
        /// we support explicit map
        static constexpr bool is_explicit = false;
        /// we support implicit projection
        static constexpr bool is_implicit = true;
        /// this is a 2D object
        static constexpr int world_dim = 3;
        /// we have 1D parameters
        static constexpr int param_dim = 1;
      };

      /**
       * \brief Sphere chart class template.
       *
       * This class represents a 3D sphere characterised by a midpoint and a radius.
       *
       * This chart implements only the implicit chart interface.
       *
       * \tparam Mesh_
       * The type of the mesh to be parameterised by this chart.
       *
       * \author Peter Zajac
       */
      template<typename Mesh_>
      class Sphere :
        public ChartCRTP<Sphere<Mesh_>, Mesh_, SphereTraits>
      {
      public:
        /// CRTP base class
        typedef ChartCRTP<Sphere<Mesh_>, Mesh_, SphereTraits> BaseClass;
        /// Floating point type
        typedef typename BaseClass::CoordType DataType;
        /// Vector type for world points, aka image points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Vector type for parameter points, aka domain points
        typedef typename BaseClass::ParamPoint ParamPoint;

      protected:
        /// the sphere's midpoint
        WorldPoint _midpoint;
        /// the sphere's radius
        DataType _radius;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] mid_x, mid_y, mid_z
         * The coordinates of the sphere midpoint.
         *
         * \param[in] radius
         * The radius of the sphere. Must be positive.
         */
        explicit Sphere(DataType mid_x, DataType mid_y, DataType mid_z, DataType radius) :
          _radius(radius)
        {
          ASSERT_(radius > DataType(0));
          _midpoint[0] = mid_x;
          _midpoint[1] = mid_y;
          _midpoint[2] = mid_z;
        }

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "sphere";
        }

        /**
         * \brief Projects a single world point
         *
         * \param[in,out] point
         * The world point to be projected
         *
         */
        void project(WorldPoint& point) const
        {
          point -= _midpoint;
          DataType dist = point.norm_euclid();
          point *= (_radius / dist);
          point += _midpoint;
        }

        /**
         * \brief Projects all mesh points identified by a meshpart
         *
         * \param[in,out] mesh
         * The mesh whose points will be projected
         *
         * \param[in] meshpart
         * The MeshPart identifying the point to be projected
         *
         */
        void project_meshpart(Mesh_& mesh, const MeshPart<Mesh_>& meshpart) const
        {
          auto& vtx = mesh.get_vertex_set();
          const auto& target_vtx = meshpart.template get_target_set<0>();

          for(Index i(0); i < meshpart.get_num_entities(0); ++i)
          {
            project(reinterpret_cast<WorldPoint&>(vtx[target_vtx[i]]));
          }
        }

        /** \copydoc ChartBase::write */
        virtual void write(std::ostream& os, const String& sindent) const override
        {
          os << sindent << "<Sphere";
          os << " radius=\"" << this->_radius << "\"";
          os << " midpoint=\"" << this->_midpoint[0] << " " << this->_midpoint[1] << " " << this->_midpoint[2] << "\"";
          os << " />" << std::endl;
        }
      };
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_SPHERE_HPP
