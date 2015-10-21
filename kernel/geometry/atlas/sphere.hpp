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

        /** \copydoc ChartBase::write_data_container */
        virtual void write_data_container(MeshStreamer::ChartContainer& chart_data) const override
        {
          chart_data.data.push_back(" <sphere>");
          chart_data.data.push_back("  radius  "+stringify(stringify_fp_sci(_radius)));
          chart_data.data.push_back("  midpoint "+stringify(stringify_fp_sci(_midpoint(0)))+" "+stringify(stringify_fp_sci(_midpoint(1)))+" "+stringify(stringify_fp_sci(_midpoint(2))));
          chart_data.data.push_back(" </sphere>");
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
        void project(Mesh_& mesh, const MeshPart<Mesh_>& meshpart) const
        {
          auto& vtx = mesh.get_vertex_set();
          const auto& target_vtx = meshpart.template get_target_set<0>();

          for(Index i(0); i < meshpart.get_num_entities(0); ++i)
          {
            project(reinterpret_cast<WorldPoint&>(vtx[target_vtx[i]]));
          }
        }

        /**
         * \brief Builds a Sphere<Mesh> object from parsed data
         *
         * \param[in] data
         * Parameters for the object in the form of Strings
         *
         * \param[in] line
         * The line in which the sphere data started in the original file
         *
         * \returns
         * A pointer to the new object.
         */
        static Sphere<Mesh_>* parse(const std::deque<String>& data, const Index line)
        {
          bool have_midpoint(false), have_radius(false);
          DataType mid_x(0), mid_y(0), mid_z(0), radius(0);
          std::deque<String> lines(data), dat_line;

          if(lines.front() != "<sphere>")
            throw InternalError("Expected '<sphere>' but found '" + lines.front() + "' in line " + stringify(line));
          lines.pop_front();

          // loop over all data chunk lines
          for(Index lno(line+1); !lines.empty(); ++lno)
          {
            String l = lines.front().trim();
            lines.pop_front();

            if(l.empty())
              continue;

            if(l == "</sphere>")
              break;

            // split the line by whitespaces
            l.split_by_charset(dat_line);
            if(dat_line[0] == "radius")
            {
              if(dat_line.size() != 2)
                throw InternalError("Invalid number of arguments in line " + stringify(lno));

              // try to parse radius
              have_radius = dat_line[1].parse(radius);
              if(!have_radius)
                throw InternalError("Failed to parse '" + dat_line[1] + "' as sphere radius in line " + stringify(lno));

              // okay
              continue;
            }
            if(dat_line[0] == "midpoint")
            {
              if(dat_line.size() != 4)
                throw InternalError("Invalid number of arguments in line " + stringify(lno));

              // try to parse midpoint
              have_midpoint = dat_line[1].parse(mid_x) && dat_line[2].parse(mid_y) && dat_line[3].parse(mid_z);
              if(!have_midpoint)
                throw InternalError("Failed to parse '" + l + "' as sphere midpoint in line " + stringify(lno));

              // okay
              continue;
            }

            // something invalid
            throw InternalError("Unexpected '" + l + "' in line " + stringify(lno));
          }

          if(!have_midpoint)
            throw  InternalError("Sphere chart has no midpoint!");
          if(!have_radius)
            throw  InternalError("Sphere chart has no radius!");

          // okay, create sphere
          return new Sphere<Mesh_>(mid_x, mid_y, mid_z, radius);
        }
      };
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_SPHERE_HPP
