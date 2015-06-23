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
        typedef ChartCRTP<Sphere<Mesh_>, Mesh_, SphereTraits> BaseClass;
        typedef typename BaseClass::CoordType DataType;
        typedef typename BaseClass::WorldPoint WorldPoint;
        typedef typename BaseClass::ParamPoint ParamPoint;

      protected:
        /// the circle's midpoint
        WorldPoint _mid_point;
        /// the circle's radius
        DataType _radius;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] mid_x, mid_y, mid_z
         * The coordinates of the sphere  midpoint.
         *
         * \param[in] radius
         * The radius of the sphere. Must be positive.
         */
        explicit Sphere(DataType mid_x, DataType mid_y, DataType mid_z, DataType radius) :
          _radius(radius)
        {
          ASSERT_(radius > DataType(0));
          _mid_point[0] = mid_x;
          _mid_point[1] = mid_y;
          _mid_point[2] = mid_z;
        }

        /** \copydoc ChartBase::project() */
        void project(WorldPoint& point) const
        {
          point -= _mid_point;
          DataType dist = point.norm_euclid();
          point *= (_radius / dist);
          point += _mid_point;
        }

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

          // okay, create tube
          return new Sphere<Mesh_>(mid_x, mid_y, mid_z, radius);
        }
      };
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_SPHERE_HPP
