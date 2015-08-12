#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_TUBE_HPP
#define KERNEL_GEOMETRY_ATLAS_TUBE_HPP 1

#include <kernel/geometry/atlas/chart.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Atlas
    {
      /// Tube chart traits
      struct TubeTraits
      {
        /// we support explicit map
        static constexpr bool is_explicit = false;
        /// we support implicit projection
        static constexpr bool is_implicit = true;
        /// this is a 2D object
        static constexpr int world_dim = 3;
        /// we have 1D parameters
        static constexpr int param_dim = 2;
      };

      /**
       * \brief Tube chart class template
       *
       * This class represents a 3D tube (aka the curved boundary part of a cylinder)
       * characterised by a midpoint, a rotation axis and a radius.
       *
       * This chart implements only the implicit chart interface.
       *
       * \tparam Mesh_
       * The type of the mesh to be parameterised by this chart.
       *
       * \author Peter Zajac
       */
      template<typename Mesh_>
      class Tube :
        public ChartCRTP<Tube<Mesh_>, Mesh_, TubeTraits>
      {
      public:
        typedef ChartCRTP<Tube<Mesh_>, Mesh_, TubeTraits> BaseClass;
        typedef typename BaseClass::CoordType DataType;
        typedef typename BaseClass::WorldPoint WorldPoint;
        typedef typename BaseClass::ParamPoint ParamPoint;

      protected:
        /// the tube's midpoint
        WorldPoint _midpoint;
        /// the tube's rotation axis
        WorldPoint _rot_axis;
        /// the tube's rotation axis inverse length
        DataType _rot_scale;
        /// the tube's radius
        DataType _radius;
        /// trafo coefficients
        DataType _trafo_x0, _trafo_x1, _trafo_y0, _trafo_y1;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] mid_x, mid_y, mid_z
         * The coordinates of the tube midpoint.
         *
         * \param[in] axis_x, axis_y, axis_z
         * The components of the rotation axis vector. At least one of those three parameters must
         * be non-zero, but the axis vector need not be normalised.
         *
         * \param[in] radius
         * The radius of the tube. Must be positive.
         */
        explicit Tube(
          DataType mid_x, DataType mid_y, DataType mid_z,
          DataType axis_x, DataType axis_y, DataType axis_z,
          DataType radius) :
          _radius(radius)
        {
          ASSERT_(radius > DataType(0));
          _midpoint[0] = mid_x;
          _midpoint[1] = mid_y;
          _midpoint[2] = mid_z;
          _rot_axis[0] = axis_x;
          _rot_axis[1] = axis_y;
          _rot_axis[2] = axis_z;
          // get inverse rotation axis length
          _rot_scale = (DataType(1) / Math::sqr(_rot_axis.norm_euclid()));
        }

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "tube";
        }

        /** \copydoc ChartBase::write_data_container */
        virtual void write_data_container(MeshStreamer::ChartContainer& chart_data) const override
        {
          chart_data.data.push_back(" <tube>");
          chart_data.data.push_back("  radius  "+stringify(scientify(_radius)));
          chart_data.data.push_back("  midpoint "+stringify(scientify(_midpoint(0)))+" "+stringify(scientify(_midpoint(1)))+" "+stringify(scientify(_midpoint(2))));
          chart_data.data.push_back("  axis"+stringify(scientify(_rot_axis(0)))+" "+stringify(scientify(_rot_axis(1)))+" "+stringify(scientify(_rot_axis(2))));
          chart_data.data.push_back(" </tube>");
        }

        /** \copydoc ChartBase::project() */
        void project(WorldPoint& point) const
        {
          // subtract tube midpoint
          point -= _midpoint;
          // project point onto axis
          WorldPoint axis_point = (_rot_scale * Tiny::dot(_rot_axis, point)) * _rot_axis;
          // compute difference of point and axis point
          WorldPoint diff_point = point - axis_point;
          // compute projected point
          point = (_midpoint + axis_point) + diff_point * (_radius / diff_point.norm_euclid());
        }

        static Tube<Mesh_>* parse(const std::deque<String>& data, const Index line)
        {
          bool have_midpoint(false), have_axis(false), have_radius(false);
          DataType mid_x(0), mid_y(0), mid_z(0), axis_x(0), axis_y(0), axis_z(0), radius(0);
          std::deque<String> lines(data), dat_line;

          if(lines.front() != "<tube>")
            throw InternalError("Expected '<tube>' but found '" + lines.front() + "' in line " + stringify(line));
          lines.pop_front();

          // loop over all data chunk lines
          for(Index lno(line+1); !lines.empty(); ++lno)
          {
            String l = lines.front().trim();
            lines.pop_front();

            if(l.empty())
              continue;

            if(l == "</tube>")
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
                throw InternalError("Failed to parse '" + dat_line[1] + "' as tube radius in line " + stringify(lno));

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
                throw InternalError("Failed to parse '" + l + "' as tube midpoint in line " + stringify(lno));

              // okay
              continue;
            }
            if(dat_line[0] == "axis")
            {
              if(dat_line.size() != 4)
                throw InternalError("Invalid number of arguments in line " + stringify(lno));

              // try to parse axis
              have_axis = dat_line[1].parse(axis_x) && dat_line[2].parse(axis_y) && dat_line[3].parse(axis_z);
              if(!have_axis)
                throw InternalError("Failed to parse '" + l + "' as tube axis in line " + stringify(lno));

              // okay
              continue;
            }

            // something invalid
            throw InternalError("Unexpected '" + l + "' in line " + stringify(lno));
          }

          if(!have_midpoint)
            throw  InternalError("Tube chart has no midpoint!");
          if(!have_radius)
            throw  InternalError("Tube chart has no radius!");
          if(!have_axis)
            throw  InternalError("Tube chart has no axis!");

          // okay, create tube
          return new Tube<Mesh_>(mid_x, mid_y, mid_z, axis_x, axis_y, axis_z, radius);
        }
      }; // class Tube<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_TUBE_HPP
