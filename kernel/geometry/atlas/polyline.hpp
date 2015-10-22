#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_POLYLINE_HPP
#define KERNEL_GEOMETRY_ATLAS_POLYLINE_HPP 1

#include <kernel/geometry/atlas/chart.hpp>

#include <deque>

namespace FEAST
{
  namespace Geometry
  {
    namespace Atlas
    {
      /// Polyline chart traits
      struct PolylineTraits
      {
        /// we support explicit map
        static constexpr bool is_explicit = true;
        /// we don't support implicit project
        static constexpr bool is_implicit = false;
        /// this is a 2D object
        static constexpr int world_dim = 2;
        /// we have 1D parameters
        static constexpr int param_dim = 1;
      };

      /**
       * \brief Polyline chart class template
       *
       * This chart represents a 2D polygonal chain (aka "polyline") characterised
       * by vector of parameter and world point coordinates.
       *
       * This chart currently implements only the explicit chart interface.
       *
       * \tparam Mesh_
       * The type of the mesh to be parameterised by this chart.
       */
      template<typename Mesh_>
      class Polyline :
        public ChartCRTP<Polyline<Mesh_>, Mesh_, PolylineTraits>
      {
      public:
        /// The CRTP base class
        typedef ChartCRTP<Polyline<Mesh_>, Mesh_, PolylineTraits> BaseClass;
        /// Floating point type for coordinates
        typedef typename BaseClass::CoordType DataType;
        /// Vector type for world points aka. image points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Vector type for parameter points aka. domain points
        typedef typename BaseClass::ParamPoint ParamPoint;

      protected:
        /// The world points making up the polygon line
        std::deque<WorldPoint> _world;
        /// The parameter values for these world points
        std::deque<ParamPoint> _param;

      public:
        /// default CTOR
        explicit Polyline()
        {
        }

        /**
         * \brief Constructor.
         *
         * \param[in] world
         * A deque of world points representing the polyline.
         *
         * \param[in] param
         * A deque of parameter points representing the parameter space.
         */
        explicit Polyline(const std::deque<WorldPoint>& world, const std::deque<ParamPoint>& param) :
          _world(world),
          _param(param)
        {
          ASSERT_(_param.size() == _world.size());
          // we need at least 2 points
          ASSERT_(_param.size() > std::size_t(1));
        }

        /**
         * \brief Maps a single parameter point
         *
         * \param[out] point
         * The image of the parameter point under the chart mapping
         *
         * \param[in] param
         * The parameter point to be mapped
         *
         */
        void map(WorldPoint& point, const ParamPoint& param) const
        {
          // get our point coord
          const DataType& x = param[0];

          // check for boundary
          if(x <= _param.front()[0])
          {
            point = _world.front();
            return;
          }
          if(x >= _param.back()[0])
          {
            point = _world.back();
            return;
          }

          /// \todo replace by binary search
          // apply linear search
          for(std::size_t i(0); (i+1) < _param.size(); ++i)
          {
            // fetch left and right interval ends
            const DataType& xl = _param[i  ][0];
            const DataType& xr = _param[i+1][0];

            // is this our interval?
            if((xl <= x) && (x <= xr))
            {
              // compute normalised interpolation factor
              DataType t1 = (x - xl) / (xr - xl);
              DataType t0 = DataType(1) - t1;

              // interpolate world point
              point = (t0 * _world[i+0]) + (t1 * _world[i+1]);
              return;
            }
          }
        }

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "polyline";
        }

        /** \copydoc ChartBase::write_data_container */
        virtual void write_data_container(MeshStreamer::ChartContainer& chart_container) const override
        {

          size_t num_points(_world.size());

          // GCC 4.9.2 complains at link time in the img_dim line otherwise
          int img_dim(BaseClass::world_dim);

          chart_container.data.push_back(" <polyline>");
          chart_container.data.push_back("  img_dim "+stringify(img_dim));
          chart_container.data.push_back("  num_points "+stringify(num_points));
          chart_container.data.push_back("  <points>");

          for(size_t i(0); i < num_points; ++i)
          {
            String tmp("  ");

            for(int j(0); j < BaseClass::param_dim; ++j)
              tmp += " "+stringify((_param[i])[j]);

            // GCC 4.9.2 does not complain here, though
            for(int j(0); j < BaseClass::world_dim; ++j)
              tmp += " "+stringify((_world[i])[j]);

            chart_container.data.push_back(tmp);
          }

          chart_container.data.push_back("  </points>");
          chart_container.data.push_back(" </polyline>");
        }

        /**
         * \brief Builds a Polyline<Mesh> object from parsed data
         *
         * \param[in] data
         * Parameters for the object in the form of Strings
         *
         * \param[in] line
         * The line in which the polyline data started in the original file
         *
         * \returns
         * A pointer to the new object.
         */
        static Polyline<Mesh_>* parse(const std::deque<String>& data, const Index line)
        {
          bool have_img_dim(false), have_num_points(false);
          std::deque<String> lines(data), dat_line;
          int img_dim = 0;
          int num_points = 0;

          if(lines.front() != "<polyline>")
            throw InternalError("Invalid marker");
          lines.pop_front();

          // loop over all data chunk lines
          Index lno(line+1);
          for(; !lines.empty(); ++lno)
          {
            // get line
            String l = lines.front().trim();

            // is it the <points> marker?
            if(l == "<points>")
              break;

            // pop line from deque
            lines.pop_front();

            // trim line and check for emptyness
            if(l.empty())
              continue;

            // split line
            l.split_by_charset(dat_line);

            if(dat_line[0] == "img_dim")
            {
              if(dat_line.size() != 2)
                throw InternalError("Invalid number of arguments");

              // try to parse image dimension
              have_img_dim = dat_line[1].parse(img_dim);
              if(!have_img_dim)
                throw InternalError("Failed to parse '" + dat_line[1] + "' as image dimension");

              continue;
            }

            if(dat_line[0] == "num_points")
            {
              if(dat_line.size() != 2)
                throw InternalError("Invalid number of arguments");

              // try to parse radius
              have_num_points = dat_line[1].parse(num_points);
              if(!have_num_points)
                throw InternalError("Failed to parse '" + dat_line[1] + "' as point count");

              continue;
            }

            // something invalid
            throw InternalError("Unrecognised line '" + l + "'");
          }

          if(!have_img_dim)
            throw  InternalError("Polyline has no image dimension!");
          if(!have_num_points)
            throw  InternalError("Polyline has no point count!");

          // ensure that we have at least 2 points
          if(num_points < 2)
            throw InternalError("Invalid number of points for Polyline!");

          // create the polyline
          switch(img_dim)
          {
          /*case 1:
            _parse_polyline1d_points(lines, num_points);
            break;*/
          case 2:
            return _parse_polyline2d_points(lines, lno, num_points);
          /*case 3:
            _parse_polyline3d_points(lines, num_points);
            break;*/
          default:
            throw InternalError("Invalid image dimension for Polyline!");
          }

          return nullptr;
        }

      protected:
        /**
         * \brief Parses the 'points' section of a 2D 'polyline' chart
         */
        static Polyline<Mesh_>* _parse_polyline2d_points(std::deque<String>& lines, const Index line, const int num_points)
        {
          // get front line
          if(lines.front() != "<points>")
            throw InternalError("Expected '<points>' but found '" + lines.front() + "'");
          lines.pop_front();

          std::deque<String> dat_line;
          std::deque<WorldPoint> world;
          std::deque<ParamPoint> param;
          world.resize(std::size_t(num_points));
          param.resize(std::size_t(num_points));

          // now loop over all lines
          for(Index lno(line+1), i(0); !lines.empty(); ++lno)
          {
            String l = lines.front().trim();
            lines.pop_front();
            if(l.empty())
              continue;

            if(l == "</points>")
              break;

            // split line
            l.split_by_charset(dat_line);

            if(dat_line.size() != 3)
              throw InternalError("Invalid number of coordinates in line " + stringify(lno));

            // parse parameters
            if(!dat_line[0].parse(param[i][0]))
              throw InternalError("Failed to parse '" + dat_line[0] + "' as parameter value in line " + stringify(lno));
            if(!dat_line[1].parse(world[i][0]))
              throw InternalError("Failed to parse '" + dat_line[1] + "' as x-coordinate in line " + stringify(lno));
            if(!dat_line[2].parse(world[i][1]))
              throw InternalError("Failed to parse '" + dat_line[2] + "' as y-coordinate in line " + stringify(lno));

            // push new point
            ++i;
          }

          // check point count
          if(int(param.size()) != num_points)
            throw InternalError("Parsed " + stringify(param.size()) + " points but expected " + stringify(num_points));

          // alright, create the chart
          return new Polyline<Mesh_>(world, param);
        }
      }; // class Polyline<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_POLYLINE_HPP
