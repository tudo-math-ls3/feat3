#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_CIRCLE_HPP
#define KERNEL_GEOMETRY_ATLAS_CIRCLE_HPP 1

#include <kernel/geometry/atlas/chart.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Atlas
    {
      /// Circle chart traits
      struct CircleTraits
      {
        /// we support explicit map
        static constexpr bool is_explicit = true;
        /// we support implicit projection
        static constexpr bool is_implicit = true;
        /// this is a 2D object
        static constexpr int world_dim = 2;
        /// we have 1D parameters
        static constexpr int param_dim = 1;
      };

      /**
       * \brief Circle chart class template.
       *
       * This class represents a 2D circle characterised by a midpoint, a radius and a 1D interval
       * for the parameter space.
       *
       * This chart implements both the implicit and explicit chart interfaces.
       *
       * The parameter domain [a,b] is mapped onto the circle, where both interval endpoints
       * are mapped onto the east point on the circle, i.e. the point with the X/Y-coordinates
       * [mid_x + radius, mid_y]. If the parameter interval is oriented positive, i.e. if it
       * holds that a < b, then the circle is parameterised in counter-clock-wise orientation.
       * If the parameter interval is oriented negative, i.e. if it holds that b < a, then
       * the circle is parameterised in clock-wise orientation.
       *
       * \tparam Mesh_
       * The type of the mesh to be parameterised by this chart.
       *
       * \author Peter Zajac
       */
      template<typename Mesh_>
      class Circle :
        public ChartCRTP<Circle<Mesh_>, Mesh_, CircleTraits>
      {
      public:
        typedef ChartCRTP<Circle<Mesh_>, Mesh_, CircleTraits> BaseClass;
        typedef typename BaseClass::CoordType DataType;
        typedef typename BaseClass::WorldPoint WorldPoint;
        typedef typename BaseClass::ParamPoint ParamPoint;

      protected:
        /// the circle's midpoint
        WorldPoint _mid_point;
        /// the circle's radius
        DataType _radius;
        /// domain transformation coefficients
        DataType _trafo_a, _trafo_b;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] mid_x, mid_y
         * The coordinates of the circle midpoint.
         *
         * \param[in] radius
         * The radius of the circle. Must be positive.
         *
         * \param[in] param_l, param_r
         * Specifies the domain of the parameter interval.
         * See the documentation of this class template for details.
         */
        explicit Circle(DataType mid_x, DataType mid_y, DataType radius,
          DataType param_l = DataType(0), DataType param_r = DataType(1)) :
          _radius(radius),
          _trafo_a(-param_l),
          _trafo_b((DataType(2) * Math::pi<DataType>()) / (param_r - param_l))
        {
          ASSERT_(radius > DataType(0));
          _mid_point[0] = mid_x;
          _mid_point[1] = mid_y;
        }

        /** \copydoc ChartBase::project() */
        void project(WorldPoint& point) const
        {
          point -= _mid_point;
          DataType dist = point.norm_euclid();
          point *= (_radius / dist);
          point += _mid_point;
        }

        /** \copydoc ChartBase::map() */
        void map(WorldPoint& point, const ParamPoint& param) const
        {
          // transform paramter to interval [0, 2*pi)
          DataType x = (param[0] + _trafo_a) * _trafo_b;
          point[0] = this->_mid_point[0] + this->_radius * Math::cos(x);
          point[1] = this->_mid_point[1] + this->_radius * Math::sin(x);
        }

        static Circle<Mesh_>* parse(const std::deque<String>& data, const Index line)
        {
          bool have_midpoint(false), have_radius(false), have_domain(false);
          DataType mid_x(0), mid_y(0), dom_l(0), dom_r(0), radius(0);
          std::deque<String> lines(data), dat_line;

          if(lines.front() != "<circle>")
            throw InternalError("Expected '<circle>' but found '" + lines.front() + "' in line " + stringify(line));
          lines.pop_front();

          // loop over all data chunk lines
          for(Index lno(line+1); !lines.empty(); ++lno)
          {
            String l = lines.front().trim();
            lines.pop_front();

            if(l.empty())
              continue;

            if(l == "</circle>")
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
                throw InternalError("Failed to parse '" + dat_line[1] + "' as circle radius in line " + stringify(lno));

              // okay
              continue;
            }
            if(dat_line[0] == "midpoint")
            {
              if(dat_line.size() != 3)
                throw InternalError("Invalid number of arguments in line " + stringify(lno));

              // try to parse midpoint
              have_midpoint = dat_line[1].parse(mid_x) && dat_line[2].parse(mid_y);
              if(!have_midpoint)
                throw InternalError("Failed to parse '" + dat_line[1] + " " + dat_line[2] + "' as circle midpoint in line " + stringify(lno));

              // okay
              continue;
            }
            if(dat_line[0] == "domain")
            {
              if(dat_line.size() != 3)
                throw InternalError("Invalid number of arguments in line " + stringify(lno));

              have_domain = dat_line[1].parse(dom_l) && dat_line[2].parse(dom_r);
              if(!have_domain)
                throw InternalError("Failed to parse '" + dat_line[1] + " " + dat_line[2] + "' as circle domain in line " + stringify(lno));

              // okay
              continue;
            }

            // something invalid
            throw InternalError("Unexpected '" + l + "' in line " + stringify(lno));
          }

          if(!have_midpoint)
            throw  InternalError("Circle chart has no midpoint!");
          if(!have_radius)
            throw  InternalError("Circle chart has no radius!");
          if(!have_domain)
            throw  InternalError("Circle chart has no domain!");

          // okay, create circle
          return new Circle<Mesh_>(mid_x, mid_y, radius, dom_l, dom_r);
        }
      }; // class Circle<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_CIRCLE_HPP