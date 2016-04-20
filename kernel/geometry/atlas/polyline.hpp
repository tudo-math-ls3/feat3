#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_POLYLINE_HPP
#define KERNEL_GEOMETRY_ATLAS_POLYLINE_HPP 1

#include <kernel/geometry/atlas/spline_base.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Atlas
    {
      /**
       * \brief Polyline chart class template
       *
       * This chart represents a 2D polygonal chain (aka "polyline") characterised
       * by vector of parameter and world point coordinates.
       *
       * This chart implements both the implicit and explicit chart interfaces.
       *
       * \tparam Mesh_
       * The type of the mesh to be parameterised by this chart.
       *
       * \author Peter Zajac
       */
      template<typename Mesh_>
      class Polyline :
        public SplineBaseCRTP<Polyline<Mesh_>, Mesh_>
      {
      public:
        /// The CRTP base class
        typedef SplineBaseCRTP<Polyline<Mesh_>, Mesh_> BaseClass;
        /// Floating point type for coordinates
        typedef typename BaseClass::CoordType DataType;
        /// Vector type for world points aka. image points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Vector type for parameter points aka. domain points
        typedef typename BaseClass::ParamPoint ParamPoint;

      public:
        /// default CTOR
        explicit Polyline(bool closed = false) :
          BaseClass(closed)
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
         *
         * \param[in] closed
         * Specifies whether the polyline is closed.
         */
        explicit Polyline(const std::deque<WorldPoint>& world, const std::deque<ParamPoint>& param, bool closed) :
          BaseClass(closed)
        {
          this->_world = world;
          this->_param = param;

          // we need at least 2 points
          ASSERT_(this->_world.size() > std::size_t(1));
          ASSERT_(this->_param.empty() || (this->_param.size() == this->_world.size()));
        }

        /**
         * \brief Maps a local segment parameter point
         *
         * \param[in] i
         * The segment index onto which to map
         *
         * \param[in] t
         * The local segment parameter
         */
        WorldPoint map_on_segment(const Index i, const DataType t) const
        {
          ASSERT_((i+1) < Index(this->_world.size()));
          return ((DataType(1) - t) * this->_world[i]) + (t * this->_world[i+1]);
        }

        /**
         * \brief Projects a point onto one segment of the polyline.
         *
         * \param[in] i
         * The index of the polyline segment onto which to project
         *
         * \param[in] point
         * The world point that is to be projected.
         *
         * \returns
         * The projected world point parameter.
         */
        DataType project_on_segment(const Index i, const WorldPoint& point) const
        {
          ASSERT_((i+1) < Index(this->_world.size()));

          // get the ends of the line segment
          const WorldPoint& x0 = this->_world.at(i);
          const WorldPoint& x1 = this->_world.at(i+1);

          // compute xe := x1-x0 and xp := p-x0
          WorldPoint xe = x1 - x0;
          WorldPoint xp = point - x0;

          // compute t = <xp,xe>/<xe,xe> and clamp it to [0,1]
          return Math::clamp(Tiny::dot(xp,xe) / Tiny::dot(xe,xe), DataType(0), DataType(1));
        }

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "polyline";
        }

        /** \copydoc ChartBase::write */
        virtual void write(std::ostream& os, const String& sindent) const override
        {
          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          os << sindent << "<Polyline dim=\"2\" size=\"" << this->_world.size() << "\"";
          os << " type=\"" << (this->_closed ? "closed" : "open") << "\">" << std::endl;
          this->write_points(os, sind);
          this->write_params(os, sind);
          os << sindent << "</Polyline>" << std::endl;
        }

      protected:
        /**
         * \brief Parses the 'points' section of a 2D 'polyline' chart
         */
        static Polyline<Mesh_>* _parse_polyline2d_points(
          std::deque<String>& lines, const Index line, const int num_points)
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
          return new Polyline<Mesh_>(world, param, false);
        }
      }; // class Polyline<...>

      template<typename Mesh_, typename ChartReturn_ = ChartBase<Mesh_>>
      class PolylineChartParser :
        public Xml::MarkupParser
      {
      private:
        typedef Polyline<Mesh_> ChartType;
        typedef typename ChartType::DataType DataType;
        ChartReturn_*& _chart;
        Polyline<Mesh_>* _polyline;
        Index _size;

      public:
        explicit PolylineChartParser(ChartReturn_*& chart) :
          _chart(chart),
          _polyline(nullptr),
          _size(0)
        {
        }

        virtual bool attribs(std::map<String,bool>& attrs) const override
        {
          attrs.emplace("dim", true);
          attrs.emplace("size", true);
          attrs.emplace("type", false);
          return true;
        }

        virtual void create(
          int iline,
          const String& sline,
          const String&,
          const std::map<String, String>& attrs,
          bool closed) override
        {
          // make sure this one isn't closed
          if(closed)
            throw Xml::GrammarError(iline, sline, "Invalid closed Polyline markup");

          Index dim(0);

          // try to parse the dimension
          if(!attrs.find("dim")->second.parse(dim))
            throw Xml::GrammarError(iline, sline, "Failed to parse polyline dimension");
          if(dim != Index(2))
            throw Xml::GrammarError(iline, sline, "Invalid polyline dimension");

          // try to parse the size
          if(!attrs.find("size")->second.parse(_size))
            throw Xml::GrammarError(iline, sline, "Failed to parse polyline size");
          if(_size < Index(2))
            throw Xml::GrammarError(iline, sline, "Invalid polyline size");

          // try to check type
          bool poly_closed(false);
          auto it = attrs.find("type");
          if(it != attrs.end())
          {
            String stype = it->second;
            if(it->second == "closed")
              poly_closed = true;
            else if (it->second != "open")
              throw Xml::ContentError(iline, sline, "Invalid polyline type; must be either 'closed' or 'open'");
          }

          // up to now, everything's fine
          _polyline = new ChartType(poly_closed);
        }

        virtual void close(int, const String&) override
        {
          // okay
          _chart = _polyline;
        }

        virtual bool content(int, const String&) override
        {
          return false;
        }

        virtual std::shared_ptr<Xml::MarkupParser> markup(int, const String&, const String& name) override
        {
          if(name == "Points") return std::make_shared<SplinePointsParser<Polyline<Mesh_>>>(*_polyline, _size);
          if(name == "Params") return std::make_shared<SplineParamsParser<Polyline<Mesh_>>>(*_polyline, _size);
          return nullptr;
        }
      }; // class PolylineChartParser<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_POLYLINE_HPP
