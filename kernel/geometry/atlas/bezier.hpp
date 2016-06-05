#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_BEZIER_HPP
#define KERNEL_GEOMETRY_ATLAS_BEZIER_HPP 1

#include <kernel/geometry/atlas/spline_base.hpp>

namespace FEAT
{
  namespace Geometry
  {
    namespace Atlas
    {
      /**
       * \brief Bezier chart class template
       *
       * This chart represents a 2D cubic Bezier spline characterised
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
      class Bezier :
        public SplineBaseCRTP<Bezier<Mesh_>, Mesh_>
      {
      public:
        /// The CRTP base class
        typedef SplineBaseCRTP<Bezier<Mesh_>, Mesh_> BaseClass;
        /// Floating point type for coordinates
        typedef typename BaseClass::CoordType DataType;
        /// Vector type for world points aka. image points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Vector type for parameter points aka. domain points
        typedef typename BaseClass::ParamPoint ParamPoint;

      protected:
        /// The control points for the Bezier curve
        std::deque<WorldPoint> _control;

      public:
        /// default CTOR
        explicit Bezier(bool closed = false) :
          BaseClass(closed),
          _control()
        {
        }

        void push_control(const WorldPoint& xc1, const WorldPoint& xc2)
        {
          _control.push_back(xc1);
          _control.push_back(xc2);
        }

        std::deque<WorldPoint>& get_control_points()
        {
          return _control;
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
          ASSERT((i+1) < Index(this->_world.size()));

          // compute normalised interpolation factor
          const DataType s = DataType(1) - t;

          // references to the four Bezier points
          const WorldPoint& p0 = this->_world[i];
          const WorldPoint& p1 = this->_control[2*i];
          const WorldPoint& p2 = this->_control[2*i+1];
          const WorldPoint& p3 = this->_world[i+1];

          // interpolated points
          const WorldPoint p01 = (s*p0)  + (t*p1);
          const WorldPoint p12 = (s*p1)  + (t*p2);
          const WorldPoint p23 = (s*p2)  + (t*p3);
          const WorldPoint p02 = (s*p01) + (t*p12);
          const WorldPoint p13 = (s*p12) + (t*p23);

          // interpolate world point
          return (s * p02) + (t *p13);
        }

        /**
         * \brief Projects a point onto one segment of the Bezier spline.
         *
         * \param[in] i
         * The index of the Bezier segment onto which to project
         *
         * \param[in] point
         * The world point that is to be projected.
         *
         * \returns
         * The projected world point parameter.
         */
        DataType project_on_segment(const Index i, const WorldPoint& point) const
        {
          ASSERT((i+1) < Index(this->_world.size()));

          // references to the four Bezier points
          const WorldPoint& p0 = this->_world[i];
          const WorldPoint& p1 = this->_control[2*i];
          const WorldPoint& p2 = this->_control[2*i+1];
          const WorldPoint& p3 = this->_world[i+1];

          //
          // Algorithm description
          // ---------------------
          // This function is meant to find the parameter value 't' of the point
          // B(t) closest to the given input point 'X' (named 'point'), i.e.
          //
          //                 argmin     { (B(t) - X)^2 }
          //             {0 <= t <= 1}
          //
          // As Bezier curves a cubic polynomials (and therefore smooth), the
          // parameter 't' we are looking for fulfills the orthogonal projection
          // property
          //
          //              f(t) :=   < B'(t), B(t) - X > = 0
          //
          // The algorithm implemented below is a simple Newton iteration applied
          // onto the non-linear (fifth degree polynomial) equation above.
          // For Newton, we require the derivative f'(t) of f(t), which is
          //
          //           f'(t) = < B"(t), B(t) - X > + < B'(t), B'(t) >
          //
          // and the Newton iteration is defined as
          //
          //               t_{k+1} := t_k - f(t_k) / f'(t_k)
          //
          // One problem is to find an appropriate initial guess t_0 for our
          // algorithm. In many cases, the initial guess t_0 = 1/2 will do the
          // job. However, if the Bezier curve is S-shaped, this algorithm may
          // fail, therefore we have to choose another initial guess.
          // To circumvent this problem, we explicitly test the distance
          // (B(t) - X)^2 for t in {1/5, 1/2, 4/5} and choose the t with
          // the minimal distance as an initial guess.
          //
          // Unfortunately, there is no guarantee that this algorithm will
          // converge.

          // choose an appropriate initial parameter value
          DataType t = DataType(0.5);
          {
            DataType d1 = (point - map_on_segment(i, DataType(0.2))).norm_euclid_sqr();
            DataType d2 = (point - map_on_segment(i, DataType(0.5))).norm_euclid_sqr();
            DataType d3 = (point - map_on_segment(i, DataType(0.8))).norm_euclid_sqr();
            if(d1 < Math::min(d2,d3)) t = DataType(0.2);
            if(d2 < Math::min(d3,d1)) t = DataType(0.5);
            if(d3 < Math::min(d1,d2)) t = DataType(0.8);
            /*
            DataType d1 = (point - p0).norm_euclid_sqr();
            DataType d2 = (point - map_on_segment(i, DataType(0.25))).norm_euclid_sqr();
            DataType d3 = (point - map_on_segment(i, DataType(0.50))).norm_euclid_sqr();
            DataType d4 = (point - map_on_segment(i, DataType(0.75))).norm_euclid_sqr();
            DataType d5 = (point - p3).norm_euclid_sqr();
            if(d1 < Math::min(Math::min(d2,d3), Math::min(d4,d5))) t = DataType(0.00);
            if(d2 < Math::min(Math::min(d3,d4), Math::min(d5,d1))) t = DataType(0.25);
            if(d3 < Math::min(Math::min(d4,d5), Math::min(d1,d2))) t = DataType(0.50);
            if(d4 < Math::min(Math::min(d5,d1), Math::min(d2,d3))) t = DataType(0.75);
            if(d5 < Math::min(Math::min(d1,d2), Math::min(d3,d4))) t = DataType(1.00);
            */
          }

          // Newton-Iteration
          for(int iter(0); iter < 10; ++iter)
          {
            // map local point and subtract world point: B(t) - X
            DataType vpx = (((-p0[0] + DataType(3)*p1[0] - DataType(3)*p2[0] + p3[0]) * t +
                            DataType(3) * (p0[0] - DataType(2)*p1[0] + p2[0])) * t +
                            DataType(3) * (-p0[0] + p1[0])) * t + p0[0] - point[0];
            DataType vpy = (((-p0[1] + DataType(3)*p1[1] - DataType(3)*p2[1] + p3[1]) * t +
                            DataType(3) * (p0[1] - DataType(2)*p1[1] + p2[1])) * t +
                            DataType(3) * (-p0[1] + p1[1])) * t + p0[1] - point[1];

            // compute first derivatives: B'(t)
            DataType d1x = DataType(3) * (((-p0[0] + DataType(3)*p1[0] - DataType(3)*p2[0] + p3[0]) * t +
                           DataType(2) * (p0[0] - DataType(2)*p1[0] + p2[0])) * t + (-p0[0] + p1[0]));
            DataType d1y = DataType(3) * (((-p0[1] + DataType(3)*p1[1] - DataType(3)*p2[1] + p3[1]) * t +
                           DataType(2) * (p0[1] - DataType(2)*p1[1] + p2[1])) * t + (-p0[1] + p1[1]));

            // compute second derivatives: B"(t)
            DataType d2x = DataType(6) * ((-p0[0] + DataType(3)*(p1[0] - p2[0]) + p3[0]) * t +
                          ( p0[0] - DataType(2)*p1[0] + p2[0]));
            DataType d2y = DataType(6) * ((-p0[1] + DataType(3)*(p1[1] - p2[1]) + p3[1]) * t +
                          ( p0[1] - DataType(2)*p1[1] + p2[1]));

            // compute function value: f(t) := < B'(t), B(t) - X >
            DataType fv = (d1x*vpx + d1y*vpy);
            if(Math::abs(fv) < DataType(1E-8))
              return t; // finished!

            // compute function derivative:  f'(t) := < B"(t), B(t) - X > + < B'(t), B'(t) >
            DataType fd = (d2x*vpx + d2y*vpy + d1x*d1x + d1y*d1y);
            if(Math::abs(fd) < DataType(1E-8))
            {
              // This should not happen...
              break;
            }

            // compute t-update:
            DataType tu = -fv / fd;

            // ensure that we do not try to run beyond our segment;
            // in this case the projected point is one of our segment ends
            if((t <= DataType(0)) && (tu <= DataType(0))) break;
            if((t >= DataType(1)) && (tu >= DataType(0))) break;

            // compute new t and clamp to [0,1]
            t = Math::clamp(t + tu, DataType(0), DataType(1));
          }

          // Maximum number of iterations reached...
          return t;
        }

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "bezier";
        }

        //virtual Mathvoid write_data_container(MeshStreamer::ChartContainer&) const override
        //{
        //  throw InternalError("Obsolete Bezier export not implemented");
        //}

        /** \copydoc ChartBase::write */
        virtual void write(std::ostream& os, const String& sindent) const override
        {
          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          os << sindent << "<Bezier dim=\"2\" size=\"" << this->_world.size() << "\"";
          os << " type=\"" << (this->_closed ? "closed" : "open") << "\">" << std::endl;
          this->write_points(os, sind);
          this->write_params(os, sind);
          this->write_controls(os, sind);
          os << sindent << "</Bezier>" << std::endl;
        }

      protected:
        void write_controls(std::ostream& os, const String& sindent) const
        {
          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          os << sindent << "<Controls>" << std::endl;
          for(std::size_t i(0); i < _control.size(); ++i)
          {
            os << sind << _control[i][0];
            for(int j(1); j < BaseClass::world_dim; ++j)
              os << " " << _control[i][j];
            ++i;
            for(int j(0); j < BaseClass::world_dim; ++j)
              os << " " << _control[i][j];
            os << std::endl;
          }
          os << sindent << "</Controls>" << std::endl;
        }
      }; // class Bezier<...>

      template<typename Mesh_>
      class BezierControlParser :
        public Xml::MarkupParser
      {
        typedef Bezier<Mesh_> ChartType;

      private:
        Bezier<Mesh_>& _spline;
        Index _size, _read;

      public:
        explicit BezierControlParser(Bezier<Mesh_>& spline, Index size) :
          _spline(spline), _size(size), _read(0)
        {
        }

        virtual bool attribs(std::map<String,bool>&) const override
        {
          return true;
        }

        virtual void create(int iline, const String& sline, const String&, const std::map<String, String>&, bool closed) override
        {
          if(closed)
            throw Xml::GrammarError(iline, sline, "Invalid closed markup");
        }

        virtual void close(int iline, const String& sline) override
        {
          // ensure that we have read all vertices
          if(_read < _size)
            throw Xml::GrammarError(iline, sline, "Invalid terminator; expected point");
        }

        virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
        {
          // no children allowed
          return nullptr;
        }

        virtual bool content(int iline, const String& sline) override
        {
          // make sure that we do not read more points than expected
          if(_read >= _size)
            throw Xml::ContentError(iline, sline, "Invalid content; exprected terminator");

          // split line by whitespaces
          std::deque<String> scoords;
          sline.split_by_charset(scoords);

          // check size
          if(scoords.size() != std::size_t(2*ChartType::world_dim))
            throw Xml::ContentError(iline, sline, "Invalid number of control point coordinates");

          typename ChartType::WorldPoint point1, point2;

          // try to parse all coords
          for(int i(0); i < ChartType::world_dim; ++i)
          {
            if(!scoords.at(std::size_t(i)).parse(point1[i]))
              throw Xml::ContentError(iline, sline, "Failed to parse point coordinate");
          }
          for(int i(0); i < ChartType::world_dim; ++i)
          {
            if(!scoords.at(std::size_t(ChartType::world_dim+i)).parse(point2[i]))
              throw Xml::ContentError(iline, sline, "Failed to parse point coordinate");
          }

          // push
          _spline.push_control(point1, point2);

          // okay, another point done
          ++_read;

          return true;
        }
      };

      template<typename Mesh_, typename ChartReturn_ = ChartBase<Mesh_>>
      class BezierChartParser :
        public Xml::MarkupParser
      {
      private:
        typedef Bezier<Mesh_> ChartType;
        typedef typename ChartType::DataType DataType;
        ChartReturn_*& _chart;
        ChartType* _spline;
        Index _size;

      public:
        explicit BezierChartParser(ChartReturn_*& chart) :
          _chart(chart),
          _spline(nullptr),
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
            throw Xml::GrammarError(iline, sline, "Invalid closed Bezier markup");

          Index dim(0);

          // try to parse the dimension
          if(!attrs.find("dim")->second.parse(dim))
            throw Xml::GrammarError(iline, sline, "Failed to parse Bezier dimension");
          if(dim != Index(2))
            throw Xml::GrammarError(iline, sline, "Invalid Bezier dimension");

          // try to parse the size
          if(!attrs.find("size")->second.parse(_size))
            throw Xml::GrammarError(iline, sline, "Failed to parse Bezier size");
          if(_size < Index(2))
            throw Xml::GrammarError(iline, sline, "Invalid Bezier size");

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
          _spline = new ChartType(poly_closed);
        }

        virtual void close(int, const String&) override
        {
          // okay
          _chart = _spline;
        }

        virtual bool content(int, const String&) override
        {
          return false;
        }

        virtual std::shared_ptr<Xml::MarkupParser> markup(int, const String&, const String& name) override
        {
          if(name == "Points")   return std::make_shared<SplinePointsParser<Bezier<Mesh_>>>(*_spline, _size);
          if(name == "Params")   return std::make_shared<SplineParamsParser<Bezier<Mesh_>>>(*_spline, _size);
          if(name == "Controls") return std::make_shared<BezierControlParser<Mesh_>>(*_spline, _size-1);
          return nullptr;
        }
      }; // class BezierChartParser<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAT
#endif // KERNEL_GEOMETRY_ATLAS_BEZIER_HPP
