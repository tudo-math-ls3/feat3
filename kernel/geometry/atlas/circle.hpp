#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_CIRCLE_HPP
#define KERNEL_GEOMETRY_ATLAS_CIRCLE_HPP 1

#include <kernel/geometry/atlas/chart.hpp>

namespace FEAT
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
        /// The CRTP base class
        typedef ChartCRTP<Circle<Mesh_>, Mesh_, CircleTraits> BaseClass;
        /// Floating point type
        typedef typename BaseClass::CoordType CoordType;
        /// Vector type for world points, aka image points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Vector type for parametrisation points, aka domain points
        typedef typename BaseClass::ParamPoint ParamPoint;

      protected:
        /// the circle's midpoint
        WorldPoint _midpoint;
        /// the circle's radius
        CoordType _radius;
        /// Specifies whether the circle mapping has a domain
        bool _have_domain;
        /// Left parametrisation domain boundary
        CoordType _trafo_a;
        /// Right parametrisation domain boundary
        CoordType _trafo_b;

      public:
        /**
         * \brief Creates a Circle chart for implicit adaption
         *
         * \param[in] mid_x, mid_y
         * The coordinates of the circle midpoint.
         *
         * \param[in] radius
         * The radius of the circle. Must be positive.
         */
        explicit Circle(CoordType mid_x, CoordType mid_y, CoordType radius) :
          _radius(radius),
          _have_domain(false),
          _trafo_a(0.0),
          _trafo_b(0.0)
        {
          XASSERTM(radius > CoordType(0), "invalid circle radius");
          _midpoint[0] = mid_x;
          _midpoint[1] = mid_y;
        }

        /**
         * \brief Creates a Circle chart for both implicit and explicit adaption
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
        explicit Circle(CoordType mid_x, CoordType mid_y, CoordType radius,
          CoordType param_l, CoordType param_r) :
          _radius(radius),
          _have_domain(true),
          _trafo_a(-param_l),
          _trafo_b((CoordType(2) * Math::pi<CoordType>()) / (param_r - param_l))
        {
          XASSERTM(radius > CoordType(0), "invalid circle radius");
          _midpoint[0] = mid_x;
          _midpoint[1] = mid_y;
        }

        bool can_explicit() const
        {
          return _have_domain;
        }

        /**
         * \brief Projects a single world point
         *
         * \param[in,out] point
         * The world point to be projected
         *
         */
        void project_point(WorldPoint& point) const
        {
          WorldPoint dist_vec(point - _midpoint);
          CoordType distance(dist_vec.norm_euclid());

          if(distance < Math::eps<CoordType>())
          {
            dist_vec(0) = _radius;
            point += dist_vec;
          }
          else
          {
            point = _midpoint + (_radius / distance)*dist_vec;
          }
        }

        /// \copydoc ChartBase::dist()
        CoordType compute_dist(const WorldPoint& point) const
        {
          return Math::abs(compute_signed_dist(point));
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point) const
        {
          return (point - _midpoint).norm_euclid() - _radius;
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
            project_point(reinterpret_cast<WorldPoint&>(vtx[target_vtx[i]]));
          }
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
          // transform parameter to interval [0, 2*pi)
          CoordType x = (param[0] + _trafo_a) * _trafo_b;
          point[0] = this->_midpoint[0] + this->_radius * Math::cos(x);
          point[1] = this->_midpoint[1] + this->_radius * Math::sin(x);
        }

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "circle";
        }

        /** \copydoc ChartBase::write */
        virtual void write(std::ostream& os, const String& sindent) const override
        {
          os << sindent << "<Circle";
          os << " radius=\"" << this->_radius << "\"";
          os << " midpoint=\"" << this->_midpoint[0] << " " << this->_midpoint[1] << "\"";
          // reconstruct domain from trafo
          if(_have_domain)
          {
            CoordType param_l(-_trafo_a);
            CoordType param_r(param_l + CoordType(2) * Math::pi<CoordType>() / _trafo_b);
            os << " domain=\"" << param_l << " " << param_r << "\"";
          }
          os << " />" << std::endl;
        }
      }; // class Circle

      template<typename Mesh_, typename ChartReturn_ = ChartBase<Mesh_>>
      class CircleChartParser :
        public Xml::MarkupParser
      {
      private:
        typedef Circle<Mesh_> ChartType;
        typedef typename ChartType::CoordType CoordType;
        ChartReturn_*& _chart;

      public:
        explicit CircleChartParser(ChartReturn_*& chart) :
          _chart(chart)
        {
        }

        virtual bool attribs(std::map<String,bool>& attrs) const override
        {
          attrs.emplace("radius", true);
          attrs.emplace("midpoint", true);
          attrs.emplace("domain", false);
          return true;
        }

        virtual void create(
          int iline,
          const String& sline,
          const String&,
          const std::map<String, String>& attrs,
          bool) override
        {
          CoordType radius = CoordType(0);
          CoordType mid_x = CoordType(0);
          CoordType mid_y = CoordType(0);
          CoordType dom_0 = CoordType(0);
          CoordType dom_1 = CoordType(1);
          bool have_domain(false);

          // try to parse the radius
          if(!attrs.find("radius")->second.parse(radius))
            throw Xml::GrammarError(iline, sline, "Failed to parse circle radius");
          if(radius < CoordType(1E-5))
            throw Xml::GrammarError(iline, sline, "Invalid circle radius");

          // try to parse midpoind
          std::deque<String> mids;
          attrs.find("midpoint")->second.split_by_charset(mids);
          if(mids.size() != std::size_t(2))
            throw Xml::GrammarError(iline, sline, "Invalid circle midpoint string");
          if(!mids.front().parse(mid_x) || !mids.back().parse(mid_y))
            throw Xml::GrammarError(iline, sline, "'Failed to parse circle midpoint");

          // do we have a domain?
          auto it = attrs.find("domain");
          if(it != attrs.end())
          {
            have_domain = true;
            std::deque<String> doms;
            it->second.split_by_charset(doms);
            if(doms.size() != std::size_t(2))
              throw Xml::GrammarError(iline, sline, "Invalid circle domain string");
            if(!doms.front().parse(dom_0) || !doms.back().parse(dom_1))
              throw Xml::GrammarError(iline, sline, "'Failed to parse circle domain");
          }

          // everything seems fine, let's create the chart then
          if(have_domain)
            _chart = new ChartType(mid_x, mid_y, radius, dom_0, dom_1);
          else
            _chart = new ChartType(mid_x, mid_y, radius);
        }

        virtual void close(int, const String&) override
        {
        }

        virtual bool content(int, const String&) override
        {
          return false;
        }

        virtual std::shared_ptr<Xml::MarkupParser> markup(int, const String&, const String&) override
        {
          return nullptr;
        }
      }; // class CircleChartParser<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAT
#endif // KERNEL_GEOMETRY_ATLAS_CIRCLE_HPP
