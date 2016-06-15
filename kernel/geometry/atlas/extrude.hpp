#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_EXTRUDE_HPP
#define KERNEL_GEOMETRY_ATLAS_EXTRUDE_HPP 1

#include <kernel/geometry/atlas/chart.hpp>
// supported sub-charts:
#include <kernel/geometry/atlas/circle.hpp>
#include <kernel/geometry/atlas/bezier.hpp>
#include <kernel/geometry/atlas/polyline.hpp>

namespace FEAT
{
  namespace Geometry
  {
    namespace Atlas
    {
      template<typename SubChart_>
      struct ExtrudeTraits
      {
        /// we support explicit map
        static constexpr bool is_explicit = SubChart_::is_explicit;
        /// we support implicit projection
        static constexpr bool is_implicit = SubChart_::is_implicit;
        /// this is a 3D object
        static constexpr int world_dim = 3;
        /// we have 2D parameters
        static constexpr int param_dim = 2;
      };

      template<typename Mesh_, typename SubChart_>
      class Extrude :
        public ChartCRTP<Extrude<Mesh_, SubChart_>, Mesh_, ExtrudeTraits<SubChart_>>
      {
      public:
        /// The CRTP base class
        typedef ChartCRTP<Extrude<Mesh_, SubChart_>, Mesh_, ExtrudeTraits<SubChart_>> BaseClass;
        /// Floating point type
        typedef typename BaseClass::CoordType CoordType;
        /// Vector type for world points, aka image points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Vector type for parametrisation points, aka domain points
        typedef typename BaseClass::ParamPoint ParamPoint;

        typedef typename SubChart_::WorldPoint SubWorldPoint;
        typedef typename SubChart_::ParamPoint SubParamPoint;

      //protected:
        /// our sub-chart object
        SubChart_* _sub_chart;
        /// \todo trafo

      public:
        explicit Extrude() :
          _sub_chart(nullptr)
        {
        }

        virtual ~Extrude()
        {
          if(_sub_chart != nullptr)
            delete _sub_chart;
        }

        void set_sub_chart(SubChart_* sub_chart)
        {
          XASSERTM(_sub_chart == nullptr, "Extrude chart already has a sub-chart");
          _sub_chart = sub_chart;
        }

        bool can_implicit() const
        {
          return _sub_chart->can_implicit();
        }

        bool can_explicit() const
        {
          return _sub_chart->can_explicit();
        }

        void project_point(WorldPoint& point) const
        {
          // extract sub-world
          SubWorldPoint sub_point;
          sub_point[0] = point[0];
          sub_point[1] = point[1];
          _sub_chart->project_point(sub_point);
          point[0] = sub_point[0];
          point[1] = sub_point[1];
        }

        void project_meshpart(Mesh_& mesh, const MeshPart<Mesh_>& meshpart) const
        {
          auto& vtx = mesh.get_vertex_set();
          const auto& target_vtx = meshpart.template get_target_set<0>();
          for(Index i(0); i < meshpart.get_num_entities(0); ++i)
          {
            project_point(reinterpret_cast<WorldPoint&>(vtx[target_vtx[i]]));
          }
        }

        void map(WorldPoint& point, const ParamPoint& param) const
        {
          SubParamPoint sub_param;
          sub_param[0] = param[0];
          SubWorldPoint sub_point;
          _sub_chart->map(sub_point, sub_param);
          point[0] = sub_point[0];
          point[1] = sub_point[1];
          point[2] = param[1]; // Z-coord
        }

        /// \copydoc ChartBase::dist()
        CoordType compute_dist(const WorldPoint& point) const
        {
          SubWorldPoint sub_point;
          sub_point[0] = point[0];
          sub_point[1] = point[1];
          return _sub_chart->compute_dist(sub_point);
        }

        /// \copydoc ChartBase::dist()
        CoordType compute_dist(const WorldPoint& point, WorldPoint& grad_dist) const
        {
          SubWorldPoint sub_point;
          sub_point[0] = point[0];
          sub_point[1] = point[1];

          SubWorldPoint sub_grad_dist;

          CoordType my_dist(_sub_chart->compute_dist(sub_point, sub_grad_dist));

          ASSERT(my_dist >= CoordType(0));

          grad_dist[0] = sub_grad_dist[0];
          grad_dist[1] = sub_grad_dist[1];
          grad_dist[2] = CoordType(0);

          return my_dist;
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point) const
        {
          SubWorldPoint sub_point;
          sub_point[0] = point[0];
          sub_point[1] = point[1];
          return _sub_chart->compute_signed_dist(sub_point);
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point, WorldPoint& grad_dist) const
        {
          SubWorldPoint sub_point;
          sub_point[0] = point[0];
          sub_point[1] = point[1];

          SubWorldPoint sub_grad_dist;

          CoordType my_dist(_sub_chart->compute_dist(sub_point, sub_grad_dist));

          grad_dist[0] = sub_grad_dist[0];
          grad_dist[1] = sub_grad_dist[1];
          grad_dist[2] = CoordType(0);

          return my_dist;
        }

        virtual String get_type() const override
        {
          return "extrude";
        }

        virtual void write(std::ostream& os, const String& sindent) const override
        {
          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          os << sindent << "<Extrude>" << std::endl;
          _sub_chart->write(os, sind);
          os << sindent << "</Extrude>" << std::endl;
        }
      };

      template<typename Mesh_, bool enable_ = (Mesh_::shape_dim > 2)>
      class ExtrudeChartParser :
        public Xml::DummyParser
      {
      public:
        explicit ExtrudeChartParser(ChartBase<Mesh_>*&)
        {
          throw InternalError("Thou shall not arrive here");
        }
      };

      template<typename Mesh_>
      class ExtrudeChartParser<Mesh_, true> :
        public Xml::MarkupParser
      {
      private:
        ChartBase<Mesh_>*& _chart;

      public:
        explicit ExtrudeChartParser(ChartBase<Mesh_>*& chart) :
          _chart(chart)
        {
        }

        virtual bool attribs(std::map<String,bool>& attrs) const override
        {
          // allowed for upwards compatibility
          attrs.emplace("origin", false);
          attrs.emplace("offset", false);
          attrs.emplace("transform", false);
          return true;
        }

        virtual void create(
          int,
          const String& ,
          const String&,
          const std::map<String, String>&,
          bool) override
        {
        }

        virtual void close(int iline, const String& sline) override
        {
          if(_chart == nullptr)
            throw Xml::GrammarError(iline, sline, "invalid empty <Extrude> markup");
        }

        virtual bool content(int, const String&) override
        {
          return false;
        }

        virtual std::shared_ptr<Xml::MarkupParser> markup(int, const String&, const String& name) override
        {
          typedef typename Mesh_::ShapeType ShapeType;
          typedef typename Mesh_::CoordType CoordType;
          typedef typename Shape::FaceTraits<ShapeType, 1>::ShapeType SubShapeType;
          typedef ConformalMesh<SubShapeType, 2, 2, CoordType> SubMeshType;

          // What have we here?
          if(name == "Bezier")
          {
            auto* ext = new Extrude<Mesh_, Bezier<SubMeshType>>();
            _chart = ext;
            return std::make_shared<Atlas::BezierChartParser<SubMeshType, Bezier<SubMeshType>>>(ext->_sub_chart);
          }
          if(name == "Circle")
          {
            auto* ext = new Extrude<Mesh_, Circle<SubMeshType>>();
            _chart = ext;
            return std::make_shared<Atlas::CircleChartParser<SubMeshType, Circle<SubMeshType>>>(ext->_sub_chart);
          }
          if(name == "Polyline")
          {
            auto* ext = new Extrude<Mesh_, Polyline<SubMeshType>>();
            _chart = ext;
            return std::make_shared<Atlas::PolylineChartParser<SubMeshType, Polyline<SubMeshType>>>(ext->_sub_chart);
          }

          return nullptr;
        }
      }; // class ExtrudeChartParser<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_ATLAS_EXTRUDE_HPP
