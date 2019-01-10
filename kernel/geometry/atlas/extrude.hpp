#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_EXTRUDE_HPP
#define KERNEL_GEOMETRY_ATLAS_EXTRUDE_HPP 1

#include <kernel/geometry/atlas/chart.hpp>
// supported sub-charts:
#include <kernel/geometry/atlas/bezier.hpp>
#include <kernel/geometry/atlas/circle.hpp>

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

      protected:
        /// the origin point
        Tiny::Vector<CoordType, 2> _origin;
        /// the offset point
        Tiny::Vector<CoordType, 3> _offset;
        /// the rotation matrix
        Tiny::Matrix<CoordType, 3, 3> _rotation;

      public:
        explicit Extrude() :
          _sub_chart(nullptr)
        {
          _origin.format();
          _offset.format();
          _rotation.set_identity();
        }

        explicit Extrude(SubChart_* sub_chart) :
          _sub_chart(sub_chart)
        {
          _origin.format();
          _offset.format();
          _rotation.set_identity();
        }

        virtual ~Extrude()
        {
          if(_sub_chart != nullptr)
            delete _sub_chart;
        }

        void set_origin(CoordType x, CoordType y)
        {
          _origin[0] = x;
          _origin[1] = y;
        }

        void set_offset(CoordType x, CoordType y, CoordType z)
        {
          _offset[0] = x;
          _offset[1] = y;
          _offset[2] = z;
        }

        void set_angles(CoordType yaw, CoordType pitch, CoordType roll)
        {
          _rotation.set_rotation_3d(yaw, pitch, roll);
        }

        void set_sub_chart(SubChart_* sub_chart)
        {
          XASSERTM(_sub_chart == nullptr, "Extrude chart already has a sub-chart");
          _sub_chart = sub_chart;
        }

        virtual bool can_implicit() const override
        {
          return _sub_chart->can_implicit();
        }

        virtual bool can_explicit() const override
        {
          return _sub_chart->can_explicit();
        }

        void transform_2d_to_3d(WorldPoint& p, const SubWorldPoint& q, const CoordType z = CoordType(0)) const
        {
          // apply our extrude transformation:
          // p = w + R * (q - v)
          WorldPoint x;
          x[0] = q[0] - _origin[0];
          x[1] = q[1] - _origin[1];
          x[2] = z;
          p.set_mat_vec_mult(_rotation, x);
          p += _offset;
        }

        void transform_3d_to_2d(const WorldPoint& p, SubWorldPoint& q, CoordType& z) const
        {
          // apply our inverse extrude transformation:
          // q = v + R^{-1} * (p - w)
          WorldPoint x = p - _offset;
          WorldPoint y;
          // rotation matrices are orthogonal, i.e. R^{-1} = R^T
          y.set_vec_mat_mult(x, _rotation);
          q[0] = y[0] - _origin[0];
          q[1] = y[1] - _origin[1];
          z    = y[2];
        }

        void transform_3d_to_2d(const WorldPoint& p, SubWorldPoint& q) const
        {
          CoordType z = CoordType(0);
          transform_3d_to_2d(p, q, z);
        }

        virtual void transform(const WorldPoint& origin, const WorldPoint& angles, const WorldPoint& offset) override
        {
          // the rotation matrix
          Tiny::Matrix<CoordType, 3, 3> rotation;
          rotation.set_rotation_3d(angles[0], angles[1], angles[2]);

          // we now have 2 transformations:
          // 1. from the extrude chart itself and
          // 2. from the parameters given to this function
          //
          // The new extrusion chart transformation is given by the concatenation of the
          // two above transformations:
          //
          //           v2 + R2 * ( (v1 + R1*(x - w1)) - w2)
          //         = v2 + R2*(v1 - w2) + R2*R1*(x - w1)
          //           \---------------/   \---/      \/
          //                  v3             R3       w3
          //
          // where
          //   vi = offset
          //   Ri = rotation
          //   wi = origin

          WorldPoint tmp1, tmp2;

          // let us first update the offset vector v3:
          // compute offset vector v3 := v2 + R1*(v1 - w2)
          tmp1 = _offset - origin;
          tmp2.set_mat_vec_mult(rotation, tmp1);
          _offset = offset + tmp2;

          // compute rotation matrix R3 := R2*R1
          Tiny::Matrix<CoordType, 3, 3> old_rot(_rotation);
          _rotation.set_mat_mat_mult(rotation, old_rot);

          // our origin vector w3 remains the same (= w1)
        }

        void project_point(WorldPoint& point) const
        {
          // transform to 2d
          SubWorldPoint sub_point;
          CoordType z = CoordType(0);
          transform_3d_to_2d(point, sub_point, z);
          // project in 2d
          _sub_chart->project_point(sub_point);
          // transform to 3d
          transform_2d_to_3d(point, sub_point, z);
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

        void map_param(WorldPoint& point, const ParamPoint& param) const
        {
          // map the parameter in 2D
          SubParamPoint sub_param;
          sub_param[0] = param[0];
          SubWorldPoint sub_point;
          _sub_chart->map_param(sub_point, sub_param);
          // transform to 3d
          transform_2d_to_3d(point, sub_point, param[1]);
        }

        /// \copydoc ChartBase::dist()
        CoordType compute_dist(const WorldPoint& point) const
        {
          // transform to 2d
          SubWorldPoint sub_point;
          transform_3d_to_2d(point, sub_point);
          return _sub_chart->compute_dist(sub_point);
        }

        /// \copydoc ChartBase::dist()
        CoordType compute_dist(const WorldPoint& point, WorldPoint& grad_dist) const
        {
          // transform to 2d
          SubWorldPoint sub_point;
          CoordType z = CoordType(0);
          transform_3d_to_2d(point, sub_point, z);

          SubWorldPoint sub_grad_dist;

          CoordType my_dist(_sub_chart->compute_dist(sub_point, sub_grad_dist));

          ASSERT(my_dist >= CoordType(0));

          // transform to 3d
          transform_2d_to_3d(grad_dist, sub_grad_dist, z);

          return my_dist;
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point) const
        {
          // transform to 2d
          SubWorldPoint sub_point;
          transform_3d_to_2d(point, sub_point);
          return _sub_chart->compute_signed_dist(sub_point);
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point, WorldPoint& grad_dist) const
        {
          // transform to 2d
          SubWorldPoint sub_point;
          CoordType z = CoordType(0);
          transform_3d_to_2d(point, sub_point, z);

          SubWorldPoint sub_grad_dist;

          CoordType my_dist(_sub_chart->compute_dist(sub_point, sub_grad_dist));

          // transform to 3d
          transform_2d_to_3d(grad_dist, sub_grad_dist, z);

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

          const CoordType tol = Math::pow(Math::eps<CoordType>(), CoordType(0.7));

          os << sindent << "<Extrude";
          if(_origin.norm_euclid_sqr() > tol)
            os << " origin=\"" << _origin[0] << " " << _origin[1] << "\"";
          if(_offset.norm_euclid_sqr() > tol)
            os << " offset=\"" << _offset[0] << " " << _offset[1] << " " << _offset[2] << "\"";

          // check whether our rotation matrix is the identitiy matrix
          if(Math::sqr(_rotation.norm_sub_id_frobenius()) > tol)
          {
            // Now that's tricky: we have to reconstruct the yaw-pitch-roll angles
            // from the rotation matrix:  see \cite Craig04, page 43, eq. (2.66) -- (2.68)
            const CoordType pi = Math::pi<CoordType>();
            const CoordType pi_2 = CoordType(0.5) * pi;
            const CoordType r11 = _rotation(0,0);
            const CoordType r12 = _rotation(0,1);
            const CoordType r21 = _rotation(1,0);
            const CoordType r22 = _rotation(1,1);
            const CoordType r31 = _rotation(2,0);
            const CoordType r32 = _rotation(2,1);
            const CoordType r33 = _rotation(2,2);
            CoordType ay = CoordType(0); // yaw   (in book: alpha)
            CoordType ap = CoordType(0); // pitch (in book: beta)
            CoordType ar = CoordType(0); // roll  (in book: gamma)

            // determine the pitch first
            ap = Math::atan2(-r31, Math::sqrt(r11*r11 + r21*r21));

            // check pitch for "gimbal lock" angles +pi/2 and -pi/2
            if(ap >= pi_2 - tol)
            {
              ap = pi_2;
              ar = Math::atan2(r12, r22);
            }
            else if(ap <= -pi_2 + tol)
            {
              ap = -pi_2;
              ar = -Math::atan2(r12, r22);
            }
            else // no gimbal lock
            {
              // Note: we can safely drop the "c\beta" denominators from the book,
              // as these always cancel out by definition of atan2
              ay = Math::atan2(r21, r11);
              ar = Math::atan2(r32, r33);
            }

            // convert angles from radians to revolutions and write out
            const CoordType mult = CoordType(0.5) / pi;
            os << " angles=\"" << (mult*ay) << " " << (mult*ap) << " " << (mult*ar) << "\"";
          }
          os << ">" << std::endl;

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
          throw InternalError(__func__,__FILE__,__LINE__,"Thou shall not arrive here");
        }
      };

      template<typename Mesh_>
      class ExtrudeChartParser<Mesh_, true> :
        public Xml::MarkupParser
      {
      private:
        typedef typename ChartBase<Mesh_>::CoordType CoordType;
        ChartBase<Mesh_>*& _chart;
        CoordType _ori_x, _ori_y;
        CoordType _off_x, _off_y, _off_z;
        CoordType _yaw, _pitch, _roll;

      public:
        explicit ExtrudeChartParser(ChartBase<Mesh_>*& chart) :
          _chart(chart),
          _ori_x(CoordType(0)),
          _ori_y(CoordType(0)),
          _off_x(CoordType(0)),
          _off_y(CoordType(0)),
          _off_z(CoordType(0)),
          _yaw(CoordType(0)),
          _pitch(CoordType(0)),
          _roll(CoordType(0))
        {
        }

        virtual bool attribs(std::map<String,bool>& attrs) const override
        {
          // allowed for upwards compatibility
          attrs.emplace("origin", false);
          attrs.emplace("offset", false);
          attrs.emplace("angles", false);
          return true;
        }

        virtual void create(
          int iline,
          const String& sline,
          const String&,
          const std::map<String, String>& attrs,
          bool) override
        {
          // parse 2D origin (if given)
          {
            auto it = attrs.find("origin");
            if(it != attrs.end())
            {
              std::deque<String> sori = it->second.split_by_whitespaces();
              if(sori.size() != std::size_t(2))
                throw Xml::GrammarError(iline, sline, "Invalid Extrude chart origin attribute");
              if(!sori.front().parse(_ori_x) || !sori.back().parse(_ori_y))
                throw Xml::GrammarError(iline, sline, "'Failed to parse extrude origin attribute");

            }
          }

          // parse 3D offset (if given)
          {
            auto it = attrs.find("offset");
            if(it != attrs.end())
            {
              std::deque<String> soff = it->second.split_by_whitespaces();
              if(soff.size() != std::size_t(3))
                throw Xml::GrammarError(iline, sline, "Invalid Extrude chart offset attribute");
              if(!soff.at(0).parse(_off_x) || !soff.at(1).parse(_off_y) || !soff.at(2).parse(_off_z))
                throw Xml::GrammarError(iline, sline, "'Failed to parse extrude offset attribute");

            }
          }

          // parse angles (if given)
          {
            auto it = attrs.find("angles");
            if(it != attrs.end())
            {
              std::deque<String> sang = it->second.split_by_whitespaces();
              if(sang.size() != std::size_t(3))
                throw Xml::GrammarError(iline, sline, "Invalid Extrude chart angles attribute");
              if(!sang.at(0).parse(_yaw) || !sang.at(1).parse(_pitch) || !sang.at(2).parse(_roll))
                throw Xml::GrammarError(iline, sline, "'Failed to parse extrude angles attribute");

              // Note: the angles are specifies in revolutions, but we need radians, so multiply by 2*pi
              const CoordType mult = CoordType(2) * Math::pi<CoordType>();
              _yaw *= mult;
              _pitch *= mult;
              _roll *= mult;
            }
          }
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
          typedef typename Shape::FaceTraits<ShapeType, 1>::ShapeType SubShapeType;
          typedef ConformalMesh<SubShapeType, 2, CoordType> SubMeshType;

          // What have we here?
          if(name == "Circle")
          {
            auto* ext = new Extrude<Mesh_, Circle<SubMeshType>>();
            ext->set_origin(_ori_x, _ori_y);
            ext->set_offset(_off_x, _off_y, _off_z);
            ext->set_angles(_yaw, _pitch, _roll);
            _chart = ext;
            return std::make_shared<Atlas::CircleChartParser<SubMeshType, Circle<SubMeshType>>>(ext->_sub_chart);
          }
          if(name == "Bezier")
          {
            auto* ext = new Extrude<Mesh_, Bezier<SubMeshType>>();
            ext->set_origin(_ori_x, _ori_y);
            ext->set_offset(_off_x, _off_y, _off_z);
            ext->set_angles(_yaw, _pitch, _roll);
            _chart = ext;
            return std::make_shared<Atlas::BezierChartParser<SubMeshType, Bezier<SubMeshType>>>(ext->_sub_chart);
          }

          return nullptr;
        }
      }; // class ExtrudeChartParser<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_ATLAS_EXTRUDE_HPP
