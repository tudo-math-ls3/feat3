// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/runtime.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include "gendie_common.hpp"
#include "materials.hpp"



namespace Gendie
{
  using namespace FEAT;

  enum InflowGeometry
  {
    unkown = 0,
    circle = 1,
    ellipse = 2,
    ring = 3,
    rectangle = 4,
    curved_rectangle = 5
  };

  String stringify_enum(InflowGeometry type)
  {
    switch(type)
    {
      case unkown:
        return String("Unknown");
      case circle:
        return String("Circle");
      case ellipse:
        return String("Ellipse");
      case ring:
        return String("Ring");
      case rectangle:
        return String("Rectangle");
      case curved_rectangle:
        return String("Curved-Rectangle");
      default:
        XABORTM("InflowGeometry not known");
    }
    return String();
  }

  enum InflowType
  {
    best = 0,
    constant = 1,
    parabolic = 2,
    curved_flat = 3
  };

  String stringify_enum(InflowType type)
  {
    switch(type)
    {
      case best:
        return String("Best");
      case constant:
        return String("Constant");
      case parabolic:
        return String("Parabolic");
      case curved_flat:
        return String("Curved-Flat");
      default:
        XABORTM("InflowType not known");
    }
    return String();
  }

  template<typename MeshType_, typename DT_ = typename MeshType_::CoordType>
  class InflowBoundary;

  namespace Intern
  {
    template<typename DT_, int dim_>
    inline Tiny::Vector<DT_, dim_> project_vector(const Tiny::Vector<DT_, dim_>& point, const Tiny::Vector<DT_, dim_>& normal, DT_ dist)
    {
      ASSERTM(normal.normalized(), "Normal not normalized");
      const DT_ orth_dist = Tiny::dot(point, normal) - dist;
      return point - orth_dist * normal;
    }

    template<typename DT_, int dim_>
    inline Tiny::Vector<DT_, dim_> project_vector(const Tiny::Vector<DT_, dim_>& point, const Tiny::Vector<DT_, dim_>& normal, const Tiny::Vector<DT_, dim_>& lot_p)
    {
      ASSERTM(normal.normalized(), "Normal not normalized");
      const DT_ orth_dist = Tiny::dot(point - lot_p, normal);
      return point - orth_dist * normal;
    }

    bool parse_inflow_option(const std::pair<String, bool>& option, InflowType& in_type)
    {
      if((!option.second) || option.first.empty() || option.first.compare_no_case("best") || option.first.compare_no_case("default"))
      {
        in_type = InflowType::best;
        return true;
      }
      if((option.first.compare_no_case("constant") == 0))
      {
        in_type = InflowType::constant;
        return true;
      }
      else if (option.first.compare_no_case("parabolic") == 0)
      {
       in_type = InflowType::parabolic;
       return true;
      }
      else if(option.first.compare_no_case("curved_flat") == 0)
      {
        in_type = InflowType::curved_flat;
        return true;
      }
      XABORTM("No known inflow type");
      return false;
    }

    bool parse_inflow_geometry(const std::pair<String, bool>& option, InflowGeometry& in_type)
    {
      if(!option.second)
        return false;
      if(option.first.empty())
        XABORTM("No inflow geometry given");
      if(option.first.compare_no_case("circle") == 0)
      {
        in_type = InflowGeometry::circle;
        return true;
      }
      else if (option.first.compare_no_case("ellipse") == 0)
      {
       in_type = InflowGeometry::ellipse;
       return true;
      }
      else if (option.first.compare_no_case("ring") == 0)
      {
       in_type = InflowGeometry::ring;
       return true;
      }
      else if (option.first.compare_no_case("rectangle") == 0)
      {
       in_type = InflowGeometry::rectangle;
       return true;
      }
      else if (option.first.compare_no_case("curved_rectangle") == 0)
      {
       in_type = InflowGeometry::curved_rectangle;
       return true;
      }
      XABORTM("No known inflow type" + option.first);
      return false;
    }
  } // namespace Intern

  template<typename MeshType_, typename DT_>
  class InflowGeometryHandler
  {
  public:
    typedef DT_ DataType;
    static constexpr int dim = MeshType_::shape_dim;
    typedef Tiny::Vector<DataType, dim> VecType;
    typedef MeshType_ MeshType;
    typedef std::unique_ptr<FEAT::Geometry::MeshPart<MeshType>> MeshPartP;
    VecType center, flow_normal;//, bounds; //TODO: necessary?
    DataType d_plane;
    InflowGeometry inflow_geometry;
    InflowType inflow_type;
    int flow_id;
    // typedef std::function<bool(const VecType&)> HitTestFunc;
    typedef Analytic::SimplifiedLambdaVectorFunction3D<std::function<VecType(VecType)>> DirichletBoundaryFunc;
    // typedef std::function<VecType(const VecType&)> DirichletBoundaryFunc;

    virtual ~InflowGeometryHandler(){}
  private:
    std::array<std::array<VecType, 3>, 6> _generate_planes(const VecType& bb_min, const VecType& bb_max)
    {
      std::array<std::array<VecType, 3>, 6> tmp;
      // x-y min plane
      tmp[0] = {VecType{bb_min[0], bb_min[1], bb_min[2]}, VecType{DataType(1), DataType(0), DataType(0)}, VecType{DataType(0), DataType(1), DataType(0)}};
      // x-z min plane
      tmp[1] = {VecType{bb_min[0], bb_min[1], bb_min[2]}, VecType{DataType(1), DataType(0), DataType(0)}, VecType{DataType(0), DataType(0), DataType(1)}};
      // y-z min plane
      tmp[2] = {VecType{bb_min[0], bb_min[1], bb_min[2]}, VecType{DataType(0), DataType(1), DataType(0)}, VecType{DataType(0), DataType(0), DataType(1)}};
      // x-y max plane
      tmp[3] = {VecType{bb_max[0], bb_max[1], bb_max[2]}, VecType{-DataType(1), DataType(0), DataType(0)}, VecType{DataType(0), -DataType(1), DataType(0)}};
      // x-z max plane
      tmp[4] = {VecType{bb_max[0], bb_max[1], bb_max[2]}, VecType{-DataType(1), DataType(0), DataType(0)}, VecType{DataType(0), DataType(0), -DataType(1)}};
      // y-z max plane
      tmp[5] = {VecType{bb_max[0], bb_max[1], bb_max[2]}, VecType{DataType(0), -DataType(1), DataType(0)}, VecType{DataType(0), DataType(0), -DataType(1)}};

      return tmp;
    }

  protected:
    //formatting variables
    static constexpr int prec = 2;
    static constexpr int width = 5;

    VecType _parse_vector(const String& entry)
    {
      VecType tmp(DataType(0));
      std::deque<String> vals = entry.split_by_charset(",");
      XASSERTM(vals.size() == 3, "Entry " + entry + " does not contain 3 coordinates");
      for(int i = 0; i < 3; i++)
      {
        XASSERTM(vals[i].parse(tmp[i]), "Could not parse entry");
      }
      return tmp;
    }
  public:

    InflowGeometryHandler(const FEAT::PropertyMap* prop, const VecType& bb_min, const VecType& bb_max, int id, InflowType inflow_type_ = InflowType::best, InflowGeometry in_geo = InflowGeometry::unkown) : inflow_geometry(in_geo), inflow_type(inflow_type_), flow_id(id)
    {
      {
        auto [val, success] = prop->query("center");
        XASSERTM(success, "Could not parse center: ");
        center = _parse_vector(val);
      }
      // { //note face_normal can be safely ignored
      //   auto [val, success] = prop->query("face_normal");
      //   XASSERTM(success, "Could not parse face normal");
      //   face_normal = _parse_vector(val);
      // }
      {
        auto [val, success] = prop->query("flow_normal");
        XASSERTM(success, "Could not parse flow normal");
        flow_normal = _parse_vector(val);
        flow_normal.normalize();
      }
      {
        // by which factor we shift our center into normal direction
        DataType linear_factor = Math::Limits<DataType>::max();
        // project the center onto the boundingbox
        for(const std::array<VecType, 3>& vecs : _generate_planes(bb_min, bb_max))
        {
          const VecType& lot_p = vecs[0];
          const VecType& span1 = vecs[1];
          const VecType& span2 = vecs[2];
          //first of all, check if normal is linear independent to plane by calculating the cross product
          {
            VecType orto;
            Tiny::cross(orto, span1, span2);
            orto.normalize();
            if(Math::abs(Tiny::dot(orto, flow_normal)) <= Math::pow(Math::eps<DataType>(), DataType(0.8)))
              continue;
          }
          // now construct our linear system to be solved
          Tiny::Matrix<DataType, 3, 3> A{{-flow_normal[0], span1[0], span2[0]},
                                         {-flow_normal[1], span1[1], span2[1]},
                                         {-flow_normal[2], span1[2], span2[2]}};
          VecType rhs = center - lot_p;
          Tiny::Matrix<DataType, 3, 3> A_inv;
          A_inv.set_inverse(A);
          VecType sol = A_inv * rhs;

          // first entry is our new factor for projection, if it is absolutly smaller
          linear_factor = Math::abs(sol[0]) < Math::abs(linear_factor) ? sol[0] : linear_factor;

        }
        // print warning if the linear_factor is larger than 1%
        if(Math::abs(linear_factor) >= DataType(0.01))
          std::printf("WARNING: Projection factor %f larger than 0.01\n", double(linear_factor));
        // and now project
        center = center + linear_factor * flow_normal;
      }
      // precalculate the distance of the 'aufpunkt'
      d_plane = Tiny::dot(center, flow_normal);
    }

    virtual String get_bnd_name() const = 0;

    virtual MeshPartP create_mesh_part(const MeshType&, const DataType tolerance=1E-6) const = 0;

    virtual DirichletBoundaryFunc get_dirichlet_boundary_function(DataType throughput, const Material<DataType>& mat, const DataType) const = 0;

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const
    {
      String s;
      s += String("Inflow Geometry").pad_back(padlen, pc) + ": " + stringify_enum(inflow_geometry) + "\n";
      s += String("Inflow Type").pad_back(padlen, pc) + ": " + stringify_enum(inflow_type) + "\n";
      if constexpr(dim==3)
      {
        s += String("Inflow center").pad_back(padlen, pc) + ": [" + stringify_fp_fix(center[0], prec, width)
                      + ", " + stringify_fp_fix(center[1], prec, width) + ", " + stringify_fp_fix(center[2], prec, width) + "]\n";
        s += String("Inflow normal").pad_back(padlen, pc) + ": [" + stringify_fp_fix(flow_normal[0], prec, width)
                      + ", " + stringify_fp_fix(flow_normal[1], prec, width) + ", " + stringify_fp_fix(flow_normal[2], prec, width) + "]\n";
        // s += String("Inflow bounds").pad_back(padlen, pc) + ": [" + stringify_fp_fix(bounds[0], prec, width)
        //               + ", " + stringify_fp_fix(bounds[1], prec, width) + ", " + stringify_fp_fix(bounds[2], prec, width) + "]\n";
      }
      else
      {
        s += String("Inflow center").pad_back(padlen, pc) + ": [" + stringify_fp_fix(center[0], prec, width)
                      + ", " + stringify_fp_fix(center[1], prec, width) + "]\n";
        s += String("Inflow normal").pad_back(padlen, pc) + ": [" + stringify_fp_fix(flow_normal[0], prec, width)
                      + ", " + stringify_fp_fix(flow_normal[1], prec, width) + "]\n";
        // s += String("Inflow bounds").pad_back(padlen, pc) + ": [" + stringify_fp_fix(bounds[0], prec, width)
        //               + ", " + stringify_fp_fix(bounds[1], prec, width) + "]\n";
      }
      s += String("Inflow plane distance").pad_back(padlen, pc) + ": " + stringify_fp_fix(d_plane, prec, width);

      return s;
    }
  }; // class InflowGeometryHandler

  template<typename MeshType_, typename DT_ = typename MeshType_::DataType>
  class UnknownGeo : public InflowGeometryHandler<MeshType_, DT_>
  {
    public:
    typedef InflowGeometryHandler<MeshType_, DT_> BaseClass;
    typedef typename BaseClass::VecType VecType;
    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::MeshType MeshType;
    typedef typename BaseClass::MeshPartP MeshPartP;
    typedef typename BaseClass::DirichletBoundaryFunc DirichletBoundaryFunc;
    static constexpr int dim = BaseClass::dim;
    virtual ~UnknownGeo(){}
    explicit UnknownGeo(const FEAT::PropertyMap* prop, const VecType& bb_min, const VecType& bb_max, int id) : BaseClass(prop, bb_min, bb_max, id) {}
    auto get_hit_test_function(DT_) const
    {
      return [](const VecType&){return false;};
    }

    virtual String get_bnd_name() const override
    {
      return String("unknown_" + stringify(this->flow_id));
    }

    virtual MeshPartP create_mesh_part(const MeshType&, DT_) const override
    {
      XABORTM("Unknown geometry");
      return MeshPartP{};
    }

    virtual DirichletBoundaryFunc get_dirichlet_boundary_function(DataType, const Material<DataType>&, const DataType) const override
    {
      return DirichletBoundaryFunc([](auto a){return decltype(a)::null();});
    }

    virtual String format_string(std::size_t , const char) const override
    {
      return String("Unknown\n");
    }
  }; // class UnknownGeo

  template<typename MeshType_, typename DT_ = typename MeshType_::DataType>
  class Circle : public InflowGeometryHandler<MeshType_, DT_>
  {
    public:
    typedef InflowGeometryHandler<MeshType_, DT_> BaseClass;
    typedef typename BaseClass::VecType VecType;
    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::MeshType MeshType;
    typedef typename BaseClass::MeshPartP MeshPartP;
    typedef typename BaseClass::DirichletBoundaryFunc DirichletBoundaryFunc;
    static constexpr int dim = BaseClass::dim;
    DataType radius, inner_radius;


    virtual ~Circle(){}

    explicit Circle(const FEAT::PropertyMap* prop, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::curved_flat, DataType rel = DataType(0.8))
    : BaseClass(prop, bb_min, bb_max, id, in_type, InflowGeometry::circle),
      radius(DataType(0)),
      inner_radius(DataType(0))
    {
      PARSE_PROP_OPTION(prop, "radius", radius, false);
      inner_radius = radius * rel;
    }

    auto get_hit_test_function(DataType tolerance = DataType(1E-6)) const
    {
      return [tolerance, flow_n =this->flow_normal, cent=this->center, d_p = this->d_plane, rad=this->radius](const Tiny::Vector<DT_, 3>& in){
        // distance in orthognal direction
        const DataType orth_dist = Tiny::dot(in, flow_n) - d_p;
        const Tiny::Vector<DataType, 3> proj = in - orth_dist * flow_n;
        const Tiny::Vector<DataType, 3> dist = proj - cent;
        return (Math::abs(orth_dist) < tolerance) && (dist.norm_euclid_sqr() <= rad*rad);
      };
    }

    virtual String get_bnd_name() const override
    {
      return String("bnd:inflow_circle_") + stringify(this->flow_id);
    }

    virtual MeshPartP create_mesh_part(const MeshType& mesh, DT_ tolerance = DT_(1E-6)) const override
    {
      Geometry::HitTestFactory hit_test(get_hit_test_function(tolerance), mesh);
      return hit_test.make_unique();
    }

    virtual DirichletBoundaryFunc get_dirichlet_boundary_function(DataType throughput, const Material<DataType>& mat, const DataType mesh_unit = DataType(1E3)) const override
    {
      // always use SI as input -> vol_flow has to be in mm <- automate this!
      const DataType unit_trans = Math::cub(mesh_unit)/DataType(3600); // e.g. m^3/h -> mm^3/s
      const DataType vol_flow = throughput * unit_trans / mat.get_density();
      switch(this->inflow_type)
      {
        case InflowType::constant:
        {
          const DataType max_velo = vol_flow / (Math::pi<DataType>() * Math::sqr(radius));
          return DirichletBoundaryFunc([max_velo, _center=this->center, _radius=this->radius, _flow_normal=this->flow_normal, d_p=this->d_plane]
                  (const VecType& point) -> VecType
          {
            const VecType proj = Intern::project_vector(point, _flow_normal, d_p);
            return ((proj - _center).norm_euclid_sqr() >= Math::sqr(_radius) ? DataType(0) : max_velo) * _flow_normal ;
          });
        }
        case InflowType::parabolic:
        {
          const DataType max_velo = vol_flow / (DataType(0.5) * Math::pi<DataType>() * Math::sqr(radius) * Math::sqr(radius));
          return DirichletBoundaryFunc([max_velo, _center=this->center, _radius=this->radius, _flow_normal=this->flow_normal, d_p=this->d_plane]
                  (const VecType& point) -> VecType
          {
            const DataType dist = (Intern::project_vector(point, _flow_normal, d_p) - _center).norm_euclid();
            return (dist >= _radius) ? VecType::null() : (max_velo * (_radius - dist) * (_radius + dist)) * _flow_normal ;
          });
        }
        case InflowType::best:
        case InflowType::curved_flat:
        {
          const DataType k = vol_flow / (DataType(0.5) * Math::pi<DataType>() * (Math::sqr(radius) * Math::sqr(radius) - Math::sqr(inner_radius) * Math::sqr(inner_radius)));
          const DataType max_velo = (radius - inner_radius) * (radius + inner_radius) * k;
          // printf("Inner radius %f, max velo %f", inner_radius, max_velo);
          return DirichletBoundaryFunc([max_velo, k, _center=this->center, _radius=this->radius, _inner_radius=this->inner_radius, _flow_normal=this->flow_normal, d_p=this->d_plane]
                  (const VecType& point) -> VecType
          {
            const DataType dist = (Intern::project_vector(point, _flow_normal, d_p) - _center).norm_euclid();
            if(dist >= _radius)
              return VecType::null();
            return ((dist <= _inner_radius) ? max_velo : k * (_radius - dist) * (_radius + dist)) * _flow_normal ;
          });
        }
        default:
        {
          XABORTM("Not implemnted");
        }
      }
    }

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Inflow Radius").pad_back(padlen, pc) + ": " + stringify_fp_fix(radius, this->prec, this->width) + "\n";
      s += String("Inflow Inner-Radius").pad_back(padlen, pc) + ": " + stringify_fp_fix(inner_radius, this->prec, this->width);
      return s;
    }
  }; // class Circle

  template<typename MeshType_, typename DT_ = typename MeshType_::DataType>
  class Rectangle : public InflowGeometryHandler<MeshType_, DT_>
  {
    public:
    typedef InflowGeometryHandler<MeshType_, DT_> BaseClass;
    typedef typename BaseClass::VecType VecType;
    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::MeshType MeshType;
    typedef typename BaseClass::MeshPartP MeshPartP;
    typedef typename BaseClass::DirichletBoundaryFunc DirichletBoundaryFunc;
    static constexpr int dim = BaseClass::dim;
    /// lower left and upper right corners
    VecType a, b;


    virtual ~Rectangle(){}

    explicit Rectangle(const FEAT::PropertyMap* prop, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::parabolic)
    : BaseClass(prop, bb_min, bb_max, id, in_type, InflowGeometry::rectangle),
      a(),
      b()
    {
      {
        auto [val, success] = prop->query("midpoint_a"); //placeholder for boundary name
        XASSERTM(success, "Could not parse midpoint_a: ");
        a = this->_parse_vector(val);
      }
      {
        auto [val, success] = prop->query("midpoint_b"); //placeholder for boundary name
        XASSERTM(success, "Could not parse midpoint_b: ");
        b = this->_parse_vector(val);
      }
    }

    auto get_hit_test_function(DataType tolerance = DataType(1E-6)) const
    {
      return [tolerance, flow_n =this->flow_normal, cent=this->center, d_p = this->d_plane, _a=this->a, _b=this->b](const Tiny::Vector<DT_, 3>& in){
        // distance in orthognal direction
        const DataType orth_dist = Tiny::dot(in, flow_n) - d_p;
        const Tiny::Vector<DataType, 3> dC = in - cent;
        const Tiny::Vector<DataType, 3> aC = _a - cent;
        const Tiny::Vector<DataType, 3> bC = _b - cent;

        // projected distance: ||A - D||^2 - (A-D,P-D)^2 , resp. for B
        return (Math::abs(orth_dist) < tolerance) && (Math::sqr(aC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, aC)) > DataType(0))
                && (Math::sqr(bC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, bC)) > DataType(0));
      };
    }

    virtual String get_bnd_name() const override
    {
      return String("bnd:inflow_rectangle_") + stringify(this->flow_id);
    }

    virtual MeshPartP create_mesh_part(const MeshType& mesh, DT_ tolerance = DT_(1E-6)) const override
    {
      Geometry::HitTestFactory hit_test(get_hit_test_function(tolerance), mesh);
      return hit_test.make_unique();
    }

    virtual DirichletBoundaryFunc get_dirichlet_boundary_function(DataType throughput, const Material<DataType>& mat, const DataType mesh_unit = DataType(1E3)) const override
    {
      // always use SI as input -> vol_flow has to be in mm <- automate this!
      const DataType unit_trans = Math::cub(mesh_unit)/DataType(3600); // e.g. m^3/h -> mm^3/s
      const DataType vol_flow = throughput * unit_trans / mat.get_density();
      switch(this->inflow_type)
      {
        case InflowType::best:
        case InflowType::parabolic:
        {
          const VecType aC = a - this->center;
          const VecType bC = b - this->center;
          const DataType scale = vol_flow * DataType(9)/ (DataType(16) * Math::cub(aC.norm_euclid()) * Math::cub(bC.norm_euclid()));
          return DirichletBoundaryFunc([scale, _center=this->center, _flow_normal=this->flow_normal, d_p=this->d_plane, aC, bC]
                  (const VecType& point) -> VecType
          {
            const VecType proj = Intern::project_vector(point, _flow_normal, d_p);
            const VecType dC = proj - _center;
            const DataType aux1 = Math::sqr(aC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, aC));
            const DataType aux2 = Math::sqr(bC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, bC));
            const DataType daux = aux1 * aux2;
            return (((Math::sqr(aC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, aC)) > DataType(0))
                && (Math::sqr(bC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, bC)) > DataType(0))) ?
                  daux * scale : DataType(0)) * _flow_normal ;
          });
        }
        default:
        {
          XABORTM("Not implemented");
        }
      }
    }

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Inflow Midpoint A").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->a, this->prec, this->width) + "\n";
      s += String("Inflow Midpoint B").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->b, this->prec, this->width) + "\n";
      return s;
    }
  }; // class Rectangle

  template<typename MeshType_, typename DT_ = typename MeshType_::DataType>
  class CurvedRectangle : public Rectangle<MeshType_, DT_>
  {
    public:
    typedef Rectangle<MeshType_, DT_> BaseClass;
    typedef typename BaseClass::VecType VecType;
    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::MeshType MeshType;
    typedef typename BaseClass::MeshPartP MeshPartP;
    typedef typename BaseClass::DirichletBoundaryFunc DirichletBoundaryFunc;
    static constexpr int dim = BaseClass::dim;
    DataType radius;


    virtual ~CurvedRectangle(){}

    explicit CurvedRectangle(const FEAT::PropertyMap* prop, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::parabolic)
    : BaseClass(prop, bb_min, bb_max, id, in_type),
     radius((this->b - this->center).norm_euclid())
    {
      this->inflow_geometry = InflowGeometry::curved_rectangle;
    }

    auto get_hit_test_function(DataType tolerance = DataType(1E-6)) const
    {
      return [tolerance, flow_n =this->flow_normal, cent=this->center, d_p = this->d_plane, _a=this->a, _b=this->b, rad=this->radius](const Tiny::Vector<DT_, 3>& in){
        // distance in orthognal direction
        const DataType orth_dist = Tiny::dot(in, flow_n) - d_p;
        const Tiny::Vector<DataType, 3> dC = in - cent;
        const Tiny::Vector<DataType, 3> aC = _a - cent;
        const Tiny::Vector<DataType, 3> bC = _b - cent;

        const bool inside_rect = (Math::sqr(aC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, aC)) > DataType(0))
                && (Math::sqr(bC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, bC)) > DataType(0));
        const bool inside_outer_circles = (Math::sqr(rad) - (in-_a).norm_euclid_sqr() > DataType(0)) || (Math::sqr(rad) - (in-DataType(2)*cent+_a).norm_euclid_sqr() > DataType(0));

        return (Math::abs(orth_dist) < tolerance) && (inside_rect || inside_outer_circles) ;
      };
    }

    virtual String get_bnd_name() const override
    {
      return String("bnd:inflow_curved_rectangle_") + stringify(this->flow_id);
    }

    virtual MeshPartP create_mesh_part(const MeshType& mesh, DT_ tolerance = DT_(1E-6)) const override
    {
      Geometry::HitTestFactory hit_test(get_hit_test_function(tolerance), mesh);
      return hit_test.make_unique();
    }


    virtual DirichletBoundaryFunc get_dirichlet_boundary_function(DataType throughput, const Material<DataType>& mat, const DataType mesh_unit = DataType(1E3)) const override
    {
      // always use SI as input -> vol_flow has to be in mm <- automate this!
      const DataType unit_trans = Math::cub(mesh_unit)/DataType(3600); // e.g. m^3/h -> mm^3/s
      const DataType vol_flow = throughput * unit_trans / mat.get_density();
      switch(this->inflow_type)
      {
        case InflowType::best:
        case InflowType::parabolic:
        {
          const VecType aC = this->a - this->center;
          const VecType bC = this->b - this->center;
          DataType scale = (DataType(8) / DataType(3)) * aC.norm_euclid() * bC.norm_euclid()
                                  + DataType(2) * Math::pi<DataType>() * Math::sqr(radius);
          scale = vol_flow / scale;
          return DirichletBoundaryFunc([scale, _center=this->center, _flow_normal=this->flow_normal, d_p=this->d_plane, _a=this->a, radi=DataType(1)/this->radius, bC]
                  (const VecType& point) -> VecType
          {
            const VecType proj = Intern::project_vector(point, _flow_normal, d_p);
            const VecType dC = proj - _center;
            const DataType aux =  DataType(1) - Math::sqr(Tiny::dot(dC, bC))/Math::sqr(bC.norm_euclid_sqr());
            const DataType aux2 = Math::max(DataType(0), DataType(1) - (proj-_a).norm_euclid_sqr()*Math::sqr(radi))
                                  + Math::max(DataType(0), DataType(1) - (proj-DataType(2)*_center+_a).norm_euclid_sqr() * Math::sqr(radi));
            return (Math::max(aux, aux2) *scale) * _flow_normal;
          });
        }
        default:
        {
          XABORTM("Not implemented");
        }
      }
    }

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Inflow Radius").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->radius, this->prec, this->width) + "\n";
      return s;
    }
  }; // class CurvedRectangle

  template<typename MeshType_, typename DT_ = typename MeshType_::DataType>
  class Ellipse : public InflowGeometryHandler<MeshType_, DT_>
  {
    public:
    typedef InflowGeometryHandler<MeshType_, DT_> BaseClass;
    typedef typename BaseClass::VecType VecType;
    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::MeshType MeshType;
    typedef typename BaseClass::MeshPartP MeshPartP;
    typedef typename BaseClass::DirichletBoundaryFunc DirichletBoundaryFunc;
    static constexpr int dim = BaseClass::dim;
    DataType radius_major, radius_minor;


    virtual ~Ellipse(){}

    explicit Ellipse(const FEAT::PropertyMap* prop, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::parabolic)
    : BaseClass(prop, bb_min, bb_max, id, in_type, InflowGeometry::rectangle),
      radius_major(DataType(0)),
      radius_minor(DataType(0))
    {
      PARSE_PROP_OPTION(prop, "radius_major", radius_major, false);
      PARSE_PROP_OPTION(prop, "radius_minor", radius_minor, false);
    }

    auto get_hit_test_function(DataType tolerance = DataType(1E-6)) const
    {
      XABORTM("Ellipse has to be reimplemented");
      // return [tolerance, flow_n =this->flow_normal, cent=this->center, d_p = this->d_plane, a=this->radius_major, b=this->radius_minor](const Tiny::Vector<DT_, 3>& in){
      return [tolerance, flow_n =this->flow_normal, cent=this->center, d_p = this->d_plane](const Tiny::Vector<DT_, 3>& in){
        // distance in orthognal direction
        const DataType orth_dist = Tiny::dot(in, flow_n) - d_p;
        const Tiny::Vector<DataType, 3> dC = in - cent;

        // TODO: not enough information to describe ellipse
        return (Math::abs(orth_dist) < tolerance) && (Math::sqr(dC[0]));
      };
    }

    virtual String get_bnd_name() const override
    {
      return String("bnd:inflow_ellipse_") + stringify(this->flow_id);
    }

    virtual MeshPartP create_mesh_part(const MeshType& mesh, DT_ tolerance = DT_(1E-6)) const override
    {
      XABORTM("Ellipse not implemented yet");
      Geometry::HitTestFactory hit_test(get_hit_test_function(tolerance), mesh);
      return hit_test.make_unique();
    }

    virtual DirichletBoundaryFunc get_dirichlet_boundary_function(DataType throughput, const Material<DataType>& mat, const DataType mesh_unit = DataType(1E3)) const override
    {
      // always use SI as input -> vol_flow has to be in mm <- automate this!
      const DataType unit_trans = Math::cub(mesh_unit)/DataType(3600); // e.g. m^3/h -> mm^3/s
      const DataType vol_flow = throughput * unit_trans / mat.get_density();
      (void)vol_flow;
      switch(this->inflow_type)
      {
        case InflowType::best:
        case InflowType::parabolic:
        {
          // return DirichletBoundaryFunc([center=this->center, flow_normal=this->flow_normal, d_p=this->d_plane]
          return DirichletBoundaryFunc([_flow_normal=this->flow_normal, d_p=this->d_plane]
                  (const VecType& point) -> VecType
          {
            const VecType proj = Intern::project_vector(point, _flow_normal, d_p);
            (void)proj;
            return _flow_normal;
          });
        }
        default:
        {
          XABORTM("Not implemented");
        }
      }
    }

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Inflow Radius Major").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->radius_major, this->prec, this->width) + "\n";
      s += String("Inflow Radius Minor").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->radius_minor, this->prec, this->width) + "\n";
      return s;
    }
  }; // class Ellipse

  template<typename MeshType_, typename DT_ = typename MeshType_::DataType>
  class Ring : public InflowGeometryHandler<MeshType_, DT_>
  {
    public:
    typedef InflowGeometryHandler<MeshType_, DT_> BaseClass;
    typedef typename BaseClass::VecType VecType;
    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::MeshType MeshType;
    typedef typename BaseClass::MeshPartP MeshPartP;
    typedef typename BaseClass::DirichletBoundaryFunc DirichletBoundaryFunc;
    static constexpr int dim = BaseClass::dim;
    /// lower left and upper right corners
    DataType radius_outer;
    DataType radius_inner;


    virtual ~Ring(){}

    explicit Ring(const FEAT::PropertyMap* prop, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::parabolic)
    : BaseClass(prop, bb_min, bb_max, id, in_type, InflowGeometry::ring),
      radius_outer(DataType(0)),
      radius_inner(DataType(0))
    {
      PARSE_PROP_OPTION(prop, "outer_radius", radius_outer, false);
      PARSE_PROP_OPTION(prop, "inner_radius", radius_inner, false);
      XASSERTM(this->radius_outer > this->radius_inner, "Outer radius must be larger than inner one");
    }

    auto get_hit_test_function(DataType tolerance = DataType(1E-6)) const
    {
      return [tolerance, flow_n =this->flow_normal, cent=this->center, d_p = this->d_plane, rad_out=this->radius_outer, rad_in=this->radius_inner](const Tiny::Vector<DT_, 3>& in){
        // distance in orthognal direction
        const DataType orth_dist = Tiny::dot(in, flow_n) - d_p;
        const DataType dist = (in - cent).norm_euclid_sqr();

        return (Math::abs(orth_dist) < tolerance) && (Math::sqr(rad_out) > dist) && (Math::sqr(rad_in) < dist);
      };
    }

    virtual String get_bnd_name() const override
    {
      return String("bnd:inflow_ring_") + stringify(this->flow_id);
    }

    virtual MeshPartP create_mesh_part(const MeshType& mesh, DT_ tolerance = DT_(1E-6)) const override
    {
      Geometry::HitTestFactory hit_test(get_hit_test_function(tolerance), mesh);
      return hit_test.make_unique();
    }

    virtual DirichletBoundaryFunc get_dirichlet_boundary_function(DataType throughput, const Material<DataType>& mat, const DataType mesh_unit = DataType(1E3)) const override
    {
      // always use SI as input -> vol_flow has to be in mm <- automate this!
      const DataType unit_trans = Math::cub(mesh_unit)/DataType(3600); // e.g. m^3/h -> mm^3/s
      const DataType vol_flow = throughput * unit_trans / mat.get_density();
      switch(this->inflow_type)
      {
        case InflowType::best:
        case InflowType::parabolic:
        {
          const DataType integral = (Math::pi<DataType>() / DataType(6)) * (Math::cub(this->radius_outer - this->radius_inner))
                                                                            * (this->radius_inner + this->radius_outer);
          const DataType scale = vol_flow / integral;
          return DirichletBoundaryFunc([scale, _center=this->center, _flow_normal=this->flow_normal, d_p=this->d_plane, rad_out=this->radius_outer, rad_in=this->radius_inner]
                  (const VecType& point) -> VecType
          {
            const VecType proj = Intern::project_vector(point, _flow_normal, d_p);
            const DataType dist = (proj - _center).norm_euclid();
            return (((rad_out > dist) && (rad_in < dist)) ? ((rad_out-dist)*(dist-rad_in)*scale) : DataType(0)) * _flow_normal;
          });
        }
        default:
        {
          XABORTM("Not implemented");
        }
      }
    }

    virtual String format_string(std::size_t padlen = std::size_t(30), const char pc = '.') const override
    {
      String s = BaseClass::format_string(padlen, pc) + "\n";
      s += String("Inflow Outer Radius").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->radius_outer, this->prec, this->width) + "\n";
      s += String("Inflow Inner Radius").pad_back(padlen, pc) + ": " + stringify_fp_fix(this->radius_inner, this->prec, this->width) + "\n";
      return s;
    }
  }; // class Rectangle

  template<typename MeshType_, typename DT_ = typename MeshType_::DataType>
  std::unique_ptr<InflowGeometryHandler<MeshType_, DT_>> create_geo_handler(const PropertyMap* prop, const Tiny::Vector<DT_, 3>& bb_min, const Tiny::Vector<DT_, 3>& bb_max, int id, InflowType in_type)
  {
    InflowGeometry inflow_type = InflowGeometry::unkown;
    if(!Intern::parse_inflow_geometry(prop->get_entry("type"), inflow_type))
      XABORTM("Something went wrong while parsing inflow_type");
    switch(inflow_type)
    {
      case Gendie::InflowGeometry::circle:
      {
        return std::make_unique<Circle<MeshType_, DT_>>(prop, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::rectangle:
      {
        return std::make_unique<Rectangle<MeshType_, DT_>>(prop, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::curved_rectangle:
      {
        return std::make_unique<CurvedRectangle<MeshType_, DT_>>(prop, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::ring:
      {
        return std::make_unique<Ring<MeshType_, DT_>>(prop, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::ellipse:
      {
        XABORTM("ERROR: Ellipse definition not complete!");
        return std::make_unique<Ellipse<MeshType_, DT_>>(prop, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::unkown:
      {
        XABORTM("Unkown inflow geometry type");
        break;
      }
      default:
      {
        XABORTM("Not Implemented " + stringify(inflow_type) + " yet");
        break;
      }
    }
    return std::make_unique<UnknownGeo<MeshType_, DT_>>(prop, bb_min, bb_max, id);
  }


  template<typename Mesh_, typename DT_>
  class InflowBoundary
  {
  public:
    typedef Mesh_ MeshType;
    static constexpr int dim = MeshType::shape_dim;
    typedef Tiny::Vector<DT_, dim> PointType;
    typedef Tiny::Vector<DT_, dim> ImageType;
    typedef DT_ DataType;
    typedef InflowGeometryHandler<MeshType, DataType> GeoHandler;
    typedef typename GeoHandler::MeshPartP MeshPartP;
    typedef typename GeoHandler::DirichletBoundaryFunc DirichletBoundaryFunc;
  private:
    std::unique_ptr<GeoHandler> _geo_handler;
    Index _material_id;
    DataType _throughput;
    DataType _temperature;
    InflowGeometry _inflow_geometry;
    InflowType _inflow_type;
    static constexpr std::size_t padlen = std::size_t(30);

  public:
    /**
      * \brief Ctor for creating a inflow from a property-map referencing a vector of materials
      *
      * \param[in] inflow_property A property map consisting of
      *                              throughput=float
      *                              material_id=int
      *                              inflow_params/
      *                                      type=STRING (InflowGeometry)
      *                                      inflow_parameters...  (depending in type)
      *
      * \param[in] bb_min The minimal x,y,z coordinates of the boundingbox
      *
      * \param[in] bb_max The maximum x,y,z coordinates of the boundingbox
      *
      * \param[in] mesh Is ignored but used to infer the correct template parameter.
      *
      * \param[in] id The id of the meshpart -> used for naming schema
      *
      * \param[in] inflow_type Inflow type to be used, if default, the type is depending on the actual geometry
      */
    InflowBoundary(const FEAT::PropertyMap* inflow_property, const Tiny::Vector<DataType, 3>& bb_min, const Tiny::Vector<DataType, 3>& bb_max, int id, InflowType inflow_type = InflowType::best)
      : _material_id(~Index(0)),
        _throughput(DataType(0)),
        _temperature(DataType(300)),
        _inflow_type(inflow_type)
    {
      PARSE_PROP_OPTION(inflow_property, "throughput", _throughput, false);
      PARSE_PROP_OPTION(inflow_property, "material_id", _material_id, false);
      //TODO: for now ignore if not parsable
      PARSE_PROP_OPTION(inflow_property, "temperature", _temperature, true);

      //get sub property map
      const PropertyMap* sub_prop = inflow_property->get_sub_section("inflow_params");
      XASSERTM(sub_prop, "Could not find inflow params section");

      XASSERTM(Intern::parse_inflow_geometry(sub_prop->query("type"), _inflow_geometry), "Could not parse inflow option");
      _geo_handler = create_geo_handler<MeshType, DataType>(sub_prop, bb_min, bb_max, id, _inflow_type);
    }

    //provide other means of init

    /**
     * \brief Creates a meshpart defined on the given mesh
     *
     * \param[in] mesh The mesh the meshpart is to be created on
     *
     * \param[in] tolerance The thickness in normal direction in which points are registered as inside.
     *                      This should be greater zero.
     */
    MeshPartP create_mesh_part(const MeshType& mesh, DataType tolerance = DataType(1E-6)) const
    {
      ASSERTM(tolerance > DataType(0), "Tolerance not positive");
      return _geo_handler->create_mesh_part(mesh, tolerance);
    }

    /**
     * \brief Returns the boundary meshpart name
     *
     * \returns String with naming schema "bnd:inflow_type_id"
     */
    String get_bnd_name() const
    {
      return _geo_handler->get_bnd_name();
    }

    /**
     * \brief Create the inflow function used for the corresponding unit filter.
     *
     * \param[in] materials A vector of materials indexed by the material id
     */
    DirichletBoundaryFunc get_diri_inflow_function(const std::vector<Material<DataType>>& materials, const DataType mesh_unit_scale = DataType(1E3)) const
    {
      XASSERTM(materials.size() > _material_id, "Material vector and material id do not fit together");
      return _geo_handler->get_dirichlet_boundary_function(this->_throughput, materials.at(_material_id), mesh_unit_scale);
    }

    String format_string() const
    {
      String s;
      const char pc = '.';

      s = "Inflow " + get_bnd_name() + " parameters:\n";
      s += String("Material ID").pad_back(padlen, pc) + ": " + stringify(_material_id) + "\n";
      s += String("Throughput").pad_back(padlen, pc) + ": " + stringify_fp_fix(_throughput, 3, 10) + " kg/h\n";
      s += String("Temperature").pad_back(padlen, pc) + ": " + stringify_fp_fix(_temperature, 3, 10) + " C\n";
      s += _geo_handler->format_string();

      return s;
    }

    DataType get_temperature() const
    {
      return _temperature;
    }

    DataType get_throughput() const
    {
      return _throughput;
    }
  };

}
