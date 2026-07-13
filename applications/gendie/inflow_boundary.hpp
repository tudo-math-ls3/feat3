// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/runtime.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/util/arg_parser.hpp>
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

  enum class InflowGeometry : std::uint8_t
  {
    unkown = 0,
    circle,
    ellipse,
    ring,
    rectangle,
    curved_rectangle
  };

  std::ostream& operator<<(std::ostream& os, InflowGeometry geo)
  {
    switch(geo)
    {
      case InflowGeometry::unkown:
        return os << "Unknown";
      case InflowGeometry::circle:
        return os << "Circle";
      case InflowGeometry::ellipse:
        return os << "Ellipse";
      case InflowGeometry::ring:
        return os << "Ring";
      case InflowGeometry::rectangle:
        return os << "Rectangle";
      case InflowGeometry::curved_rectangle:
        return os << "Curved-Rectangle";
      default:
        XABORTM("InflowGeometry not known");
    }
    return os;
  }

  enum class InflowType : std::uint8_t
  {
    best = 0,
    constant,
    parabolic,
    curved_flat
  };

  std::ostream& operator<<(std::ostream& os, InflowType type)
  {
    switch(type)
    {
      case InflowType::best:
        return os << "Best";
      case InflowType::constant:
        return os << "Constant";
      case InflowType::parabolic:
        return os << "Parabolic";
      case InflowType::curved_flat:
        return os << "Curved-Flat";
      default:
        XABORTM("InflowType not known");
    }
    return os;
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

    bool try_parse_inflow_option(const String& option, InflowType& in_type)
    {
      if(option.empty() || (option.compare_no_case("best") == 0) || (option.compare_no_case("default") == 0))
      {
        in_type = InflowType::best;
        return true;
      }
      if((option.compare_no_case("constant") == 0))
      {
        in_type = InflowType::constant;
        return true;
      }
      else if (option.compare_no_case("parabolic") == 0)
      {
       in_type = InflowType::parabolic;
       return true;
      }
      else if(option.compare_no_case("curved_flat") == 0)
      {
        in_type = InflowType::curved_flat;
        return true;
      }
      return false;
    }

    bool try_parse_inflow_geometry(const String& option, InflowGeometry& in_type)
    {
      if(option.compare_no_case("circle") == 0)
      {
        in_type = InflowGeometry::circle;
        return true;
      }
      else if (option.compare_no_case("ellipse") == 0)
      {
       in_type = InflowGeometry::ellipse;
       return true;
      }
      else if (option.compare_no_case("ring") == 0)
      {
       in_type = InflowGeometry::ring;
       return true;
      }
      else if (option.compare_no_case("rectangle") == 0)
      {
       in_type = InflowGeometry::rectangle;
       return true;
      }
      else if (option.compare_no_case("curved_rectangle") == 0)
      {
       in_type = InflowGeometry::curved_rectangle;
       return true;
      }
      return false;
    }

    static bool is_inflow_profile(const String& s)
    {
      InflowType type = InflowType::best;
      return try_parse_inflow_option(s, type);
    }

    static bool is_inflow_geometry(const String& s)
    {
      InflowGeometry type = InflowGeometry::unkown;
      return try_parse_inflow_geometry(s, type);
    }

    static bool is_circle(const String& s)
    {
      return s.compare_no_case("circle") == 0;
    }

    static bool is_rectangle(const String& s)
    {
      return s.compare_no_case("rectangle") == 0 || s.compare_no_case("curved_rectangle") == 0;
    }

    static bool is_ring(const String& s)
    {
      return s.compare_no_case("ring") == 0;
    }

    static bool is_ellipse(const String& s)
    {
      return s.compare_no_case("ellipse") == 0;
    }

    template<typename DT_>
    struct InflowGeometryParser : public ArgParser
    {
      Parameter<DT_> throughput = parameter(DT_(0.0)).required().property("throughput");
      Parameter<DT_> temperature = parameter(DT_(300.0)).property("temperature");
      Parameter<Index> material_id = parameter(Index(0)).required().property("material_id");
      Parameter<String> inflow_profile = parameter(String("BEST"))
                                              .property("inflow_params/profile")
                                              .name("Inflow profile type")
                                              .help_string("Valid options are: BEST, DEFAULT, CONSTANT, PARABOLIC, CURVED_FLAT")
                                              .validator(is_inflow_profile);
      Parameter<String> center = parameter(String("")).property("inflow_params/center").required();
      // Parameter<String> face_normal = parameter(String("")).property("inflow_params/face_normal");
      Parameter<String> flow_normal = parameter(String("")).property("inflow_params/flow_normal").required();
      Parameter<bool> projection_warning = parameter(false).property("inflow_params/projection_warning");
      Parameter<DT_> radius = parameter(DT_(0.0)).property("inflow_params/radius");
      Parameter<String> midpoint_a = parameter(String("")).property("inflow_params/midpoint_a");
      Parameter<String> midpoint_b = parameter(String("")).property("inflow_params/midpoint_b");
      Parameter<DT_> outer_radius = parameter(DT_(0.0)).property("inflow_params/outer_radius");
      Parameter<DT_> inner_radius = parameter(DT_(0.0)).property("inflow_params/inner_radius");
      Parameter<DT_> radius_major = parameter(DT_(0.0)).property("inflow_params/radius_major");
      Parameter<DT_> radius_minor = parameter(DT_(0.0)).property("inflow_params/radius_minor");
      Parameter<String> inflow_type = parameter(String(""))
                                              .property("inflow_params/type")
                                              .required()
                                              .name("Inflow geometry type")
                                              .help_string("Valid options are: CIRCLE, ELLIPSE, RING, RECTANGLE, CURVED_RECTANGLE")
                                              .validator(is_inflow_geometry)
                                              .needs_if(radius, is_circle)
                                              .needs_if(midpoint_a, is_rectangle)
                                              .needs_if(midpoint_b, is_rectangle)
                                              .needs_if(outer_radius, is_ring)
                                              .needs_if(inner_radius, is_ring)
                                              .needs_if(radius_major, is_ellipse)
                                              .needs_if(radius_minor, is_ellipse);
      // Parameter<String> test_point = parameter(String("")).property("inflow_params/test_point");

      InflowGeometry geometry_type() const
      {
        InflowGeometry tmp = InflowGeometry::unkown;
        const bool success = try_parse_inflow_geometry(*inflow_type, tmp);
        XASSERTM(success, "Could not parse inflow geometry");
        return tmp;
      }

      InflowType profile_type() const
      {
        InflowType tmp = InflowType::best;
        const bool success = try_parse_inflow_option(*inflow_profile, tmp);
        XASSERTM(success, "Could not parse inflow profile");
        return tmp;
      }
    };
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
    bool projection_warning;
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
        XASSERTM(vals[std::size_t(i)].parse(tmp[i]), "Could not parse entry");
      }
      return tmp;
    }
  public:

    InflowGeometryHandler(const Intern::InflowGeometryParser<DataType>& parser, const VecType& bb_min,
                          const VecType& bb_max, int id, InflowType inflow_type_ = InflowType::best,
                          InflowGeometry in_geo = InflowGeometry::unkown)
    : inflow_geometry(in_geo), inflow_type(inflow_type_), flow_id(id), projection_warning(*parser.projection_warning)
    {
      center = _parse_vector(*parser.center);
      flow_normal = _parse_vector(*parser.flow_normal);
      flow_normal.normalize();
      {
        // by which factor we shift our center into normal direction
        DataType linear_factor = Math::Limits<DataType>::max();
        // project the center onto the boundingbox
        for(const std::array<VecType, 3>& vecs : _generate_planes(bb_min, bb_max))
        {
          const VecType& lot_p = vecs[0];
          const VecType& span1 = vecs[1];
          const VecType& span2 = vecs[2];
          {
            VecType orto;
            Tiny::cross(orto, span1, span2);
            orto.normalize();
            if(Math::abs(Tiny::dot(orto, flow_normal)) <= Math::pow(Math::eps<DataType>(), DataType(0.8)))
              continue;
          }
          Tiny::Matrix<DataType, 3, 3> A{{-flow_normal[0], span1[0], span2[0]},
                                         {-flow_normal[1], span1[1], span2[1]},
                                         {-flow_normal[2], span1[2], span2[2]}};
          VecType rhs = center - lot_p;
          Tiny::Matrix<DataType, 3, 3> A_inv;
          A_inv.set_inverse(A);
          VecType sol = A_inv * rhs;
          linear_factor = Math::abs(sol[0]) < Math::abs(linear_factor) ? sol[0] : linear_factor;
        }
        if(this->projection_warning && Math::abs(linear_factor) >= DataType(0.01))
          std::printf("WARNING: Projection factor %f larger than 0.01\n", double(linear_factor));
        center = center + linear_factor * flow_normal;
      }
      d_plane = Tiny::dot(center, flow_normal);
    }

    virtual String get_bnd_name() const = 0;

    virtual MeshPartP create_mesh_part(const MeshType&, const DataType tolerance=1E-6) const = 0;

    virtual DirichletBoundaryFunc get_dirichlet_boundary_function(DataType throughput, const Material<DataType>& mat, const DataType) const = 0;

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
    explicit UnknownGeo(const Intern::InflowGeometryParser<DataType>& parser, const VecType& bb_min, const VecType& bb_max, int id) : BaseClass(parser, bb_min, bb_max, id) {}
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

    explicit Circle(const Intern::InflowGeometryParser<DataType>& parser, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::curved_flat, DataType rel = DataType(0.8))
    : BaseClass(parser, bb_min, bb_max, id, in_type, InflowGeometry::circle),
      radius(*parser.radius),
      inner_radius(radius * rel)
    {
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

    explicit Rectangle(const Intern::InflowGeometryParser<DataType>& parser, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::parabolic)
    : BaseClass(parser, bb_min, bb_max, id, in_type, InflowGeometry::rectangle),
      a(),
      b()
    {
      a = this->_parse_vector(*parser.midpoint_a);
      a = Intern::project_vector(a, this->flow_normal, this->d_plane);
      b = this->_parse_vector(*parser.midpoint_b);
      b = Intern::project_vector(b, this->flow_normal, this->d_plane);
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
        // parabolic function with x at sides 0, i.e. f(x) = (|a-c|^2-|(x-c)_|_(a-c)|^2)* (|b-c|^2-|(x-c)_|_(b-c)|^2)
        // additional, we have to normalize, .i.e. integrate above function across the rectangle
        case InflowType::parabolic:
        {
          const VecType aC = a - this->center;
          const VecType bC = b - this->center;
          const DataType euc_a = aC.norm_euclid_sqr();
          const DataType euc_b = bC.norm_euclid_sqr();
          const DataType scale = vol_flow * DataType(9)/ (DataType(16) * Math::cub(aC.norm_euclid()) * Math::cub(bC.norm_euclid()));
          return DirichletBoundaryFunc([scale, _center=this->center, _flow_normal=this->flow_normal, d_p=this->d_plane, aC, bC, euc_a, euc_b]
                  (const VecType& point) -> VecType
          {
            const VecType proj = Intern::project_vector(point, _flow_normal, d_p);
            const VecType dC = proj - _center;
            const DataType aux1 = Math::max(DataType(0), euc_a - Math::sqr(Tiny::dot(dC, aC))/euc_a);
            const DataType aux2 = Math::max(DataType(0), euc_b - Math::sqr(Tiny::dot(dC, bC))/euc_b);
            const DataType daux = aux1 * aux2;
            return (daux * scale) * _flow_normal ;
          });
        }
        default:
        {
          XABORTM("Not implemented");
        }
      }
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

    explicit CurvedRectangle(const Intern::InflowGeometryParser<DataType>& parser, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::parabolic)
    : BaseClass(parser, bb_min, bb_max, id, in_type),
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
          const DataType euc_b = bC.norm_euclid_sqr();
          DataType scale = (DataType(8) / DataType(3)) * aC.norm_euclid() * bC.norm_euclid()
                                  + DataType(2) * Math::pi<DataType>() * Math::sqr(radius);
          scale = vol_flow / scale;
          return DirichletBoundaryFunc([scale, _center=this->center, _flow_normal=this->flow_normal, d_p=this->d_plane, _a=this->a, radi=DataType(1)/this->radius, aC, bC, euc_b]
                  (const VecType& point) -> VecType
          {
            const VecType proj = Intern::project_vector(point, _flow_normal, d_p);
            const VecType dC = proj - _center;
            const DataType aux =  Math::max(DataType(0), DataType(1) - Math::sqr(Tiny::dot(dC, bC)/euc_b));
            const DataType aux2 = Math::max(DataType(0), DataType(1) - (proj-_a).norm_euclid_sqr()*Math::sqr(radi))
                                  + Math::max(DataType(0), DataType(1) - (proj-DataType(2)*_center+_a).norm_euclid_sqr() * Math::sqr(radi));
            const bool outside_rec = (Math::sqr(aC.norm_euclid_sqr()) - Math::sqr(Tiny::dot(dC, aC))) <= DataType(0);
            return ((outside_rec ? aux2 : aux) *scale) * _flow_normal;
          });
        }
        default:
        {
          XABORTM("Not implemented");
        }
      }
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

    explicit Ellipse(const Intern::InflowGeometryParser<DataType>& parser, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::parabolic)
    : BaseClass(parser, bb_min, bb_max, id, in_type, InflowGeometry::ellipse),
      radius_major(*parser.radius_major),
      radius_minor(*parser.radius_minor)
    {
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
        return (Math::abs(orth_dist) < tolerance) && bool(Math::sqr(dC[0]));
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

    explicit Ring(const Intern::InflowGeometryParser<DataType>& parser, const VecType& bb_min, const VecType& bb_max, int id, InflowType in_type = InflowType::parabolic)
    : BaseClass(parser, bb_min, bb_max, id, in_type, InflowGeometry::ring),
      radius_outer(*parser.outer_radius),
      radius_inner(*parser.inner_radius)
    {
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

  }; // class Rectangle

  template<typename MeshType_, typename DT_ = typename MeshType_::DataType>
  std::unique_ptr<InflowGeometryHandler<MeshType_, DT_>> create_geo_handler(const Intern::InflowGeometryParser<DT_>& parser, const Tiny::Vector<DT_, 3>& bb_min, const Tiny::Vector<DT_, 3>& bb_max, int id, InflowType in_type)
  {
    switch(parser.geometry_type())
    {
      case Gendie::InflowGeometry::circle:
      {
        return std::make_unique<Circle<MeshType_, DT_>>(parser, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::rectangle:
      {
        return std::make_unique<Rectangle<MeshType_, DT_>>(parser, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::curved_rectangle:
      {
        return std::make_unique<CurvedRectangle<MeshType_, DT_>>(parser, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::ring:
      {
        return std::make_unique<Ring<MeshType_, DT_>>(parser, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::ellipse:
      {
        XABORTM("ERROR: Ellipse definition not complete!");
        return std::make_unique<Ellipse<MeshType_, DT_>>(parser, bb_min, bb_max, id, in_type);
      }
      case Gendie::InflowGeometry::unkown:
      {
        XABORTM("Unkown inflow geometry type");
        break;
      }
      default:
      {
        XABORTM("Not Implemented " + stringify(parser.geometry_type()) + " yet");
        break;
      }
    }
    return std::make_unique<UnknownGeo<MeshType_, DT_>>(parser, bb_min, bb_max, id);
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

    void _init_from_parser(const Intern::InflowGeometryParser<DataType>& parser, const Tiny::Vector<DataType, 3>& bb_min, const Tiny::Vector<DataType, 3>& bb_max, int id)
    {
      _throughput = *parser.throughput;
      _material_id = *parser.material_id;
      _temperature = *parser.temperature;
      _inflow_geometry = parser.geometry_type();
      _inflow_type = parser.profile_type();
      _geo_handler = create_geo_handler<MeshType, DataType>(parser, bb_min, bb_max, id, _inflow_type);
    }

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
      Intern::InflowGeometryParser<DataType> parser;
      XASSERTM(inflow_property != nullptr, "Inflow property map is nullptr");
      XASSERTM(parser.parse(*inflow_property), "Could not parse inflow section");
      _init_from_parser(parser, bb_min, bb_max, id);
    }

    InflowBoundary(const Intern::InflowGeometryParser<DataType>& parser, const Tiny::Vector<DataType, 3>& bb_min, const Tiny::Vector<DataType, 3>& bb_max, int id)
      : _material_id(~Index(0)),
        _throughput(DataType(0)),
        _temperature(DataType(300)),
        _inflow_geometry(InflowGeometry::unkown),
        _inflow_type(InflowType::best)
    {
      _init_from_parser(parser, bb_min, bb_max, id);
    }

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

    DataType get_temperature() const
    {
      return _temperature;
    }

    DataType get_throughput() const
    {
      return _throughput;
    }
  };

  enum class InflowParserError : std::uint8_t
  {
    success = 0,
    no_inflow_config = 1,
    material_id_not_found = 2,
    invalid_inflow = 3
  };

  static inline const char* parse_error(InflowParserError t)
  {
    switch(t)
    {
      case Gendie::InflowParserError::success:
        return "Success";
      case Gendie::InflowParserError::no_inflow_config:
        return "No Inflow Config Found";
      case Gendie::InflowParserError::material_id_not_found:
        return "Material id not found in materials";
      case Gendie::InflowParserError::invalid_inflow:
        return "Failed to parse inflow entry";
    }
    return "Unknown";
  }

  template<typename InflowBound_, typename DT_, typename Domain_>
  Gendie::InflowParserError parse_inflow_boundaries(std::vector<InflowBound_>& inflow_boundaries, const FEAT::PropertyMap* inflow_list_config, const std::vector<Material<DT_>>& materials,
                              const Domain_& domain, String& inflow_info)
  {
    XASSERTM(inflow_list_config, "Inflow config not found");
    if(!inflow_list_config)
      return Gendie::InflowParserError::no_inflow_config;

    const auto bb_pair = FEAT::Control::Domain::DomainGlobalBoundingBoxHelper<Domain_>::calculate_global_boundingbox(domain, Index(0));
    // get boundingboxes from our voxel domain
    const auto& mesh_bb_min = bb_pair.first;
    const auto& mesh_bb_max = bb_pair.second;

    inflow_boundaries.clear();

    for(auto it = inflow_list_config->begin_section(); it != inflow_list_config->end_section(); ++it)
    {
      int id = -1;
      it->first.parse(id);
      Intern::InflowGeometryParser<DT_> parser;
      if(!parser.parse(*it->second))
        return InflowParserError::invalid_inflow;
      if(*parser.material_id >= materials.size())
        return InflowParserError::material_id_not_found;
      inflow_boundaries.emplace_back(parser, mesh_bb_min, mesh_bb_max, id);
      inflow_info+="Inflow " + inflow_boundaries[std::size_t(id)].get_bnd_name()  + " parameters:\n";
      inflow_info+=parser.display();
    }

    return InflowParserError::success;
  }


  template<typename Domain_, typename InflowBound_, int bnd_axis_ = 2>
  static inline void create_boundaries(Domain_& domain, const std::vector<InflowBound_>& inflow_boundaries,
                                       const char wall_name[] = "bnd:wall", const char out_name[] = "bnd:outflow")
  {
    typedef typename InflowBound_::DataType DataType_;
    typedef typename Domain_::LevelType::MeshType MeshType_;
    typedef typename Domain_::LevelType::PartType MeshPartType_;

    const auto bb_pair = FEAT::Control::Domain::DomainGlobalBoundingBoxHelper<Domain_>::calculate_global_boundingbox(domain, Index(0));
    // get boundingboxes from our voxel domain
    // const auto& mesh_bb_min = bb_pair.first;
    const auto& mesh_bb_max = bb_pair.second;

    const DataType_ tol_adp = FEAT::Math::pow(DataType_(2), -DataType_(domain.get_desired_level_max()));

    const DataType_ bnd_max = mesh_bb_max[bnd_axis_]- DataType_(0.001) * tol_adp;
    auto hit_r = [bnd_max](auto p) {return p[bnd_axis_] > bnd_max;};

    // create boundary meshparts on each level
    for(Index i(0); i < domain.size_physical(); ++i)
    {
      auto& mesh_node = *domain.at(i)->get_mesh_node();
      // if(param_holder.use_base_mesh)
      mesh_node.remove_all_mesh_parts();
      // create mesh-part from chart by using a hit-test factory
      std::vector<std::pair<FEAT::String, std::unique_ptr<MeshPartType_>>> bnd_inflows;
      bnd_inflows.reserve(inflow_boundaries.size());
      for(const auto& b : inflow_boundaries)
      {
        //TODO: tolerance adaptive to mesh size?
        bnd_inflows.push_back(std::make_pair(b.get_bnd_name(), b.create_mesh_part(domain.at(i)->get_mesh(), DataType_(1E-2) * tol_adp)));
      }
      Geometry::HitTestFactory<decltype(hit_r), MeshType_> factory_outflow(hit_r, domain.at(i)->get_mesh());

      // auto bnd_inflow = factory_inflow.make_unique();
      auto bnd_outflow = factory_outflow.make_unique();

      FEAT::Geometry::GlobalMaskedBoundaryFactory<MeshType_> boundary_factory(*mesh_node.get_mesh());
      // boundary_factory.add_mask_meshpart(*bnd_inflow);
      for(auto& [a,b] : bnd_inflows)
      {
        boundary_factory.add_mask_meshpart(*b);
      }

      boundary_factory.add_mask_meshpart(*bnd_outflow);
      for(const auto& v : mesh_node.get_halo_map())
        boundary_factory.add_halo(v.first,*v.second);
      boundary_factory.compile(domain.at(i).layer().comm());
      mesh_node.add_mesh_part(wall_name, boundary_factory.make_unique());
      for(auto& [a,b] : bnd_inflows)
      {
        mesh_node.add_mesh_part(a, std::move(b));
      }
      mesh_node.add_mesh_part(out_name, std::move(bnd_outflow));
    }
  }

}
