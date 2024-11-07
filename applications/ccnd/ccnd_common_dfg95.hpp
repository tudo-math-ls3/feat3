// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include "ccnd_common.hpp"

namespace CCND
{
  namespace DFG95
  {
    // template for the steady-state inflow function used in bench1 and bench2
    template<int dim>
    class SteadyInflowFunction;

    // 2D steady inflow function used in bench1 and bench2
    template<>
    class SteadyInflowFunction<2> :
      public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;

    public:
      explicit SteadyInflowFunction(DataType vmax) :
        _vmax(vmax)
      {
      }

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _d, _den;

      public:
        explicit Evaluator(const SteadyInflowFunction& function) :
          _vmax(function._vmax),
          _d(DataType(0.41)),
          _den(_d*_d)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType y = point[1];

          ValueType val;
          val[0] = (_vmax * DataType(4) * y * (_d - y)) / _den;
          val[1] = DataType(0);
          return val;
        }
      };
    }; // class SteadyInflowFunction<2>

    // 3D steady inflow function used in bench1 and bench2
    template<>
    class SteadyInflowFunction<3> :
      public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;

    public:
      explicit SteadyInflowFunction(DataType vmax) :
        _vmax(vmax)
      {
      }

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _d, _den;

      public:
        explicit Evaluator(const SteadyInflowFunction& function) :
          _vmax(function._vmax),
          _d(DataType(0.41)),
          _den(_d*_d*_d*_d)
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType y = point[1];
          const DataType z = point[2];

          ValueType val;
          val[0] = (_vmax * DataType(16) * y * (_d - y) * z * (_d - z)) / _den;
          val[1] = val[2] = DataType(0);
          return val;
        }
      };
    }; // class SteadyInflowFunction

    // template for the unsteady inflow function used in bench3
    template<int dim>
    class UnsteadyInflowFunction;

    // 2D steady inflow function used in bench3
    template<>
    class UnsteadyInflowFunction<2> :
      public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 2;
      typedef Analytic::Image::Vector<2> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _time;

    public:
      explicit UnsteadyInflowFunction(DataType vmax, DataType t) :
        _vmax(vmax),
        _time(t)
      {
      }

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _t, _pi, _d, _scale;

      public:
        explicit Evaluator(const UnsteadyInflowFunction& function) :
          _vmax(function._vmax),
          _t(function._time),
          _pi(Math::pi<DataType>()),
          _d(DataType(0.41)),
          _scale(Math::sin(_pi*_t/DataType(8))/(_d*_d))
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType y = point[1];

          ValueType val;
          val[0] = (_vmax * DataType(4) * y * (_d - y)) * _scale;
          val[1] = DataType(0);
          return val;
        }
      };
    }; // class UnsteadyInflowFunction<2>

    // 3D steady inflow function used in bench3
    template<>
    class UnsteadyInflowFunction<3> :
      public Analytic::Function
    {
    public:
      static constexpr int domain_dim = 3;
      typedef Analytic::Image::Vector<3> ImageType;
      static constexpr bool can_value = true;

    protected:
      DataType _vmax;
      DataType _time;

    public:
      explicit UnsteadyInflowFunction(DataType vmax, DataType t) :
        _vmax(vmax),
        _time(t)
      {
      }

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax, _t, _pi, _d, _scale;

      public:
        explicit Evaluator(const UnsteadyInflowFunction& function) :
          _vmax(function._vmax),
          _t(function._time),
          _pi(Math::pi<DataType>()),
          _d(DataType(0.41)),
          _scale(Math::sin(_pi*_t/DataType(8))/(_d*_d*_d*_d))
        {
        }

        ValueType value(const PointType& point)
        {
          const DataType y = point[1];
          const DataType z = point[2];

          ValueType val;
          val[0] = (_vmax * DataType(16) * y * (_d - y) * z * (_d - z)) * _scale;
          val[1] = val[2] = DataType(0);
          return val;
        }
      };
    }; // class UnsteadyInflowFunction


    // accumulator for benchmark body forces (i.e. drag and lift)
    // this class is used for the computation of the 'line integration' variants of drag and lift
    // for details see Giles et al: "Adaptive error control for finite element approximations of the lift and drag coefficients in viscous flow"
    class BenchBodyForceAccumulator
    {
    public:
      Tiny::Matrix<DataType, 3, 2> raw;

      BenchBodyForceAccumulator()
      {
        raw.format();
      }

      /// 2D variant
      template<typename T_>
      void operator()(
        const T_ omega,
        const Tiny::Vector<T_, 2, 2>& /*pt*/,
        const Tiny::Matrix<T_, 2, 1, 2, 1>& jac,
        const Tiny::Vector<T_, 2, 2>& /*val_v*/,
        const Tiny::Matrix<T_, 2, 2, 2, 2>& grad_v,
        const T_ val_p)
      {
#if 0
        // Note: This is the old FeatFlow 1.x / FeatFlow 2 way of computing the 2D body forces

        // compute normal and tangential
        const T_ n2 = T_(1) / Math::sqrt(jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0));
        const T_ tx = jac(0,0) * n2;
        const T_ ty = jac(1,0) * n2;
        const T_ nx = -ty;
        const T_ ny =  tx;

        Tiny::Matrix<T_, 2, 2, 2, 2> nt;
        nt(0,0) = tx * nx;
        nt(0,1) = tx * ny;
        nt(1,0) = ty * nx;
        nt(1,1) = ty * ny;

        const T_ dut = Tiny::dot(nt, grad_v);

        raw(0, 0) += omega * dut * ny;
        raw(0, 1) -= omega * val_p * nx;
        raw(1, 0) -= omega * dut * nx;
        raw(1, 1) -= omega * val_p * ny;
#endif
        // compute normal vector
        const Tiny::Vector<T_, 2> n = Tiny::orthogonal(jac).normalize();
        raw(0, 0) -= omega * (T_(2) * grad_v(0,0) * n[0] + (grad_v(0, 1) + grad_v(1, 0)) * n[1]);
        raw(1, 0) -= omega * (T_(2) * grad_v(1,1) * n[1] + (grad_v(1, 0) + grad_v(0, 1)) * n[0]);
        raw(0, 1) -= omega * val_p * n[0];
        raw(1, 1) -= omega * val_p * n[1];
      }

      /// 3D variant
      template<typename T_>
      void operator()(
        const T_ omega,
        const Tiny::Vector<T_, 3, 3>& /*pt*/,
        const Tiny::Matrix<T_, 3, 2, 3, 2>& jac,
        const Tiny::Vector<T_, 3, 3>& /*val_v*/,
        const Tiny::Matrix<T_, 3, 3, 3, 3>& grad_v,
        const T_ val_p)
      {
#if 0
        // Note: This is the old FeatFlow 1.x / FeatFlow 2 way of computing the 3D body forces

        // compute normal and tangential
        const T_ n2 = T_(1) / Math::sqrt(jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0));
        const T_ tx = jac(0,0) * n2;
        const T_ ty = jac(1,0) * n2;
        const T_ nx = -ty;
        const T_ ny =  tx;

        Tiny::Matrix<T_, 3, 3, 3, 3> nt;
        nt.format();
        nt(0,0) = tx * nx;
        nt(0,1) = tx * ny;
        nt(1,0) = ty * nx;
        nt(1,1) = ty * ny;

        const T_ dut = Tiny::dot(nt, grad_v);
        //const T_ dpf1 = _nu;
        //const T_ dpf2 = (2.0 / (0.1*Math::sqr(_v_max*(4.0/9.0))* 0.41)); // = 2 / (rho * U^2 * D * H)

        raw(0, 0) += omega * dut * ny;
        raw(0, 1) -= omega * val_p * nx;
        raw(1, 0) -= omega * dut * nx;
        raw(1, 1) -= omega * val_p * ny;
#endif

        // compute normal vector
        const Tiny::Vector<T_, 3> n = Tiny::orthogonal(jac).normalize();
        raw(0, 0) -= omega * (T_(2) * grad_v(0,0) * n[0] + (grad_v(0, 1) + grad_v(1, 0)) * n[1] + (grad_v(0, 2) + grad_v(2, 0)) * n[2]);
        raw(1, 0) -= omega * (T_(2) * grad_v(1,1) * n[1] + (grad_v(1, 2) + grad_v(2, 1)) * n[2] + (grad_v(1, 0) + grad_v(0, 1)) * n[0]);
        raw(2, 0) -= omega * (T_(2) * grad_v(2,2) * n[2] + (grad_v(2, 0) + grad_v(0, 2)) * n[0] + (grad_v(2, 1) + grad_v(1, 2)) * n[1]);
        raw(0, 1) -= omega * val_p * n[0];
        raw(1, 1) -= omega * val_p * n[1];
        raw(2, 1) -= omega * val_p * n[2];
      }

      DataType drag(DataType nu) const
      {
        return nu * raw(0, 0) - raw(0, 1);
      }

      DataType lift(DataType nu) const
      {
        return nu * raw(1, 0) - raw(1, 1);
      }

      DataType side(DataType nu) const
      {
        return nu * raw(2, 0) - raw(2, 1);
      }

      void sync(const Dist::Comm& comm)
      {
        comm.allreduce(&raw.v[0][0], &raw.v[0][0], std::size_t(2*dim), Dist::op_sum);
      }
    }; // class BenchBodyForces<...,2>


    // class for collecting benchmark post-processing summary data
    class BenchmarkAnalysis
    {
    public:
      /// chart for circle obstacle
      const ChartBaseType* chart_circle = nullptr;
      /// chart for sphere obstacle
      const ChartBaseType* chart_sphere = nullptr;
      /// chart for outer pipe of bench 7
      const ChartBaseType* chart_pipe = nullptr;

      /// meshpart for FBM boundary
      const MeshPartType* mesh_part_fbm = nullptr;
      /// meshpart for circle boundary
      const MeshPartType* mesh_part_bnd_c = nullptr;
      /// meshpart for sphere boundary
      const MeshPartType* mesh_part_bnd_s = nullptr;
      /// meshpart for left edge
      const MeshPartType* mesh_part_bnd_l = nullptr;
      /// meshpart for right edge
      const MeshPartType* mesh_part_bnd_r = nullptr;
      /// meshpart for inflow
      const MeshPartType* mesh_part_bnd_in = nullptr;
      /// meshpart for outflow
      const MeshPartType* mesh_part_bnd_out = nullptr;
      /// meshpart for upper inner edge
      const MeshPartType* mesh_part_inner_u = nullptr;
      /// meshpart for lower inner edge
      const MeshPartType* mesh_part_inner_l = nullptr;

      /// trace assemblers for obstacle body force and fluxes
      std::unique_ptr<Assembly::TraceAssembler<TrafoType>> body_force_asm, flux_asm_u, flux_asm_l, flux_asm_in, flux_asm_out;

      /// inverse mapping data for pressure evaluation
      Trafo::InverseMappingData<DataType, dim> point_iv_a, point_iv_e;

      /// characteristic vector of obstacle body
      LAFEM::DenseVector<DataType, IndexType> vec_char_obstacle;

      /// cubature rule for surface integration of body forces
      String cubature_body_forces_surf = "gauss-legendre:3";

      /// cubature rule for fluxes
      String cubature_flux = "gauss-legendre:3";

      /// cubature rule for general post processing
      String cubature_postproc = "gauss-legendre:3";

      // post-processing data

      // reference values for benchmark values
      DataType ref_drag = DataType(dim == 2 ? 5.57953523384 : 6.18533);
      DataType ref_lift = DataType(dim == 2 ? 0.010618948146 : 0.009401);
      DataType ref_side = DataType(0.0);
      DataType ref_p_diff = DataType(dim == 2 ? 0.11752016697 : 0.170826996);
      DataType ref_p_drop = DataType(dim == 2 ? 0.06943442 : 0.126759);

      // drag/lift forces by surface integration
      DataType drag_coeff_surf, lift_coeff_surf, side_coeff_surf, drag_err_surf, lift_err_surf, side_err_surf;

      // drag/lift forces by volume integration
      DataType drag_coeff_vol, lift_coeff_vol, side_coeff_vol, drag_err_vol, lift_err_vol, side_err_vol;

      // pressure values
      DataType pres_diff, pres_diff_err;

      // pressure drop from inflow to outflow
      DataType pres_drop, pres_drop_err;

      // flow through upper/lower region and outflow flux
      DataType flux_upper, flux_lower, flux_out;

      // velocity field information
      Assembly::VelocityInfo<DataType, dim> velo_info;

      static constexpr std::size_t padlen = std::size_t(30);

      BenchmarkAnalysis() :
        drag_coeff_surf(0.0), lift_coeff_surf(0.0), side_coeff_surf(0.0), drag_err_surf(0.0), lift_err_surf(0.0), side_err_surf(0.0),
        drag_coeff_vol(0.0), lift_coeff_vol(0.0), side_coeff_vol(0.0), drag_err_vol(0.0), lift_err_vol(0.0), side_err_vol(0.0),
        pres_diff(0.0), pres_diff_err(0.0),
        pres_drop(0.0), pres_drop_err(0.0),
        flux_upper(0.0), flux_lower(0.0), flux_out(0.0)
      {
      }

      virtual ~BenchmarkAnalysis()
      {
      }

      virtual String format(bool print_errors = false, int prec = fp_num_digs+5) const
      {
        String s;
        const char pc = '.';
        // append coefficients and velocity info
        s += "Solution Analysis:";
        s += String("\nDrag Coefficient (Surface)").pad_back(padlen, pc) + ".: " + stringify_fp_fix(drag_coeff_surf, prec);
        if(print_errors)
          s += "   [ Error: " + stringify_fp_sci(drag_err_surf, prec) + " ]";
        s += String("\nLift Coefficient (Surface)").pad_back(padlen, pc) + ".: " + stringify_fp_fix(lift_coeff_surf, prec);
        if(print_errors)
          s += "   [ Error: " + stringify_fp_sci(lift_err_surf, prec) + " ]";
        s += String("\nSide Coefficient (Surface)").pad_back(padlen, pc) + ".: " + stringify_fp_fix(side_coeff_surf, prec);
        if(print_errors)
          s += "   [ Error: " + stringify_fp_sci(side_err_surf, prec) + " ]";
        s += String("\nDrag Coefficient (Volume)").pad_back(padlen, pc) + ".: " + stringify_fp_fix(drag_coeff_vol, prec);
        if(print_errors)
          s += "   [ Error: " + stringify_fp_sci(drag_err_vol, prec) + " ]";
        s += String("\nLift Coefficient (Volume)").pad_back(padlen, pc) + ".: " + stringify_fp_fix(lift_coeff_vol, prec);
        if(print_errors)
          s += "   [ Error: " + stringify_fp_sci(lift_err_vol, prec) + " ]";
        s += String("\nSide Coefficient (Volume)").pad_back(padlen, pc) + ".: " + stringify_fp_fix(side_coeff_vol, prec);
        if(print_errors)
          s += "   [ Error: " + stringify_fp_sci(side_err_vol, prec) + " ]";
        s += String("\nPressure Difference").pad_back(padlen, pc) + ".: " + stringify_fp_fix(pres_diff, prec);
        if(print_errors)
          s += "   [ Error: " + stringify_fp_sci(pres_diff_err, prec) + " ]";
        s += String("\nPressure Drop").pad_back(padlen, pc) + ".: " + stringify_fp_fix(pres_drop, prec);
        if(print_errors)
          s += "   [ Error: " + stringify_fp_sci(pres_drop_err, prec) + " ]";
        s += String("\nUpper Flux").pad_back(padlen, pc) + ".: " + stringify_fp_fix(flux_upper, prec);
        s += String("\nLower Flux").pad_back(padlen, pc) + ".: " + stringify_fp_fix(flux_lower, prec);
        s += String("\nOut Flux").pad_back(padlen, pc) + ".: " + stringify_fp_fix(flux_out, prec);
        s += String("\n") + velo_info.format_string(prec, padlen, pc);

        return s;
      }

      virtual String format_compact(String prefix, int prec = fp_num_digs) const
      {
        String s;
        //s += prefix + "DCL/LCL/DCV/LCV/SCV/PD: ";
        s += prefix + "Surface Forces (D/L/S): ";
        s += stringify_fp_fix(drag_coeff_surf, prec) + " " + stringify_fp_fix(lift_coeff_surf, prec) + " " + stringify_fp_fix(side_coeff_surf, prec) + " \n";
        s += prefix + "Volume Forces (D/L/S).: ";
        s += stringify_fp_fix(drag_coeff_vol, prec) + " " + stringify_fp_fix(lift_coeff_vol, prec) + " " + stringify_fp_fix(side_coeff_vol, prec) + "\n";
        s += prefix + "Pressure Difference...: ";
        s += stringify_fp_fix(pres_diff, prec) + "\n";
        s += prefix + "Pressure Drop.........: ";
        s += stringify_fp_fix(pres_drop, prec) + "\n";
        s += prefix + "Fluxes Upper/Lower/Out: ";
        s += stringify_fp_fix(flux_upper, prec) + " " + stringify_fp_fix(flux_lower, prec) + " " + stringify_fp_fix(flux_out, prec) + "\n";
        s += prefix + "Velocity H0/H1/Vor/Div: ";
        s += stringify_fp_fix(velo_info.norm_h0, prec) + " " + stringify_fp_fix(velo_info.norm_h1, prec) + " ";
        s += stringify_fp_fix(velo_info.vorticity, prec) + " " + stringify_fp_fix(velo_info.divergence, prec);
        return s;
      }

      virtual void create(const MeshAtlasType& atlas, const MeshNodeType& mesh_node, const SpaceVeloType& space_v)
      {
        // get charts
        chart_circle = atlas.find_mesh_chart("circle");
        chart_sphere = atlas.find_mesh_chart("sphere");
        chart_pipe   = atlas.find_mesh_chart("pipe");   // bench7 only

        // get meshparts on finest level
        mesh_part_fbm     = mesh_node.find_mesh_part("fbm");
        mesh_part_bnd_c   = mesh_node.find_mesh_part("bnd:c");
        mesh_part_bnd_s   = mesh_node.find_mesh_part("bnd:sphere");
        mesh_part_bnd_l   = mesh_node.find_mesh_part("bnd:l");
        mesh_part_bnd_r   = mesh_node.find_mesh_part("bnd:r");
        mesh_part_bnd_in  = mesh_node.find_mesh_part("bnd:in");
        mesh_part_bnd_out = mesh_node.find_mesh_part("bnd:out");
        mesh_part_inner_u = mesh_node.find_mesh_part("inner:u");
        mesh_part_inner_l = mesh_node.find_mesh_part("inner:l");

        // create trace assembler for body force assembly
        body_force_asm.reset(new Assembly::TraceAssembler<TrafoType>(space_v.get_trafo()));
        if(mesh_part_fbm != nullptr)
          body_force_asm->add_mesh_part(*mesh_part_fbm);
        if(mesh_part_bnd_c != nullptr)
          body_force_asm->add_mesh_part(*mesh_part_bnd_c);
        if(mesh_part_bnd_s != nullptr)
          body_force_asm->add_mesh_part(*mesh_part_bnd_s);
        body_force_asm->compile();

        // create trace assembler for upper flux
        flux_asm_u.reset(new Assembly::TraceAssembler<TrafoType>(space_v.get_trafo()));
        if(mesh_part_inner_u != nullptr)
          flux_asm_u->add_mesh_part(*mesh_part_inner_u);
        flux_asm_u->compile();

        // create trace assembler for lower flux
        flux_asm_l.reset(new Assembly::TraceAssembler<TrafoType>(space_v.get_trafo()));
        if(mesh_part_inner_l != nullptr)
          flux_asm_l->add_mesh_part(*mesh_part_inner_l);
        flux_asm_l->compile();

        // create trace assembler for in flux
        flux_asm_in.reset(new Assembly::TraceAssembler<TrafoType>(space_v.get_trafo()));
        if(mesh_part_bnd_l != nullptr)
          flux_asm_in->add_mesh_part(*mesh_part_bnd_l);
        if(mesh_part_bnd_in != nullptr)
          flux_asm_in->add_mesh_part(*mesh_part_bnd_in);
        flux_asm_in->compile();

        // create trace assembler for out flux
        flux_asm_out.reset(new Assembly::TraceAssembler<TrafoType>(space_v.get_trafo()));
        if(mesh_part_bnd_r != nullptr)
          flux_asm_out->add_mesh_part(*mesh_part_bnd_r);
        if(mesh_part_bnd_out != nullptr)
          flux_asm_out->add_mesh_part(*mesh_part_bnd_out);
        flux_asm_out->compile();

        // unmap pressure evaluation points p_a and p_e
        typedef Trafo::InverseMapping<TrafoType, DataType> InvMappingType;
        InvMappingType inv_mapping(space_v.get_trafo());

        // reference pressure points
        typename InvMappingType::ImagePointType v_a, v_e;
        if(dim == 2)
        {
          v_a[0] = DataType(0.15);
          v_e[0] = DataType(0.25);
          v_a[1] = v_e[1] = DataType(0.2);
        }
        else
        {
          v_a[0] = DataType(0.45);
          v_e[0] = DataType(0.55);
          v_a[1] = v_e[1] = DataType(0.2);
          v_a[2] = v_e[2] = DataType(0.205);
        }

        // unmap points
        point_iv_a = inv_mapping.unmap_point(v_a, true);
        point_iv_e = inv_mapping.unmap_point(v_e, true);

        // create characteristic function vector for circle boundary
        // this is needed for the volumetric drag/lift computation
        vec_char_obstacle = LAFEM::DenseVector<DataType, IndexType>(space_v.get_num_dofs(), DataType(0));
        if(mesh_part_bnd_c != nullptr)
        {
          LAFEM::UnitFilter<DataType, IndexType> filter_char;
          Assembly::UnitFilterAssembler<MeshType> unit_asm;
          unit_asm.add_mesh_part(*mesh_part_bnd_c);
          unit_asm.assemble(filter_char, space_v);
          filter_char.get_filter_vector().format(DataType(1));
          filter_char.filter_sol(vec_char_obstacle);
        }
        if(mesh_part_bnd_s != nullptr)
        {
          LAFEM::UnitFilter<DataType, IndexType> filter_char;
          Assembly::UnitFilterAssembler<MeshType> unit_asm;
          unit_asm.add_mesh_part(*mesh_part_bnd_s);
          unit_asm.assemble(filter_char, space_v);
          filter_char.get_filter_vector().format(DataType(1));
          filter_char.filter_sol(vec_char_obstacle);
        }
        if(mesh_part_fbm != nullptr)
        {
          LAFEM::UnitFilter<DataType, IndexType> filter_char;
          Assembly::UnitFilterAssembler<MeshType> unit_asm;
          unit_asm.add_mesh_part(*mesh_part_fbm);
          unit_asm.assemble(filter_char, space_v);
          filter_char.get_filter_vector().format(DataType(1));
          filter_char.filter_sol(vec_char_obstacle);
        }
      }

      virtual DataType bench_body_factor(int bench, DataType v_max)
      {
        if(bench == 7)
          return DataType(2) / (DataType(0.01)*Math::sqr(v_max*(DataType(0.5)))); // = 2 / (rho * U^2 * D^2)
        else if(dim == 2)
          return DataType(2) / (DataType(0.100)*Math::sqr(v_max*(DataType(2)/DataType(3)))); // = 2 / (rho * U^2 * D)
        else if(dim == 3)
          return DataType(2) / (DataType(0.041)*Math::sqr(v_max*(DataType(4)/DataType(9)))); // = 2 / (rho * U^2 * D * H)
        else
          return DataType(1);
      }

      virtual void compute_body_forces_surf(const Dist::Comm& comm, const LocalVeloVector& vec_sol_v, const LocalPresVector& vec_sol_p,
        const SpaceVeloType& space_velo, const SpacePresType& space_pres, DataType nu, DataType v_max, int bench)
      {
        Cubature::DynamicFactory cubature_factory(cubature_body_forces_surf);
        BenchBodyForceAccumulator body_force_accum;
        body_force_asm->assemble_flow_accum(body_force_accum, vec_sol_v, vec_sol_p, space_velo, space_pres, cubature_factory);
        body_force_accum.sync(comm);

        DataType bbf = bench_body_factor(bench, v_max);
        drag_coeff_surf = bbf * body_force_accum.drag(nu);
        lift_coeff_surf = bbf * body_force_accum.lift(nu);
        side_coeff_surf = bbf * body_force_accum.side(nu);
        if(Math::abs(ref_drag) > DataType(1E-7))
          drag_err_surf = Math::abs((drag_coeff_surf - ref_drag) / ref_drag);
        if(Math::abs(ref_lift) > DataType(1E-7))
          lift_err_surf = Math::abs((lift_coeff_surf - ref_lift) / ref_lift);
        if(Math::abs(ref_side) > DataType(1E-7))
          side_err_surf = Math::abs((side_coeff_surf - ref_side) / ref_side);
      }

      virtual void compute_body_forces_vol(const Dist::Comm& comm, const LocalVeloVector& vec_def_unsynced, DataType DOXY(nu), DataType v_max, int bench)
      {
        Tiny::Vector<DataType, 3> forces(DataType(0));
        Tiny::Vector<DataType, dim, 3>& frc = forces.template size_cast<dim>();

        XASSERT(vec_def_unsynced.size() == vec_char_obstacle.size());

        // get the vector arrays
        const Index n = vec_def_unsynced.size();
        const auto* vdef = vec_def_unsynced.elements();
        const auto* vchr = vec_char_obstacle.elements();
        for(Index i(0); i < n; ++i)
        {
          frc.axpy(vchr[i], vdef[i]);
        }

        // sum up over all processes
        comm.allreduce(forces.v, forces.v, 3u, Dist::op_sum);

        // save the values
        DataType bbf = bench_body_factor(bench, v_max);
        drag_coeff_vol = bbf * forces[0];
        lift_coeff_vol = bbf * forces[1];
        side_coeff_vol = bbf * forces[2];
        if(Math::abs(ref_drag) > 1E-7)
          drag_err_vol = Math::abs((drag_coeff_vol - ref_drag) / ref_drag);
        if(Math::abs(ref_lift) > 1E-7)
          lift_err_vol = Math::abs((lift_coeff_vol - ref_lift) / ref_lift);
        if(Math::abs(ref_side) > 1E-7)
          side_err_vol = Math::abs((side_coeff_vol - ref_side) / ref_side);
      }

      virtual void compute_pressure_difference(const Dist::Comm& comm, const LocalPresVector& vec_sol_p, const SpacePresType& space_pres)
      {
        // evaluate pressure
        auto pval_a = Assembly::DiscreteEvaluator::eval_fe_function(point_iv_a, vec_sol_p, space_pres);
        auto pval_e = Assembly::DiscreteEvaluator::eval_fe_function(point_iv_e, vec_sol_p, space_pres);

        // compute pressure mean
        const DataType p_a = pval_a.mean_value_dist(comm);
        const DataType p_e = pval_e.mean_value_dist(comm);
        const DataType d_p = p_a - p_e;

        // compute error to reference values
        pres_diff = d_p;
        if(Math::abs(ref_p_diff) > 1E-7)
          pres_diff_err = Math::abs((d_p - ref_p_diff) / ref_p_diff);
      }

      virtual void compute_pressure_drop(const Dist::Comm& comm, const LocalPresVector& vec_sol_p, const SpacePresType& space_pres, int bench)
      {
        Cubature::DynamicFactory cubature_factory(cubature_flux);

        DataType pv[2] =
        {
          flux_asm_in->assemble_discrete_integral(vec_sol_p, space_pres, cubature_factory),
          flux_asm_out->assemble_discrete_integral(vec_sol_p, space_pres, cubature_factory)
        };

        comm.allreduce(pv, pv, 2u, Dist::op_sum);

        // compute error to reference values
        if(bench == 7)
          pres_drop = (pv[0] - pv[1]) / (DataType(0.205*0.205) * Math::pi<DataType>());
        else if(dim == 2)
          pres_drop = (pv[0] - pv[1]) / DataType(0.41);
        else if(dim == 3)
          pres_drop = (pv[0] - pv[1]) / DataType(0.41*0.41);
        if(Math::abs(ref_p_drop) > 1E-7)
          pres_drop_err = Math::abs((pres_drop - ref_p_drop) / ref_p_drop);
      }

      virtual void compute_fluxes(const Dist::Comm& comm, const LocalVeloVector& vec_sol_v, const SpaceVeloType& space_velo)
      {
        Cubature::DynamicFactory cubature_factory(cubature_flux);

        DataType fx[3] =
        {
          flux_asm_u->assemble_discrete_integral(vec_sol_v, space_velo, cubature_factory)[0],
          flux_asm_l->assemble_discrete_integral(vec_sol_v, space_velo, cubature_factory)[0],
          flux_asm_out->assemble_discrete_integral(vec_sol_v, space_velo, cubature_factory)[0]
        };

        comm.allreduce(fx, fx, 3u, Dist::op_sum);

        flux_upper = fx[0];
        flux_lower = fx[1];
        flux_out = fx[2];
      }

      virtual void compute_velocity_info(const Dist::Comm& comm, const LocalVeloVector& vec_sol_v, const SpaceVeloType& space_velo)
      {
        velo_info = Assembly::VelocityAnalyser::compute(vec_sol_v, space_velo, cubature_postproc);
        velo_info.synchronize(comm);
      }
    }; // BenchmarkAnalysis
  } // namespace DFG95
} // namespace CCND
